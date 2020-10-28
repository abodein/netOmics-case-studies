Zscore <- function(x){
    (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

FClog <- function(x){
    x <- x + abs(min(x))
    max(x) - min(x)
}


remove.fold.change <- function(X, cutoff, hist.title = "", plot = TRUE){
    X <- as.data.frame(X)
    fc.func <- function(x){
        B <- max(x)
        A <- min(x)
        return(B-A)
    }
    fc <- lapply(X, function(x) fc.func(x)) %>% unlist
    if(plot){
        hist(fc, breaks = 20, main = hist.title)
        abline(v=log(2), col = "blue")
        abline(v=log(3), col = "red")
        legend("topright", lty = c(1,1), col = c("blue", "red"), legend = c("log(2)", "log(3)"))
    }
    X.filtered <- X[,which(fc>cutoff)]
    return(X.filtered)
}

wrapper.lmms <- function(X, time, sampleID){
    RESULT = list()
    tryCatch(
        expr = {
            RESULT[["pspline"]] <- lmms::lmmSpline(data = X, time = time,
                                                   sampleID = sampleID, deri = FALSE,
                                                   basis = "p-spline", timePredict = 1:3,
                                                   numCores = 4,
                                                   keepModels = TRUE)
        },
        error = function(e){ 
            # (Optional)
            # Do this if an error is caught...
        }
    )
    
    
    tryCatch(
        expr = {
            RESULT[["cubic"]] <- lmms::lmmSpline(data = X, time = time,
                                                 sampleID = sampleID, deri = FALSE,
                                                 basis = "cubic", timePredict = 1:3,
                                                 numCores = 4,
                                                 keepModels = TRUE)
        },
        error = function(e){ 
            # (Optional)
            # Do this if an error is caught...
        }
    )
    
    tryCatch(
        expr = {
            RESULT[["cubic.pspline"]] <- lmms::lmmSpline(data = X, time = time,
                                                         sampleID = sampleID, deri = FALSE,
                                                         basis = "cubic p-spline", timePredict = 1:3,
                                                         numCores = 4,
                                                         keepModels = TRUE)
        },
        error = function(e){ 
            # (Optional)
            # Do this if an error is caught...
        }
    )
    
    
    RESULT[["summary"]] <- imap(RESULT, ~{.x@modelsUsed %>%
            table() %>%
            as.data.frame() %>% 
            set_names(c("Model", .y))}) %>%
        plyr::join_all() %>% column_to_rownames("Model") %>%
        t
    RESULT[["summary"]][is.na(RESULT[["summary"]])] <- 0 
    RESULT[["X"]] <- X
    
    return(RESULT)
}

summary_cluster <- function(cluster.info){
    cluster.info %>% group_by(cluster, block) %>% summarise(N = n()) %>%
        spread(block, N) %>% column_to_rownames("cluster")
}

cluster_homo_block <- function(cluster.info){
    d <- cluster.info %>% group_by(cluster, block) %>%
        summarise(N = n())
    d.sum <- d %>% group_by(block) %>% summarise(total = sum(N))
    left_join(d,d.sum) %>%
        mutate(dist =  sqrt((N - (total/length(unique(d$cluster))))^2)) %>%
        pull(dist) %>% mean
}
get_stats_graph <- function(graph){
    node.c <- sum(degree(graph) != 0)
    node.i <- sum(degree(graph) == 0)
    edge <- ecount(graph)
    density <- igraph::edge_density(graph = graph)
    return(list(node.c = node.c, node.i = node.i, edge = edge, density = density))
}

quick_silhouette <- function(cb.data, cluster.info){
    # timeOmics
    D <- cor(scale(cb.data), method = "spearman")
    D <- (-D+1)^2
    cluster = as.data.frame(list("molecule" = rownames(D))) %>%
        left_join(cluster.info) %>% pull(cluster) %>%
        as.factor() %>% as.numeric()
    
    sil <- cluster::silhouette(x = cluster, dmatrix = D) 
    summary(sil)
}
remove.isolated.nodes <- function(graph, simplify = FALSE){
    stopifnot(is(graph, "igraph"))
    if(simplify == TRUE){
        graph <- igraph::simplify(graph)
    }
    Isolated <-  which(degree(graph) == 0)
    graph.2 = delete.vertices(graph, Isolated)
    return(graph.2)
}

get_stats_graph_2 <- function(graph){
    va <- vertex_attr(graph) %>% as.data.frame()
    ea <- as_long_data_frame(graph) %>% dplyr::select(from_type, to_type)
    ea[is.na(ea)] <- "0"
    levels.ea <- factor(c(ea$from_type, ea$to_type)) %>% levels
    # ea <- ea %>% mutate(from_type = factor(from_type, levels = levels.ea),
    #                     to_type = factor(to_type, levels = levels.ea)) %>%
    #     mutate(from = ifelse(as.numeric(from_type) < as.numeric(to_type), from_type, to_type),
    #            to = ifelse(as.numeric(to_type) < as.numeric(from_type), to_type, from_type)) %>%
    #     dplyr::select(from, to)
    table(dplyr::select(ea, c(from_type, to_type)))
    
}

get_go_info <- function(go){
    go.info <- GOTERM[[go]]
    term <- ontology <- ""
    if("Term" %in% slotNames(go.info)){
        term = slot(go.info, "Term")
    }
    if("Ontology" %in% slotNames(go.info)){
        ontology = slot(go.info, "Ontology")
    }
    return(list(name = go,
                lab = go,
                description = paste0(term, " (", ontology, ")")) %>% 
               as.data.frame())
}

add_go_term_layer <- function(graph){
    go_enrich_graph_tmp <- signif_goterm_map.all %>%
        na.omit() %>%
        filter(input %in% names(V(graph))) %>%
        igraph::graph.data.frame(directed = FALSE) 
    graph_go <- igraph::union(graph, go_enrich_graph_tmp)
    vertex_attr(graph_go)$type <- ifelse(str_detect(vertex_attr(graph_go)$name, "GO:"), "GO", vertex_attr(graph_go)$type)
    return(graph_go)
}

add_go_term_annot <- function(graph, go.graph){
    interaction.graph <- goterm_map %>%
        dplyr::select(input, target) %>%
        na.omit() %>%
        filter(input %in% names(V(graph)), 
               target %in% V(go.graph)$name) %>%
        igraph::graph.data.frame(directed = FALSE) 
    
    graph_tmp <- igraph::union(graph, interaction.graph)
    full_graph <- igraph::union(graph_tmp, go.graph)
    vertex_attr(full_graph) <- vertex_attr(full_graph) %>% as.data.frame(stringsAsFactors = FALSE) %>% 
        mutate(type = ifelse(!is.na(type_1), type_1, type_2)) %>%
        dplyr::select(-c(type_1,type_2)) %>% 
        as.list()
    return(full_graph)
}
add_metabo_layer <- function(graph){
    # graph
    interaction.tmp <- interaction_met_prot.df %>% filter(UNIPROT %in% V(graph)$name) 
    interaction.graph <- igraph::graph_from_data_frame(interaction.tmp, directed = FALSE)
    metab.graph.tmp <- igraph::subgraph(METAB.graph, v = interaction.tmp$met)
    union.tmp <- igraph::union(metab.graph.tmp, interaction.graph)
    graph_w_metab <- igraph::union(graph, union.tmp)
    vertex_attr(graph_w_metab)$type <- ifelse(vertex_attr(graph_w_metab)$name %in% interaction.tmp$met, "Metabolite", vertex_attr(graph_w_metab)$type)
    return(graph_w_metab)
}

add_go_metabo_edges <- function(graph){
    va <- vertex_attr(graph) %>% as.data.frame()
    go.list <- va %>% filter(type == "GO") %>% pull(name) %>% as.character()
    cpd.list <- va %>% filter(type == "Metabolite") %>% pull(name) %>% as.character()
    tmp <- ko2go %>% na.omit() %>% filter(cpd %in% cpd.list, go %in% go.list)
    print(lapply(tmp, function(x) length(unique(x))) %>% unlist())
    if(nrow(tmp) > 0){
        tmp.graph <- igraph::graph_from_data_frame(tmp, directed = FALSE)
        return(igraph::union(tmp.graph, graph))
    }
    return(graph)
}

reformate_graphs <- function(graph, cluster){
    graph <- simplify(graph)
    tmp.cluster = cluster
    
    va <-  vertex_attr(graph) %>% as.data.frame() %>%
        mutate(type = case_when(
            name %in% gene_only ~ "Gene",
            name %in% (prot_only %>% filter(cluster ==tmp.cluster) %>% pull(name)) ~ "Prot",
            name %in% prot_ego ~ "Prot.ego",
            name %in% metab_only ~ "Metabolite",
            name %in% go_only ~ "go_term",
            str_detect(name, "GO:") ~ "go_term")) %>% 
        mutate(color = case_when(
            type == "Gene" ~ "lightblue",
            type == "Prot" ~ "orange",
            type == "Prot.ego" ~ "yellow",
            type == "Metabolite" ~ "green",
            type == "go_term" ~ "grey")) %>%
        
        mutate(size = case_when(
            color %in% c("lightblue", "orange", "grey") ~ 5,
            color %in% c("green", "yellow") ~ 2) %>% as.numeric()) %>%
        mutate(lab = case_when(
            color == "lightblue" ~ "Gene",
            color == "orange" ~ "Protein", 
            color == "yellow" ~ "Prot. 1D.N.", 
            color == "green" ~ "Metabolite", 
            color == "grey" ~ "Signif. GO terms")) %>% 
        mutate(cluster = cluster) %>%
        dplyr::select(type, name, color, cluster, size, lab) %>%
        purrr::map(as.character) %>% 
        as.data.frame(stringsAsFactors = FALSE) %>%
        mutate(size = as.numeric(size)) %>%
        mutate(type = factor(type, levels = c("Gene", "Prot", "Prot.ego","Metabolite", "go_term")))
    
    tmp.va <- va %>% arrange(type)
    
    va <- va %>% mutate(molecule = factor(name, levels = unique(tmp.va$name)),
                        lab = factor(lab, levels = unique(tmp.va$lab)))
    vertex_attr(graph) <- va
    return(graph)
}

homo_network_RWR <- function(graph, seed, AdjMatrixNorm = NULL, k = 15){
    seed <- seed[seed %in% V(graph)$name]
    if(is_empty(seed)){return(NULL)}
    network <- extract_component(graph, seed)[[1]]
    
    MultiplexObject <- RandomWalkRestartMH::create.multiplex(network, Layers_Name = c("component"))
    if(is.null(AdjMatrixNorm)){
        AdjMatrix <- compute.adjacency.matrix(MultiplexObject)
        AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)   
    }
    
    RWR_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm, MultiplexObject,seed)
    TopResults  <- create.multiplexNetwork.topResults(RWR_Results,MultiplexObject,k=k)
    return(TopResults)
}

extract_component <- function(graph, ids){
    log <- igraph::decompose(graph)
    index <- lapply(log, function(x) any(ids %in% V(x)$name)) %>% unlist()
    return(log[index])
}

RWR_with_go <- function(graph, list_of_go, k = 15, va.all){
    result_go <- list()
    for(go in list_go){
        if(go %in% names(V(graph))){
            result_homo <- homo_network_RWR(graph, seed = go, k = k)
            va <- vertex_attr(result_homo) %>% as.data.frame(stringsAsFactors = FALSE) %>% 
                left_join(va.all, by = "name") %>%
                mutate(color = ifelse(name == go, "red", color))
            vertex_attr(result_homo) <- va
            print(go)
            # if("metabo" %in% va$molecule){
            #     print(go)
            #     result_go[[go]] <- result_homo
            # }   
        }
        result_go[[go]] <- result_homo
        # } else {
        #     result_go[[go]] <- FALSE
        # } 
    }
    return(result_go)
}
mechanisms_annot <- function(graph, seed_name){
    if(class(graph) == "try-error"){
        return(NULL)
    }
    if(is.null(graph)){
        return(NULL)
    }
    vertex_attr(graph) %>% as.data.frame() %>%
        left_join(va.all) %>%
        mutate(seed  = seed_name)
}