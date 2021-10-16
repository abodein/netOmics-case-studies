remove.low.cv <- function(X, cutoff = 0.5){
    # var.coef
    X <- as.data.frame(X)
    cv <- unlist(lapply(X, function(x) abs(sd(x)/mean(x))))
    return(list("X" = X[,which(cv>cutoff)], "cv" = cv, "cutoff" = cutoff ))
}

remove.fold.change <- function(X, cutoff, hist.title = "", plot = TRUE){
    X <- as.data.frame(X)
    fc.func <- function(x){
        B <- max(x+1)
        A <- min(x+1)
        return(log(B)-log(A))
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

wrapper.lmms.5tp <- function(X, time, sampleID){
    pspline <- lmms::lmmSpline(data = X, time = time,
                               sampleID = sampleID, deri = FALSE,
                               basis = "p-spline", timePredict = seq(2,48, by =2),
                               numCores = 4,
                               keepModels = TRUE)
    
    cubic <- lmms::lmmSpline(data = X, time = time,
                             sampleID = sampleID, deri = FALSE,
                             basis = "cubic", timePredict = seq(2,48, by =2),
                             numCores = 4,
                             keepModels = TRUE)
    
    cubic.pspline <- lmms::lmmSpline(data = X, time = time,
                                     sampleID = sampleID, deri = FALSE,
                                     basis = "cubic", timePredict = seq(2,48, by =2),
                                     numCores = 4,
                                     keepModels = TRUE)
    
    RESULT = list("cubic"= cubic, "pspline"=pspline, "cubic.pspline" =cubic.pspline)
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


fit.spline.predict <- function(X, spar, time, time.pred){
    fit <- smooth.spline(y = X, x=time, spar = spar)
    fit.predict <- predict(fit, time.pred)
    return(fit.predict$y)
}
get_graph <- function(X){
    mim <- build.mim(X,estimator="spearman")
    adj <- aracne(mim)
    G <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
    dg <- degree(G)
    
    Isolated = which(degree(G)==0)
    G2 = delete.vertices(G, Isolated)
    
    return(list(dg = dg, G = G, G2= G2))
}
get_stats_graph <- function(graph){
    node.c <- sum(degree(graph) != 0)
    node.i <- sum(degree(graph) == 0)
    edge <- ecount(graph)
    density <- igraph::edge_density(graph = graph)
    return(list(node.c = node.c, node.i = node.i, edge = edge, density = density))
}


plot_grn_ppi_neighbors <- function(sub, PPIM, order = 1, plot.title, plot = TRUE){
    # 1. get node from sub connected to PPIM
    node <- intersect(names(V(PPIM)), names(V(sub)))
    PPIM_sub <- subgraph(PPIM, v = node)
    
    
    # 2. make ego graph from that nodes
    ego_G <- make_ego_graph(graph = PPIM, order = order, nodes = node, mode = "all")
    ego_G <- do.call(igraph::union, ego_G)
    
    # 3. merge sub and PPIM neighbors
    super <- igraph::union(sub, ego_G)
    
    # 4. change attribute
    # -------------
    # n1 = only node from GRN
    n1 <- setdiff(names(V(sub)), names(V(PPIM_sub)))
    # n2 = node from GRN connected to PPIM
    n2 <- names(V(PPIM_sub))
    # n3 only neighbors
    n3 <- setdiff(names(V(super)), names(V(sub)))
    
    # node color
    vertex_attr(super)$color <- case_when(
        vertex_attr(super)$name %in% n1 ~ "lightblue",
        vertex_attr(super)$name %in% n2 ~ "orange",
        vertex_attr(super)$name %in% n3 ~ "yellow")
    # node size
    vertex_attr(super)$size <- case_when(
        vertex_attr(super)$name %in% n1 ~ 4,
        vertex_attr(super)$name %in% n2 ~ 5,
        vertex_attr(super)$name %in% n3 ~ 2)
    
    # remove unconnected nodes
    super.simplified = remove.isolated.nodes(super, simplify = TRUE)
    
    # 5. statistics
    graphs <- list("GRN"=sub, "PPIM_in_GRN"=PPIM_sub, "Neighbors"=super) 
    stats <- imap(graphs, ~{
        list("Conn.Node" = length(degree(.x) != 0),
             "Edges" = length(E(.x)),
             "Isol.Node" = sum(degree(.x) == 0)) %>%
            as.data.frame()
    }) %>% do.call(what = "rbind") %>%
        rownames_to_column("title") %>%
        mutate(position = c("topleft","bottomleft","bottomright"),
               color = c("lightblue", "orange", "yellow"))
    
    
    # 6. plot
    if(plot == TRUE){
        plot(super.simplified, vertex.label = NA)
        title(plot.title)
        for(i in 1:3){
            legend(stats$position[i], 
                   legend = c(paste0("Conn.Nodes = ", stats$Conn.Node[i]),
                              paste0("Edges = ", stats$Edges[i]),
                              paste0("Isol.Nodes = ", stats$Isol.Node[i])),
                   # col = stats$color[i], 
                   # pch = c(19,NA,NA),
                   title = stats$title[i],
                   bg=stats$color[i])
        }
        legend("topright", legend = c(paste0("Conn.Nodes = ", length(V(super.simplified))),
                                      paste0("Edges = ", length(E(super.simplified))))) 
    }
    return(list("graph" = super, "simplify" = super.simplified, "stat" = stats))
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

get_metabo_graph <- function(clu, metabo_w_compound_id){
    print(clu)
    cpd_of_interest <- metabo_w_compound_id %>% filter(cluster == clu) %>% pull(Compound_ID) %>% str_split(";") %>% unlist %>% na.omit
    metabo.graph.clu <- reaction_mapformula_df_filtered %>% filter(p1 %in% cpd_of_interest | p2 %in% cpd_of_interest) %>%
        graph_from_data_frame(directed = FALSE)
    vertex_attr(metabo.graph.clu)$color <- case_when(
        vertex_attr(metabo.graph.clu)$name %in%  cpd_of_interest ~ "purple",
        !(vertex_attr(metabo.graph.clu)$name %in%  cpd_of_interest) ~ "green"
    )
    if(vcount(metabo.graph.clu) != 0)
    {
        plot(metabo.graph.clu, vertex.label = NA, vertex.size = 6)
        title(main = paste0("Cluster ", clu))
        legend("bottomleft", 
               legend = c(paste0("V(ours): ", sum(vertex_attr(metabo.graph.clu)$name %in%  cpd_of_interest)),
                          paste0("V(reactions): ", sum(!(vertex_attr(metabo.graph.clu)$name %in% cpd_of_interest))),
                          paste0("E: ", length(E(metabo.graph.clu)))),
               col = c("purple", "green", "grey"),
               pch = c(19,19,NA),
               lty = c(NA,NA,1))
    }
    return(metabo.graph.clu)
}

add_metabo_layer <- function(GRN, GRN_P, METABO_M, interaction_df, cluster){
    # GRN = Gene + protein => list of node to connect
    #GRN <- GRN_PPIM$`1`$graph
    # full.metab.G <- full metabo layer, internal connection
    # interaction_df = compound <-> zma
    
    node.GRN.P <- names(V(GRN_P))
    interaction_df_filtered <- interaction_df %>% filter(v3_gene_model %in% node.GRN.P)
    
    # build a junction graph from interaction
    junction_graph <- interaction_df_filtered %>% 
        dplyr::select(v3_gene_model, p1,p2) %>% gather(p, node, -v3_gene_model) %>%
        dplyr::select(v3_gene_model, node) %>% 
        filter(v3_gene_model %in% node.GRN.P) %>%
        filter(node %in% names(V(METABO_M))) %>%
        igraph::graph_from_data_frame(directed = FALSE)
    
    # union GRN, full.metabo.G
    full_graph <- igraph::union(junction_graph, METABO_M) %>%
        igraph::union(GRN_P)
    
    # stats
    stats = list(
        stat.gene.gene = vertex_attr(GRN_P) %>% as.data.frame() %>% filter(color == "lightblue") %>% pull(name) %>% as.character(),
        stat.gene.prot = vertex_attr(GRN_P) %>% as.data.frame() %>% filter(color == "orange") %>% pull(name) %>% as.character(),
        stat.prot.ego = vertex_attr(GRN_P) %>% as.data.frame() %>% filter(color == "yellow") %>% pull(name) %>%  as.character(),
        stats.metabo.M = names(V(METABO_M))[names(V(METABO_M)) %in% metabo_nodes],
        stats.metabo.R = names(V(METABO_M))[!(names(V(METABO_M)) %in% metabo_nodes)],
        stats.gene.M = interaction_df_filtered$v3_gene_model)
    
    vertex_attr(full_graph)$color <- case_when(
        vertex_attr(full_graph)$name %in% stats$stat.gene.gene ~ "lightblue",  # only genes
        vertex_attr(full_graph)$name %in% stats$stat.gene.prot ~ "orange",  # gene/prot
        vertex_attr(full_graph)$name %in% stats$stat.prot.ego ~ "yellow", # prot only ego 1
        vertex_attr(full_graph)$name %in% stats$stats.metabo.M ~ "purple",
        vertex_attr(full_graph)$name %in% stats$stats.metabo.R ~ "green"
    )
    
    #vertex_attr(full_graph) %>% as.data.frame() %>% filter(is.na(color)) %>% dplyr::select(name)
    
    full_graph.simplified = remove.isolated.nodes(full_graph, simplify = TRUE)
    plot(full_graph.simplified, vertex.size = 4, vertex.label = NA)
    legend("bottomleft", 
           legend = c(paste0("Gene/Gene: ",sum(names(V(full_graph.simplified)) %in% stats$stat.gene.gene) ,"(",length(stats$stat.gene.gene),")" ),
                      paste0("Gene/Prot: ",sum(names(V(full_graph.simplified)) %in% stats$stat.gene.prot) ,"(",length(stats$stat.gene.prot),")"),
                      paste0("Prot/Ego: ",sum(names(V(full_graph.simplified)) %in%  stats$stat.prot.ego) ,"(",length( stats$stat.prot.ego),")"),
                      paste0("Metabo(M): ",sum(names(V(full_graph.simplified)) %in% stats$stats.metabo.M) ,"(",length(stats$stats.metabo.M),")"),
                      paste0("Metabo(R): ",sum(names(V(full_graph.simplified)) %in% stats$stats.metabo.R) ,"(",length(stats$stats.metabo.R),")"),
                      paste0("E: ",length(E(full_graph.simplified)))),
           pch = c(19,19,19,19,19, NA),
           col = c("lightblue", "orange", "yellow", "purple", "green"),
           lty = c(NA, NA, NA, NA, NA, 1))
    title(main = c(paste0("Cluster ", cluster)))
    
    #stats[["decompose"]] <- decompose(full_graph, mode = "strong")
    return(list(graph = full_graph, simplify =  full_graph.simplified, stats = stats))
}

get_enrichment <- function(graph){
    query_raw = vertex_attr(graph) %>% as.data.frame() %>% 
        filter(color %in% c("lightblue", "orange")) %>% 
        left_join(genes_v3v4, by = c("name"="v3_gene_model")) %>% dplyr::select("v4_gene_model") %>% pull %>%
        as.character() %>%
        na.omit
    res_raw <- gprofiler2::gost(query = query_raw, organism = "zmays")
    
    query_extended = vertex_attr(graph) %>% as.data.frame() %>% 
        filter(color %in% c("lightblue", "orange", "yellow")) %>% 
        left_join(genes_v3v4, by = c("name"="v3_gene_model")) %>% dplyr::select("v4_gene_model") %>% pull %>%
        as.character() %>%
        na.omit
    res_extended <- gprofiler2::gost(query = query_extended, organism = "zmays")
    return(list("query_raw" = query_raw, "query_extended"=query_extended, "res_raw" = res_raw, "res_extended" = res_extended))
}

add_go_term_layer <- function(graph){
    go_enrich_graph_tmp <- 
        # goterm_map.signif %>%
        goterm_map_all.signif %>% # choose all go terms
        dplyr::select(c(target, v3_gene_model)) %>% 
        na.omit() %>%
        filter(v3_gene_model %in% names(V(graph))) %>%
        igraph::graph.data.frame(directed = FALSE) 
    vertex_attr(go_enrich_graph_tmp)$color <- ifelse(str_detect(vertex_attr(go_enrich_graph_tmp)$name, "GO:"), "grey", "lightblue")
    graph_4_layers <- igraph::union(graph, go_enrich_graph_tmp)
    return(graph_4_layers)
}

reformate_graphs <- function(graph, cluster){
    graph <- simplify(graph)
    va <- vertex_attr(graph) %>% as.data.frame()%>%
        mutate(color_1 = as.character(color_1)) %>%
        mutate(color_2 = as.character(color_2)) %>%
        mutate(color = ifelse(is.na(color_1), "grey", color_1)) %>%  # GO terms
        mutate(type = case_when(
            color == "lightblue" ~ "gene_gene",
            color == "orange" ~ "gene_prot",
            color == "yellow" ~ "prot_ego",
            color == "purple" ~ "metabo_targeted",
            color == "green" ~  "metabo_reaction",
            color == "grey" ~ "go_term"
        )) %>%
        mutate(molecule = case_when(
            color %in% c("lightblue", "orange") ~ "gene",
            color %in% c("yellow") ~ "prot",
            color %in% c("purple", "green") ~ "metabo",
            color %in% c("grey") ~ "go_term",
        )) %>%
        mutate(size = case_when(
            color %in% c("lightblue", "orange", "purple", "grey") ~ 5,
            color %in% c("green", "yellow") ~ 2) %>% as.numeric()) %>%
        mutate(lab = case_when(
            color == "lightblue" ~ "Gene",
            color == "orange" ~ "Gene/Prot", 
            color == "yellow" ~ "Prot. 1D.N.", 
            color == "purple" ~ "Metabolite", 
            color == "green" ~ "Met. in reaction", 
            color == "grey" ~ "Signif. GO terms"
        )) %>%
        mutate(cluster = cluster) %>%
        dplyr::select(type, name, color, molecule, cluster, size, lab) %>%
        purrr::map(as.character) %>% 
        as.data.frame(stringsAsFactors = FALSE) %>%
        mutate(size = as.numeric(size)) %>%
        mutate(type = factor(type, levels = c("gene_gene", "gene_prot", "prot_ego","metabo_targeted", "metabo_reaction", "go_term")))
    
    tmp.va <- va %>% arrange(type)
    
    va <- va %>% mutate(molecule = factor(molecule, levels = unique(tmp.va$molecule)),
                        lab = factor(lab, levels = unique(tmp.va$lab)))
    vertex_attr(graph) <- va
    return(graph)
}
get_component <- function(graph, elt){
    composition <- igraph::decompose(graph)
    index <- lapply(composition, function(x) {elt %in% names(V(x))}) %>% unlist %>% which
    if(!is_null(index)){
        return(composition[[i]])
    }
    return()
}
extract_component <- function(graph, ids){
    log <- igraph::decompose(graph)
    index <- lapply(log, function(x) any(ids %in% V(x)$name)) %>% unlist()
    return(log[index])
}

homo_network_RWR <- function(graph, seed, AdjMatrixNorm = NULL, k = 15){
    seed <- seed[seed %in% V(graph)$name]
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

RWR_with_go <- function(graph, list_of_go, k = 15){
    va.all <- vertex_attr(graph) %>% as.data.frame()
    result_go <- list()
    for(go in list_of_go){
        if(go %in% names(V(graph))){
            result_homo <- homo_network_RWR(graph, seed = go, k = k)
            va <- vertex_attr(result_homo) %>% as.data.frame(stringsAsFactors = FALSE) %>% 
                left_join(va.all, by = "name") %>%
                mutate(color = ifelse(name == go, "red", color))
            vertex_attr(result_homo) <- va
            result_go[[go]] <- result_homo
        }
    }
    return(result_go)
}

stat_table <- function(list_g, va.all){
    res <- imap_dfr(list_g, ~.x %>% 
                        vertex_attr %>%
                        mutate(type = factor(type, levels = levels(va.all$type))) %>% 
                        pull(type) %>% table %>% as.data.frame() %>% 
                        column_to_rownames(".") %>% set_names(.y) %>%
                        t %>% as.data.frame() %>% rownames_to_column("seed")) %>%
        left_join(map_dfr(.$seed %>% unique, ~get_go_info(.x)), by = c("seed"="name")) 
    return(res)
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
                description = paste0(term, " (", ontology, ")"),
                source = ontology) %>% 
               as.data.frame())
}

find_annotation <- function(graph, seed){
    seed <- seed[seed %in% V(graph)$name]
    if(is_empty(seed)){
        return(NULL)
    }
    network <- extract_component(graph, seed)[[1]]
    if(vcount(network) == 1){
        return(NULL)
    }
    
    MultiplexObject <- RandomWalkRestartMH::create.multiplex(network, Layers_Name = c("component"))
    AdjMatrix <- compute.adjacency.matrix(MultiplexObject)
    AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)   
    
    
    RWR_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm, MultiplexObject,seed)
    return(RWR_Results)
}
in_between <- function(graph, seed, k = 10){
    seed <- seed[seed %in% V(graph)$name]
    if(is_empty(seed)){
        return(NULL)
    }
    network <- extract_component(graph, seed)[[1]]
    if(vcount(network) == 1){
        return(NULL)
    }
    MultiplexObject <- RandomWalkRestartMH::create.multiplex(network, Layers_Name = c("component"))
    AdjMatrix <- compute.adjacency.matrix(MultiplexObject)
    AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)   
    
    
    RWR_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm, MultiplexObject,seed)
    TopResults  <- create.multiplexNetwork.topResults(RWR_Results,MultiplexObject,k=k)
    return(TopResults)
}
annot_RWR_res <- function(RWR){
    
    if(is(RWR,"try-error") || is.null(RWR$RWRM_Results)){
        return(as.data.frame(list(NodeNames = NA, Score = NA, position = NA)))
    }
    RWR$RWRM_Results %>% #left_join(va.all, by = c("NodeNames"= "name")) %>%
        mutate(position = seq(1,nrow(.))) %>% #pull(type) %>% unique
        #filter(type == "go_term") %>%
        filter(NodeNames %>% str_detect("GO:")) %>%
        dplyr::top_n(n = -1, wt = position) %>%
        dplyr::select(NodeNames, Score, position)
    # pull(position)
}

# # 1) go
# go.info <- map_dfr(read_tsv("../analysis/gene_goterm_map.tsv")$target %>% unique, ~get_go_info(.x)) %>% 
#     dplyr::select(-source)
# 
# 
# # 2) Gene
# gene.info <- read_tsv("../ressources/genes_all.txt")[c(2,3,4)] %>% na.omit %>%
#     purrr::set_names(c("name", "lab", "description"))
# 
# # 3) Met
# met.info <- read_tsv("../../vrac/kegg_compound.tsv") %>%
#     dplyr::select(c(ENTRY, NAME)) %>%
#     mutate(NAME = str_remove(NAME, ";;.*")) %>% 
#     set_names(set_names(c("name", "lab"))) %>%
#     mutate(description = lab)
# 
# rbind.info <- rbind(go.info, met.info, gene.info)


detect_intersection <- function(graph, seed_name, va.all){
    # va.all <- vertex_attr(graph) %>% as.data.frame()
    va <- vertex_attr(graph) %>% as.data.frame() %>%
        left_join(cluster.info, by = c("name"="molecule")) %>%
        na.omit  # only sequenced molecule
    
    cluster_seed <- va %>% filter(name == seed_name) %>% pull(cluster)
    tmp <- va %>% filter(!is.na(cluster)) %>%
        mutate(cluster.seed = cluster_seed) %>%
        mutate(cluster.mol = ifelse(name == seed_name, cluster.seed, cluster)) %>%
        mutate(between = (cluster.seed != cluster.mol)) 
    flag <- tmp %>% pull(between) %>% any
    if(flag){
        to_return <- vertex_attr(graph) %>% as.data.frame() %>% 
            left_join(va.all, by = c("name"="name")) %>%
            filter(type == "go_term")
        # add go lab
        if(any(to_return$type == "go_term")){
            to_return <- left_join(to_return, map_dfr(to_return$name, ~get_go_info(.x)), by = c("name"="name")) %>%
                dplyr::select(name, description)  %>%
                mutate(seed = seed_name, cluster.seed = cluster_seed, cluster.mol = NA)
        }
        to_return <- tmp %>% filter(between) %>% mutate(seed = seed_name) %>%
            dplyr::select(seed, name, cluster.seed, cluster.mol) %>%
            mutate(description = NA) %>%
            rbind(to_return)
        return(to_return)
    }
    return(NULL)
}

# graph <- RWR_meca_res$all$`GO:0009698`
plot_graph_with_legend_2 <- function(graph, ...){
    layout(matrix(c(1,2), ncol =2, byrow = TRUE))
    # vertex attribute
    va <- vertex_attr(graph) %>% as.data.frame()
    va.tmp1 <- va %>% dplyr::select(type, color, size, lab) %>% group_by(type, color, size, lab) %>%
        summarise(N = n())
    legend.text <- va.tmp1 %>% mutate(legend.text = paste0(lab, ": ", N)) %>% pull(legend.text)
    legend.col <- va.tmp1$color %>% as.character()
    legend.size <- va.tmp1$size / max(va.tmp1$size) * 2
    
    # change lab 
    va.tmp2 <-  va %>% left_join(rbind.info, by = "name") %>% mutate(ID = name) %>%
        mutate(name = ifelse(!is.na(lab.y),as.character(lab.y), ID))
    vertex_attr(graph) <- va.tmp2
    
    # plot
    plot(x = graph, ...)
    #plot(graph, layout=layout_with_fr)
    
    
    # title
    seed = va %>% filter(color == "red") %>% left_join(go.info, by = "name")
    title(main = seed$name, sub = seed$description, cex.main = 2, cex.sub = 1.5)
    
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    
    # legend 1
    legend("topright", legend = legend.text, col = legend.col, pch = 19, pt.cex = legend.size)
    
    # legend 2
    lab2 <- va.tmp2 %>% dplyr::select(color, size, lab.y, ID, description) %>%
        mutate(lab.y = as.character(lab.y)) %>%
        mutate(lab.y = ifelse(lab.y == ID, "", lab.y)) %>%
        mutate(lab.y = ifelse(is.na(lab.y), "", paste0(lab.y, " "))) %>%
        mutate(description = as.character(description)) %>%
        mutate(description = ifelse(is.na(description), "", description)) %>%
        filter(description != "") %>%
        mutate(lab = paste0(lab.y, ID, " ", description)) %>%
        dplyr::select(color, size, lab)
    
    legend2.col <- lab2$color
    legend2.text <- lab2$lab
    legend2.size <- lab2$size / max(lab2$size)
    legend("bottomright", legend = legend2.text, 
           col = legend2.col, 
           pch = 19, 
           pt.cex = legend2.size, 
           cex = 0.7)
}
