library(tidyverse)

compound_reaction = read_tsv("/home/antoine/work/KEGGBD/ligand/compound/links/compound_reaction.list", col_names = c("cpd", "rn")) %>%
    mutate(cpd = str_remove(cpd, "cpd:"))  %>%#compound
    mutate(rn = str_remove(rn, "rn:"))

reaction_mapformula = read.table("/home/antoine/work/KEGGBD/ligand/reaction/reaction_mapformula.lst", sep = ":", col.names = c("rn", "N", "formula")) %>%
    mutate(p1 = formula %>% str_trim() %>% str_split('=') %>% map_chr(~.x[1]) %>% str_remove(" <") %>% str_trim %>% str_split(" \\+ ")) %>%
    mutate(p2 = formula %>% str_trim() %>% str_split('=') %>% map_chr(~.x[2]) %>% str_remove("> ") %>% str_trim %>% str_split(" \\+ "))

# reaction_mapformula_df <- matrix(ncol = ncol(reaction_mapformula), dimnames = list("",names(reaction_mapformula))) %>%
#     as.data.frame() %>% na.omit()
# reaction_mapformula_df = matrix(ncol = 3, dimnames = list("", c("p1", "p2", "rn"))) %>% as.data.frame() %>% na.omit()
# for( i in 1:nrow(reaction_mapformula)){
#     # p1
#     p1 <- reaction_mapformula[i,"p1", drop = TRUE][[1]]
#     p2 <- reaction_mapformula[i,"p2", drop = TRUE][[1]]
#     rn <- reaction_mapformula[i,"rn", drop = TRUE][[1]]
#     combinations <- rbind(expand.grid(p1,p1), expand.grid(p1,p2), expand.grid(p2,p2)) %>% set_names(c("p1", "p2")) %>%
#         filter(p1 != p2) %>%
#         mutate(rn = rn)
#     # for(j in nrow(combinations)){
#     #     tmp <- reaction_mapformula[i,] %>% dplyr::select(-c(p1,p2))
#     #     tmp$p1 <- as.character(combinations$p1[j])
#     #     tmp$p2 <- as.character(combinations$p2[j])
#     # }
#     # rbind(reaction_mapformula_df, tmp)
#     reaction_mapformula_df <- rbind(reaction_mapformula_df, combinations)
# }
# save(reaction_mapformula_df, file = "../vrac/reaction_mapformula_df.rda")
load("/home/antoine/Documents/TO2//vrac/reaction_mapformula_df.rda")
    # graph, make_ego_graph(compound)

ko_reaction <- read_tsv("/home/antoine/work/KEGGBD/ko/ko_reaction.list", col_names = c("ko","rn"))

kegg_compound <- read_tsv("~/Documents/TO2/vrac/kegg_compound.tsv")

# from compound_id to reaction id
compound = "C00805"
get_reaction_component <- function(compound){
    # return all compound connected to compound (ex P1) by reaction
    # P1 => P2: (P1,P2)
    # P1 + P2 => P3: (P1,P3)
    # P1 + P2 => P3 + P4: (P1,P3), (P1,P4)
    reaction.map_filtered <- compound_reaction %>% filter(cpd == compound) %>%
        left_join(reaction_mapformula, by = c("rn"='rn'))
    if(nrow(reaction.map_filtered) == 0){return(FALSE)}
    
    TMP <- matrix(ncol = ncol(reaction.map_filtered), dimnames = list("",names(reaction.map_filtered))) %>%
        as.data.frame() %>% na.omit()
    
    for( i in 1:nrow(reaction.map_filtered)){
        if(compound %in% tmp[1,"p1", drop = TRUE][[1]]){
            for(j in seq_along(reaction.map_filtered[i,"p2", drop = TRUE][[1]])){
                tmp <- reaction.map_filtered[i,]  %>% dplyr::select(-c(p1,p2))
                tmp[1, "p1"] <- compound
                tmp[1, "p2"] <- reaction.map_filtered[i,"p2", drop = TRUE][[1]][j]
                TMP <- rbind(TMP, tmp)
            }
        }
        else{ # compound in p2 so combine each p1 with compound
            for(j in seq_along(reaction.map_filtered[i,"p1", drop = TRUE][[1]])){
                tmp <- reaction.map_filtered[i,]  %>% dplyr::select(-c(p1,p2))
                tmp[1, "p1"] <- reaction.map_filtered[i,"p1", drop = TRUE][[1]][j]
                tmp[1, "p2"] <- compound
                TMP <- rbind(TMP, tmp)
            }
        }
    }
    return(TMP)
}

#get_reaction_component(compound = 'C22202')

get_compound_ID <- function(name){
    pattern = paste0("\\b(?i)", name, "\\b")
    kegg_compound_short <- kegg_compound %>% dplyr::select("ENTRY", "NAME") 
    kegg_compound_short %>% filter(str_detect(NAME, pattern))
}

# get_compound_ID("241")

# kegg ko <-> rn <-> zma
ko_reaction.list <- read_tsv("~/work/KEGGBD/ko/ko_reaction.list", col_names = c("ko", "rn"))
ko_enzyme.list <- read_tsv("~/work/KEGGBD/ko/ko_enzyme.list", col_names = c("ko", "ec"))
ko_gene.list <- read_tsv("~/work/KEGGBD/ko/ko_genes.list", col_names = c("ko", "gn"))

ko_reaction_enzyme_gene <- ko_reaction.list %>% full_join(ko_enzyme.list) %>% full_join(ko_gene.list)
tmp <- ko_reaction_enzyme_gene %>% filter(str_detect(gn, "zma:")) %>% # filter only zea mais and 
    mutate(gn = str_remove(gn,"zma:") %>% as.numeric)  %>% # and keep only the 6 last digit == unique identifier
    mutate(rn = str_remove(rn,"^rn:"))

# tmp %>% left_join(genes_v3v4)

# change name of G -> GG
genes_v3v4 <- read_tsv("/home/antoine/Documents/TO2/CS2/ressources/gene_model_xref_v4.txt", skip = 4)

lapply(genes_v3v4, function(x) sum(tmp$gn %in% x)/length(tmp$gn))  # Entrez
genes_v3v4 <- genes_v3v4 %>% dplyr::select(c("v4_gene_model", "v3_gene_model", "Entrez"))
rn_to_gene <- tmp %>% left_join(genes_v3v4, by = c("gn"="Entrez"))

# # add our genes
# as.data.frame(names(V(G))) %>% set_names("G") %>% 
#     left_join(genes_v3v4, by=c("G"="v3_gene_model"))  %>%
#     left_join(tmp, by = c("Entrez"="gn")) %>%
#     filter(!is.na(ec))
    

# add go to ko
ko_go.list <- read_tsv(file = "~/work/KEGGBD/ko/ko_go.list", col_names = c("ko", "go"))
# ko2go <- ko_reaction_enzyme_gene %>% full_join(ko_go.list) %>% 
#     mutate(rn = str_remove(rn, "rn:")) %>%
#     left_join(reaction_mapformula_df) %>% dplyr::select(go, p1, p2) %>%
#     gather(p, cpd, -c(go)) %>%
#     dplyr::select(go, cpd) %>%
#     unique %>%
#     mutate(go = str_replace(go, "go:", "GO:"))
#save(ko2go, file = "../vrac/ko2go.rda")
    
save(compound_reaction, reaction_mapformula, reaction_mapformula_df, ko_reaction, kegg_compound, get_reaction_component, get_compound_ID, ko_reaction.list, ko_enzyme.list, ko_gene.list, ko_reaction_enzyme_gene, genes_v3v4, rn_to_gene, ko_go.list, file = "./analysis/kegg_explore_res.rda")
