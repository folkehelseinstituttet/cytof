sti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Panel 1")
mapper <- list.files(sti)
mapper <- mapper[!grepl(".docx", mapper)]
mapper <- mapper[!grepl(".xlsx", mapper)]

seeds <- list.files(fs::path(sti,mapper[1]))
temp_unike <- read.csv2(fs::path(sti,mapper[1], seeds[1], paste0("unike_klustre", gsub("seed ", "", seeds[1]), ".csv")))
info <- temp_unike[,c("kort_filnavn", "status", "age", "sex", "statusTid", "tid")]

klustre <- data.frame(id = info[,1])

qTemp <- read.csv2(fs::path(sti,mapper[1], seeds[1], paste0("q5_per_kluster_k_", 10, "_seed", gsub("seed ", "", seeds[1]), mapper[1], ".csv")))
qTom2 <- data.frame(matrix(NA, ncol = ncol(qTemp), nrow = 2))
colnames(qTom2) <- colnames(qTemp)
rownames(qTom2) <- c("slett 1", "slett 2")
qTom2 <- qTom2[,-1]
qTom <- t(qTom2)
qTom <- data.frame(qTom)
qTom$marker <- rownames(qTom)

q5_trans <- qTom
q10_trans <- qTom
q25_trans <- qTom
q50_trans <- qTom
q75_trans <- qTom
q90_trans <- qTom
q95_trans <- qTom


strsplit4 <- function(x,y){
  return(strsplit(x,y)[[1]][4])
}

strsplit6 <- function(x,y){
  return(strsplit(x,y)[[1]][6])
}


for(i in 1:length(mapper)){
  seeds <- list.files(fs::path(sti,mapper[i]))
  for(j in 1:length(seeds)){
    temp_unike <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("unike_klustre", gsub("seed ", "", seeds[j]), ".csv")))
    if(file.info(fs::path(sti,mapper[i], seeds[j], paste0("sign_klustre_", gsub("seed ", "", seeds[j]), ".csv")))$size > 5){
      temp_sign <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("sign_klustre_", gsub("seed ", "", seeds[j]), ".csv")))
      colnames(temp_sign) <- c("X", "sign")
      temp <- temp_unike[, as.character(temp_sign$sign)] 
      if(!is.null(dim(temp))){
        colnames(temp) <- paste0(mapper[i], "_", colnames(temp))
        antall <- matrix(rep(temp_unike$antall_cells_totalt, ncol(temp_sign)), byrow = FALSE, ncol = ncol(temp_sign))
        klustre <- cbind(klustre, temp/antall)
      } else {
        antall <- temp_unike$antall_cells_totalt
        klustre <- cbind(klustre, temp/antall)
        colnames(klustre)[ncol(klustre)] <- paste0(mapper[i], "_", as.character(temp_sign$sign))
      }
      
      unike_k <-  unique(unlist(lapply(as.character(temp_sign$sign),strsplit4, y = "_")))
      for(k in unike_k){
        #q5
        temp_q5 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("q5_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q5 <- temp_q5[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q5_temp <- t(temp_q5)
        q5_temp <- data.frame(q5_temp)
        colnames(q5_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q5_temp)))
        q5_temp$marker <- rownames(q5_temp)
        n_marker <- which(colnames(q5_temp) == "marker")
        q5_trans <- merge(q5_trans, q5_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
        #q10
        temp_q10 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("q10_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q10 <- temp_q10[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q10_temp <- t(temp_q10)
        q10_temp <- data.frame(q10_temp)
        colnames(q10_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q10_temp)))
        q10_temp$marker <- rownames(q10_temp)
        n_marker <- which(colnames(q10_temp) == "marker")
        q10_trans <- merge(q10_trans, q10_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
        #q25
        temp_q25 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("q25_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q25 <- temp_q25[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q25_temp <- t(temp_q25)
        q25_temp <- data.frame(q25_temp)
        colnames(q25_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q25_temp)))
        q25_temp$marker <- rownames(q25_temp)
        n_marker <- which(colnames(q25_temp) == "marker")
        q25_trans <- merge(q25_trans, q25_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
        #q50
        temp_q50 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("medians_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q50 <- temp_q50[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q50_temp <- t(temp_q50)
        q50_temp <- data.frame(q50_temp)
        colnames(q50_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q50_temp)))
        q50_temp$marker <- rownames(q50_temp)
        n_marker <- which(colnames(q50_temp) == "marker")
        q50_trans <- merge(q50_trans, q50_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
        #q5
        temp_q75 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("q75_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q75 <- temp_q75[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q75_temp <- t(temp_q75)
        q75_temp <- data.frame(q75_temp)
        colnames(q75_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q75_temp)))
        q75_temp$marker <- rownames(q75_temp)
        n_marker <- which(colnames(q75_temp) == "marker")
        q75_trans <- merge(q75_trans, q75_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
        #q5
        temp_q90 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("q90_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q90 <- temp_q90[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q90_temp <- t(temp_q90)
        q90_temp <- data.frame(q90_temp)
        colnames(q90_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q90_temp)))
        q90_temp$marker <- rownames(q90_temp)
        n_marker <- which(colnames(q90_temp) == "marker")
        q90_trans <- merge(q90_trans, q90_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
        #q5
        temp_q95 <- read.csv2(fs::path(sti,mapper[i], seeds[j], paste0("q95_per_kluster_k_", k, "_seed", gsub("seed ", "", seeds[j]), mapper[i], ".csv")))
        temp_q95 <- temp_q95[,-1]
        unike_temp_sign_k <- as.character(temp_sign$sign)[grepl(paste0("k_", k), as.character(temp_sign$sign))]
        unike_k_k <- as.integer(unique(unlist(lapply(as.character(unike_temp_sign_k),strsplit6, y = "_"))))
        q95_temp <- t(temp_q95)
        q95_temp <- data.frame(q95_temp)
        colnames(q95_temp) <- paste0(mapper[i], "_seed_", gsub("seed ", "", seeds[j]), "_k_", k, "_kluster_", gsub("X", "", colnames(q95_temp)))
        q95_temp$marker <- rownames(q95_temp)
        n_marker <- which(colnames(q95_temp) == "marker")
        q95_trans <- merge(q95_trans, q95_temp[, c(unike_k_k, n_marker)], by = "marker",all = T)
      }
    }
  }
}

rownames(q5_trans) <- q5_trans$marker
q5_trans <- q5_trans[,-c(1,2,3)]
q5_all <- t(q5_trans)

rownames(q10_trans) <- q10_trans$marker
q10_trans <- q10_trans[,-c(1,2,3)]
q10_all <- t(q10_trans)

rownames(q25_trans) <- q25_trans$marker
q25_trans <- q25_trans[,-c(1,2,3)]
q25_all <- t(q25_trans)

rownames(q50_trans) <- q50_trans$marker
q50_trans <- q50_trans[,-c(1,2,3)]
q50_all <- t(q50_trans)

rownames(q75_trans) <- q75_trans$marker
q75_trans <- q75_trans[,-c(1,2,3)]
q75_all <- t(q75_trans)

rownames(q90_trans) <- q90_trans$marker
q90_trans <- q90_trans[,-c(1,2,3)]
q90_all <- t(q90_trans)

rownames(q95_trans) <- q95_trans$marker
q95_trans <- q95_trans[,-c(1,2,3)]
q95_all <- t(q95_trans)







temp0 <- cor(klustre[,-1])
temp <- temp0

like <- list()


ii <- 0
brukt <- NULL

while(length(brukt) < ncol(temp0)){
  ii <- ii + 1
  temp <- temp0[!colnames(temp0) %in% brukt, !colnames(temp0) %in% brukt]
  temp1 <- which(temp[,1] > 0.95)
  temp2 <- NULL
  for(j in 1:length(temp1)){
    temp2 <- c(temp2, which( temp[,temp1[j]] > 0.95))
  }
  temp3 <- NULL
  for(j in 1:length(temp1)){
    temp3 <- c(temp3, which( temp[,temp2[j]] > 0.95))
  }
  print(ii)
  print(unique(names(temp1)))
  print(unique(names(temp2)))
  print(unique(names(temp3)))
  like[[ii]] <- unique(c(names(temp1), names(temp2), names(temp3)))
  brukt <- c(brukt, like[[ii]])
}







scriptPath_UnsupAnalysisUt <- fs::path(scriptPath_UnsupAnalysis, "resultat_panel1_cor095")

meanAbundsLike <- matrix(NA, ncol = length(like), nrow = nrow(klustre))



for(i in 1:length(like)){
  if(length(like[[i]]) > 1){
    meanAbundsLike[,i] <- apply(klustre[, like[[i]]],1, mean)
  } else {
    meanAbundsLike[,i] <- klustre[, like[[i]]]
  }
}


cor(meanAbundsLike)
write.csv2(round(cor(meanAbundsLike),3), fs::path(scriptPath_UnsupAnalysisUt, "korr_like.csv"))

saveRDS(like, fs::path(scriptPath_UnsupAnalysisUt, "like.rds"))

for(i in 1:length(like)){
  if(length(like[[i]]) > 1){
    print(like[[i]])
  }
}


gates <- read.csv2(fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Data", "gater.csv"))
navn <- read.csv2( fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Data", "navn.csv" ))
navn$marker <- paste0("X", navn$marker)

temp <- read.csv2(fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Panel 1", "All", "seed 1234", "antallPerClusterOgPat_k_10_seed1234All.csv"))
antall_tot_All <- data.frame(id = temp$X, antall_tot = apply(temp[,2:ncol(temp)], 1, sum))
temp <- read.csv2(fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Panel 1", "CD4cells", "seed 4444", "antallPerClusterOgPat_k_5_seed4444CD4cells.csv"))
antall_tot_CD4 <- data.frame(id = temp$X, antall_tot = apply(temp[,2:ncol(temp)], 1, sum))
temp <- read.csv2(fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Panel 1", "CD8cells", "seed 8888", "antallPerClusterOgPat_k_5_seed8888CD8cells.csv"))
antall_tot_CD8 <- data.frame(id = temp$X, antall_tot = apply(temp[,2:ncol(temp)], 1, sum))
temp <- read.csv2(fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "Analyse", "Endelig", "Cytof unsupperviced", "Panel 1", "NKcells", "seed 7444", "antallPerClusterOgPat_k_5_seed7444NKcells.csv"))
antall_tot_NK <- data.frame(id = temp$X, antall_tot = apply(temp[,2:ncol(temp)], 1, sum))

#marker_plot <- function(rader, gates = NULL){


for(i_like in 1:length(like)){
  # i_like <- 23
  Cluster <- like[[i_like]]
  write.table(Cluster, fs::path(scriptPath_UnsupAnalysisUt, paste0("like_klustre_panel1_", i_like, ".txt")))
  if(length(Cluster) > 1){
    params <- list()
    params$panel <- "Panel 1"
    params$d <- klustre[, c("id", Cluster)]
    params$utSti <- fs::path(scriptPath_UnsupAnalysisUt)
    params$antall_tot_all <- antall_tot_All
    params$antall_tot_CD4 <- antall_tot_CD4
    params$antall_tot_CD8 <- antall_tot_CD8
    params$antall_tot_NK <- antall_tot_NK
    params$dInfo <- info
    params$posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", paste0("gating_results_Panel", gsub("Panel ", "", params$panel),"_mars2022"), "posNeg", "Data")
    params$nivaa <- c("Severe T1", "Moderate T1", "Severe T2", "Moderate T2", "Control" )
    rmarkdown::render(fs:::path(scriptPath_UnsupAnalysis , "resultater_Regression_negbin_rapporter_alle.Rmd"), output_file = fs::path(params$utSti,paste0("like_klustre_panel1_", i_like, ".docx")))
    
    
    
    q5 <- data.frame(Cluster = Cluster,  q5_all[Cluster,])
    
    q10 <- data.frame(Cluster = Cluster,  q10_all[Cluster,])
    q25 <- data.frame(Cluster = Cluster,  q25_all[Cluster,])
    q50 <- data.frame(Cluster = Cluster,  q50_all[Cluster,])
    q75 <- data.frame(Cluster = Cluster,  q75_all[Cluster,])
    q90 <- data.frame(Cluster = Cluster,  q90_all[Cluster,])
    q95 <- data.frame(Cluster = Cluster,  q95_all[Cluster,])
    
    
    
    q5long <- data.table::melt(data.table::setDT(q5),id.vars = "Cluster", variable.names = "marker")
    colnames(q5long)[3] <- "q5"
    q10long <- data.table::melt(data.table::setDT(q10),id.vars = "Cluster", variable.names = "marker")
    colnames(q10long)[3] <- "q10"
    q25long <- data.table::melt(data.table::setDT(q25),id.vars = "Cluster", variable.names = "marker")
    colnames(q25long)[3] <- "q25"
    q50long <- data.table::melt(data.table::setDT(q50),id.vars = "Cluster", variable.names = "marker")
    colnames(q50long)[3] <- "q50"
    q75long <- data.table::melt(data.table::setDT(q75),id.vars = "Cluster", variable.names = "marker")
    colnames(q75long)[3] <- "q75"
    q90long <- data.table::melt(data.table::setDT(q90),id.vars = "Cluster", variable.names = "marker")
    colnames(q90long)[3] <- "q90"
    q95long <- data.table::melt(data.table::setDT(q95),id.vars = "Cluster", variable.names = "marker")
    colnames(q95long)[3] <- "q95"
    
    d <-merge(q5long, q10long)
    d <-merge(d, q25long)
    d <-merge(d, q50long)
    d <-merge(d, q75long)
    d <- merge(d, q90long)
    d <- merge(d, q95long)
    d <- merge(d, navn, by.x= "variable", by.y = "marker")
    if(!is.null(gates)){
      d <- merge(d, gates, by = "marker_short_name")
    }
    
    
    # if(!is.null(order_marker_shortname[1])){
    #   if(!is.na(order_marker_shortname[1])){
    #     if(all(d$marker_short_name %in%  order_marker_shortname )){
    #       d$marker_short_name <- factor(d$marker_short_name, levels = order_marker_shortname)
    #     } else {
    #       d$marker_short_name <- factor(d$marker_short_name)
    #       print("could not change order of markers in plot due to some markers not included in list (or misspelled)")
    #     }
    #   } else {
    #     d$marker_short_name <- factor(d$marker_short_name)
    #   }
    # } else {
    d$marker_short_name <- factor(d$marker_short_name)
    # }
    
    d$col <- "0"
    # if(!is.null(highlight_cluster[1])){
    #   if(!is.na(highlight_cluster[1])){
    #     for(i in highlight_cluster){
    #       d$col[d$Cluster == i] <- i
    #     }
    #   }
    # }
    #
    
    unique_col <- unique(d$col)
    
    our_colors <- c("red", c(col50, col50))
    names(our_colors) <- 0:100
    used_colors <- our_colors[as.character(unique_col)]
    
    g <- ggplot(d, aes(x = q50, y = Cluster, xmin = q10, xmax = q90, col = col)) +
      geom_point(size = 2) +
      geom_errorbar() +
      geom_errorbar(data = d, aes(x = q50, y = Cluster, xmin = q5, xmax = q95, col = col))+
      geom_errorbar(data = d, aes(x = q50, y = Cluster, xmin = q25, xmax = q75, col = col))+
      facet_wrap(vars(marker_short_name), ncol = 11, scale = "free_x") +
      scale_color_manual(values = used_colors) 
    
    if(!is.null(gates)){
      g <- g + geom_vline(aes(xintercept  = low), linetype = "dashed", size = 1.1) +
        geom_vline(aes(xintercept  = high), linetype = "dashed", size = 1.1)
    }
    
    tiff( fs::path(scriptPath_UnsupAnalysisUt, paste0("like_klustre_panel1_", i_like, ".tiff")), width = 1350, height = 900)
    print(g)
    dev.off()
    
  }
}




Cluster <- NULL
med <- NULL

for(i in 1:length(like)){
  if(length(like[[i]]) > 1){
    Cluster <- c(Cluster, like[[i]][1])
    med <- c(med, i)
  }
}  




q5 <- data.frame(Cluster = paste0("C", med),  q5_all[Cluster,])
q10 <- data.frame(Cluster = paste0("C", med),  q10_all[Cluster,])
q25 <- data.frame(Cluster = paste0("C", med),  q25_all[Cluster,])
q50 <- data.frame(Cluster = paste0("C", med),  q50_all[Cluster,])
q75 <- data.frame(Cluster = paste0("C", med),  q75_all[Cluster,])
q90 <- data.frame(Cluster = paste0("C", med),  q90_all[Cluster,])
q95 <- data.frame(Cluster = paste0("C", med),  q95_all[Cluster,])



q5long <- data.table::melt(data.table::setDT(q5),id.vars = "Cluster", variable.names = "marker")
colnames(q5long)[3] <- "q5"
q10long <- data.table::melt(data.table::setDT(q10),id.vars = "Cluster", variable.names = "marker")
colnames(q10long)[3] <- "q10"
q25long <- data.table::melt(data.table::setDT(q25),id.vars = "Cluster", variable.names = "marker")
colnames(q25long)[3] <- "q25"
q50long <- data.table::melt(data.table::setDT(q50),id.vars = "Cluster", variable.names = "marker")
colnames(q50long)[3] <- "q50"
q75long <- data.table::melt(data.table::setDT(q75),id.vars = "Cluster", variable.names = "marker")
colnames(q75long)[3] <- "q75"
q90long <- data.table::melt(data.table::setDT(q90),id.vars = "Cluster", variable.names = "marker")
colnames(q90long)[3] <- "q90"
q95long <- data.table::melt(data.table::setDT(q95),id.vars = "Cluster", variable.names = "marker")
colnames(q95long)[3] <- "q95"

d <-merge(q5long, q10long)
d <-merge(d, q25long)
d <-merge(d, q50long)
d <-merge(d, q75long)
d <- merge(d, q90long)
d <- merge(d, q95long)
d <- merge(d, navn, by.x= "variable", by.y = "marker")
if(!is.null(gates)){
  d <- merge(d, gates, by = "marker_short_name")
}


# if(!is.null(order_marker_shortname[1])){
#   if(!is.na(order_marker_shortname[1])){
#     if(all(d$marker_short_name %in%  order_marker_shortname )){
#       d$marker_short_name <- factor(d$marker_short_name, levels = order_marker_shortname)
#     } else {
#       d$marker_short_name <- factor(d$marker_short_name)
#       print("could not change order of markers in plot due to some markers not included in list (or misspelled)")
#     }
#   } else {
#     d$marker_short_name <- factor(d$marker_short_name)
#   }
# } else {
d$marker_short_name <- factor(d$marker_short_name)
# }

for(i in med){
  d$col[d$Cluster %in% paste0("C", i)] <- paste0("C", i)
}


d$Cluster <- factor(d$Cluster, levels = paste0("C", med))
d$col <- factor(d$col, levels = paste0("C", med))

# if(!is.null(highlight_cluster[1])){
#   if(!is.na(highlight_cluster[1])){
#     for(i in highlight_cluster){
#       d$col[d$Cluster == i] <- i
#     }
#   }
# }
#

unique_col <- levels(d$col)

our_colors <- c("red", c(col50, col50))
names(our_colors) <- paste0("C", 0:100)
used_colors <- our_colors[as.character(unique_col)]

g <- ggplot(d, aes(x = q50, y = Cluster, xmin = q10, xmax = q90, col = col)) +
  geom_point(size = 2) +
  geom_errorbar() +
  geom_errorbar(data = d, aes(x = q50, y = Cluster, xmin = q5, xmax = q95, col = col))+
  geom_errorbar(data = d, aes(x = q50, y = Cluster, xmin = q25, xmax = q75, col = col))+
  facet_wrap(vars(marker_short_name), ncol = 11, scale = "free_x") +
  scale_color_manual(values = used_colors)

if(!is.null(gates)){
  g <- g + geom_vline(aes(xintercept  = low), linetype = "dashed", size = 1.1) +
    geom_vline(aes(xintercept  = high), linetype = "dashed", size = 1.1)
}


scriptPath_UnsupAnalysisUt <- fs::path(scriptPath_UnsupAnalysis, "resultat_panel1_cor095")
tiff( fs::path(scriptPath_UnsupAnalysisUt, paste0("like_klustre_panel1_", "_1_fra_hver", ".tiff")), width = 1350, height = 900)
print(g)
dev.off()





#MAIT, like 11, 12, 20, 22, 23, 29



Cluster <- NULL
med <- NULL

for(i in c(11, 12, 20, 22, 23, 29)){
  if(length(like[[i]]) > 1){
    Cluster <- c(Cluster, like[[i]][1])
    med <- c(med, i)
  }
}  




q5 <- data.frame(Cluster = paste0("C", med),  q5_all[Cluster,])
q10 <- data.frame(Cluster = paste0("C", med),  q10_all[Cluster,])
q25 <- data.frame(Cluster = paste0("C", med),  q25_all[Cluster,])
q50 <- data.frame(Cluster = paste0("C", med),  q50_all[Cluster,])
q75 <- data.frame(Cluster = paste0("C", med),  q75_all[Cluster,])
q90 <- data.frame(Cluster = paste0("C", med),  q90_all[Cluster,])
q95 <- data.frame(Cluster = paste0("C", med),  q95_all[Cluster,])



q5long <- data.table::melt(data.table::setDT(q5),id.vars = "Cluster", variable.names = "marker")
colnames(q5long)[3] <- "q5"
q10long <- data.table::melt(data.table::setDT(q10),id.vars = "Cluster", variable.names = "marker")
colnames(q10long)[3] <- "q10"
q25long <- data.table::melt(data.table::setDT(q25),id.vars = "Cluster", variable.names = "marker")
colnames(q25long)[3] <- "q25"
q50long <- data.table::melt(data.table::setDT(q50),id.vars = "Cluster", variable.names = "marker")
colnames(q50long)[3] <- "q50"
q75long <- data.table::melt(data.table::setDT(q75),id.vars = "Cluster", variable.names = "marker")
colnames(q75long)[3] <- "q75"
q90long <- data.table::melt(data.table::setDT(q90),id.vars = "Cluster", variable.names = "marker")
colnames(q90long)[3] <- "q90"
q95long <- data.table::melt(data.table::setDT(q95),id.vars = "Cluster", variable.names = "marker")
colnames(q95long)[3] <- "q95"

d <-merge(q5long, q10long)
d <-merge(d, q25long)
d <-merge(d, q50long)
d <-merge(d, q75long)
d <- merge(d, q90long)
d <- merge(d, q95long)
d <- merge(d, navn, by.x= "variable", by.y = "marker")
if(!is.null(gates)){
  d <- merge(d, gates, by = "marker_short_name")
}


# if(!is.null(order_marker_shortname[1])){
#   if(!is.na(order_marker_shortname[1])){
#     if(all(d$marker_short_name %in%  order_marker_shortname )){
#       d$marker_short_name <- factor(d$marker_short_name, levels = order_marker_shortname)
#     } else {
#       d$marker_short_name <- factor(d$marker_short_name)
#       print("could not change order of markers in plot due to some markers not included in list (or misspelled)")
#     }
#   } else {
#     d$marker_short_name <- factor(d$marker_short_name)
#   }
# } else {
d$marker_short_name <- factor(d$marker_short_name)
# }

for(i in med){
  d$col[d$Cluster %in% paste0("C", i)] <- paste0("C", i)
}


d$Cluster <- factor(d$Cluster, levels = paste0("C", med))
d$col <- factor(d$col, levels = paste0("C", med))

# if(!is.null(highlight_cluster[1])){
#   if(!is.na(highlight_cluster[1])){
#     for(i in highlight_cluster){
#       d$col[d$Cluster == i] <- i
#     }
#   }
# }
#

unique_col <- levels(d$col)

our_colors <- c("red", c(col50, col50))
names(our_colors) <- paste0("C", 0:100)
used_colors <- our_colors[as.character(unique_col)]

g <- ggplot(d, aes(x = q50, y = Cluster, xmin = q10, xmax = q90, col = col)) +
  geom_point(size = 2) +
  geom_errorbar() +
  geom_errorbar(data = d, aes(x = q50, y = Cluster, xmin = q5, xmax = q95, col = col))+
  geom_errorbar(data = d, aes(x = q50, y = Cluster, xmin = q25, xmax = q75, col = col))+
  facet_wrap(vars(marker_short_name), ncol = 11, scale = "free_x") +
  scale_color_manual(values = used_colors)

if(!is.null(gates)){
  g <- g + geom_vline(aes(xintercept  = low), linetype = "dashed", size = 1.1) +
    geom_vline(aes(xintercept  = high), linetype = "dashed", size = 1.1)
}


scriptPath_UnsupAnalysisUt <- fs::path(scriptPath_UnsupAnalysis, "resultat_panel1_cor095")
tiff( fs::path(scriptPath_UnsupAnalysisUt, paste0("like_klustre_panel1_", "_1_fra_hver_MAIT", ".tiff")), width = 1350, height = 900)
print(g)
dev.off()






