
#Manually correct the NAs ----

# dir.create("Figures")

index_na <- which(is.na(point_change))
tiff(paste(getwd(),"/Figures/point_change.tiff",sep=""), width=2000, height=1000)
par(mfrow=c(5,7))
for(i in index_na){
  plot(time_Androgen, VCapAndrogen_passed_mean[i,], type="l", main=i)
  lines(time_Androgen, VCapAndrogen_passed_mean_loess[i,], col="blue")
  lines(time_Control, VCapControl_passed_mean[i,], col="red")
}
dev.off()

index_androgen_na <- which(is.na(androgen_point_change))
tiff(paste(getwd(),"/Figures/androgen_point_change.tiff",sep=""), width=2000, height=1000)
par(mfrow=c(5,4))
for(i in index_androgen_na){
  plot(time_Androgen, VCapAndrogen_passed_mean[i,], type="l", main=i)
  lines(time_Androgen, VCapAndrogen_passed_mean_loess[i,], col="blue")
  # lines(time_Control, VCapControl_passed_mean[i,], col="red")
}
dev.off()

# ind_na <- index_na[c(5,11)]
# tiff(paste(getwd(),"/Figures/point_change2.tiff",sep=""), width=2000, height=1000)
# par(mfrow=c(2,2))
# for(i in ind_na){
#   plot(time_VCapAndrogen, VCapAndrogen_passed_mean[i,], type="l", main=i)
#   lines(time_VCapAndrogen, VCapAndrogen2[i,], col="blue")
#   lines(time_VCapControl, VCapControl_passed_mean_loess[i,], col="red")
# }
# dev.off()
# 
# substitute_na <- rep(NA, length(index_na))
# substitute_na[5] <- 5
# substitute_na[11] <- 9
# 
# point_change[index_na] <- substitute_na










tfs_time_array <- sapply(tfs_time, function(i) paste(i, collapse=","))

for(k in 1:length(tfs_time)){
  
  tiff(paste(getwd(),"/Figures/acf/time_change",time_change[k],".tiff",sep=""), width=2000, height=1000)
  n_row <- ceiling(length(tfs_time[[k]])/6)
  par(mfrow=c(6,n_row))
  
  for(i in 1:length(tfs_time[[k]])){
    a<-which(rownames(VCapAndrogen_passed_mean)==tfs_time[[k]][i])
    
    plot(time_Androgen,VCapAndrogen_passed_mean[a,],type="l", 
         main=tfs_time[[k]][i],
         cex.main=3, xaxt='n', xlab='', ylab='', cex.axis=2,
         # ylim=c(7,15))
         ylim=c(min(c(as.numeric(VCapAndrogen_passed_mean[a,]),as.numeric(VCapControl_passed_mean[a,]))),
                max(c(as.numeric(VCapAndrogen_passed_mean[a,]),as.numeric(VCapControl_passed_mean[a,])))))
    lines(time_Androgen,VCapAndrogen_passed_mean_loess[a,])
    lines(time_Androgen,VCapControl_passed_mean_loess[a,], col="red")
    lines(time_Control,VCapControl_passed_mean[a,], col="red")
    abline(v=time_change[k])
    axis(1, at=time_Androgen, lab=time_Androgen,cex.axis=2)  
  }
  dev.off()
}

tf_point <- cbind(rownames(VCapAndrogen_passed_mean), point_change)

# Which one decreases the most 
# max- min = highest
stripchart(time_change, vertical=T) #method="jitter"
text(rep(1,length(time_change)),
     time_change,tfs_time_array,pos=4, font=1)
















# Point change calculated using ONLY androgen ----

androgen_time_change <- sort(unique(androgen_point_change))
# Indeces of the TFs that change at each time point
androgen_ind_change <- sapply(androgen_time_change, function(i) which(androgen_point_change==i), simplify=FALSE)
# Symbol_ProbeID of the TFs that change at each time point
androgen_tfs_time <- lapply(androgen_ind_change, function(i) rownames(VCapAndrogen_passed_mean)[i])

#Write the list of genes that change in a text file
for(i in 1:length(androgen_tfs_time)){
  write.table(androgen_tfs_time[[i]], file=paste("androgen_TFs_change_", androgen_time_change[i], ".txt", sep=""),
              quote=FALSE, col.names=FALSE, row.names=FALSE)
}


androgen_tfs_time_array <- sapply(androgen_tfs_time, function(i) paste(i, collapse=","))

for(k in 2:length(androgen_tfs_time)){
  
  tiff(paste(getwd(),"/Figures/acf/androgen_time_change",androgen_time_change[k],".tiff",sep=""), width=2000, height=1000)
  n_row <- ceiling(length(androgen_tfs_time[[k]])/6)
  par(mfrow=c(6,n_row))
  
  for(i in 1:length(androgen_tfs_time[[k]])){
    a<-which(rownames(VCapAndrogen_passed_mean)==androgen_tfs_time[[k]][i])
    
    plot(time_Androgen,VCapAndrogen_passed_mean[a,],type="l", 
         main=androgen_tfs_time[[k]][i],
         cex.main=3, xaxt='n', xlab='', ylab='', cex.axis=2,
         # ylim=c(7,15))
         ylim=c(min(c(as.numeric(VCapAndrogen_passed_mean[a,]),as.numeric(VCapControl_passed_mean[a,]))),
                max(c(as.numeric(VCapAndrogen_passed_mean[a,]),as.numeric(VCapControl_passed_mean[a,])))))
    lines(time_Androgen,VCapAndrogen_passed_mean_loess[a,])
    lines(time_Androgen,VCapControl_passed_mean_loess[a,], col="red")
    lines(time_Control,VCapControl_passed_mean[a,], col="red")
    abline(v=androgen_time_change[k])
    axis(1, at=time_Androgen, lab=time_Androgen,cex.axis=2)  
  }
  dev.off()
}

androgen_tf_point <- cbind(rownames(VCapAndrogen_passed_mean), androgen_point_change)

















#LIMMA and RCade
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 2- Get the corresponding names ----

# "Go to Biomart -> Choose DB -> choose ds -> Homo sapiens -> Gene type -> Protein_coding (protein_coding gene) -> 
# Attributes -> Gene -> Select the attributes you are interested in.
# The click on Count - should get about 20,000
# and then finally on results and you'll sve it."
# gene_ensembl <- read.table("mart_export.txt", header=TRUE, sep="\t")
# 
# tfs_ensembl_position <- NULL
# for(i in 1:nrow(ordered_tfs)){
#   tfs_ensembl_position <- rbind(tfs_ensembl_position, 
#                                 gene_ensembl[which(as.character(gene_ensembl[,1])==as.character(ordered_tfs[i,1])),-6])
# }
# 
# which(!ordered_tfs[,1]%in%tfs_ensembl_position[,1]) #270
# # ordered_tfs[270,]
# # ensembl_gene_id hgnc_symbol entrezgene chromosome_name start_position end_position strand illumina_humanwg_6_v2
# # 2376 ENSG00000242779     ZNF702P      79986              19       52968251     53037898     -1          ILMN_1663281
# # ZNF702P is a pseudogene so I'll ignore it.
# 
# ordered_tfs2 <- ordered_tfs[-270,] 
# 
# 
# all.equal(as.character(tfs_ensembl_position[,1]),as.character(ordered_tfs2[,1]))
# all.equal(as.character(tfs_ensembl_position[,3]),as.character(ordered_tfs2[,5]))
# all.equal(as.character(tfs_ensembl_position[,4]),as.character(ordered_tfs2[,6]))
# # Ensembl, start and end are the equivalent. I can directly use ordered_tfs2. 
# # I used gene_ensembl to spot gene sthat are not coding, such as pseudogenes.









# I am only interested in AR motifs

# library(JASPAR2016)
# library(TFBSTools)
# 
# opts <- list()
# # opts[["species"]] <- 9606 #human 10090 #mouse
# opts[["tax_group"]] <- "vertebrates"
# opts[["collection"]] <- "CORE"
# opts[["matrixtype"]] <- "ICM"
# PFMatrixList_pfm <- getMatrixSet(JASPAR2016, opts)
# 
# pfm_jaspar <- Matrix(PFMatrixList_pfm)
# names(pfm_jaspar) <- toupper(name(PFMatrixList_pfm))
# 
# ARmotif <- pfm_jaspar[[which(names(pfm_jaspar)=="AR")]]
# 
# write.table(">AR", file="ARmotif.pfm", col.names=FALSE, row.names=FALSE)
# write.table(ARmotif, file="ARmotif.pfm", append=TRUE, col.names=FALSE)

# jaspar2meme </Volumes/Valeria/CRUK/AR_downstream_analysis/VCAP/ARmotif.txt | /Volumes/Valeria/CRUK/AR_downstream_analysis/VCAP>

#Scan sequences using the function FIMO from the MEMEsuite ----








# library(biomaRt)
# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# filters=listFilters(ensembl)
# Seq_gene_exon_upstream[i] <- getSequence(chromosome=as.integer(genes_corr_position$chr[i]), 
#                                         start=genes_corr_position$start[i], end=genes_corr_position$end[i],
#                                         seqType="gene_flank", type= "ensembl_transcript_id",
#                                         upstream=25000, mart=ensembl, verbose=FALSE)$gene_flank

