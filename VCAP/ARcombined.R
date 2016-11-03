# Check if a gene that is AR regulated but changes late in time is activated by a combination of AR with another TF

# I check the neighbourhood of the AR peaks and of the TSS of the gene under consideration for motifs of other TFs
#   I need to explore all the TFs not only those that change in the datasets. So check the sequences with all the 
#   motifs I have from all the databases.

# 

genes_corr_illumina_position <- androgen_probes_timechange[androgen_probes_timechange$EnsemblReannotated%in%rcade_results$geneID,]
dim(genes_corr_illumina_position)
genes_corr_position <- ordered_tfs_position[ordered_tfs_position$Ensembl.Gene.ID%in%rcade_results$geneID,]
dim(genes_corr_position)
#Get the IlluminaID
VCapAndrogen_passed_mean_illumina <- sapply(strsplit(rownames(VCapAndrogen_passed_mean),"[.]"),function(i) i[[1]])
#Select those that are in the Rcade output (Genes activated by AR)
VCapAndrogen_passed_mean_Rcade <- VCapAndrogen_passed_mean[which(VCapAndrogen_passed_mean_illumina%in%genes_corr_illumina_position$IlluminaID),]

#Percentage matrix
VCapAndrogen_pm_Rcade_percentage <- percentage_matrix(VCapAndrogen_passed_mean_Rcade)

#KML

library(kml)
kml_VCapAndrogen_pm_Rcade_percentage <- cld(VCapAndrogen_pm_Rcade_percentage)
parWithEuclidean <- parALGO(distanceName="euclidean")
kml(kml_VCapAndrogen_pm_Rcade_percentage,nbClusters=13:15,nbRedrawing=10,toPlot='none',parAlgo=parWithEuclidean)

#commands to choose the partition that best suits me.
# It will return popup window which we can use to select the options we like.
X11(type="Xlib")
choice(kml_VCapAndrogen_pm_Rcade_percentage, typeGraph = "bmp")
# use the arrows to see the partitions - see the menu for other changes
# use the space bar to select the partition I want
# use 'm' to export al the result in the current directory. cvs and bmp files

cc<-15
dir.create(paste("ARcombined_KML_",cc,sep=""))
system(paste0("mv kml* ARcombined_KML_",cc, sep=""))

cluster_file<-read.table(paste(getwd(),"/ARcombined_KML_",cc,"/kml_VCapAndrogen_pm_Rcade_percentage-C",cc,"-1-Clusters.csv",sep=""),header=TRUE, sep=";")

cluster_indeces <- lapply(1:cc, function(i) which(cluster_file[,2]==i))

#Plot all genes in all clusters.
tiff(paste(getwd(),"/ARcombined_KML_",cc,"/all_",cc,"_clusters.tiff",sep=""), width=3000, height=2000)
Col <- ceiling(length(cluster_indeces)/4)
par(mfrow=c(4,Col))
for(cl in 1:length(cluster_indeces)){
  plot(time_Androgen,VCapAndrogen_pm_Rcade_percentage[cluster_indeces[[cl]][1],],type="l", 
       main=paste("Cluster", cl, sep=" "),
       cex.main=3, xaxt='n', xlab='', ylab='', cex.axis=2, ylim=c(0,1))
  for(i in cluster_indeces[[cl]]){  
    lines(time_Androgen,VCapAndrogen_pm_Rcade_percentage[i,], col=i)
  }
  axis(1, at=time_Androgen, lab=time_Androgen,cex.axis=2)    
}
dev.off()


# Get the sequences for each cluster ----

#Get the corresponding IlluminaID
VCapAndrogen_pm_Rcade_illumina <- sapply(strsplit(rownames(VCapAndrogen_passed_mean_Rcade),"[.]"),function(i) i[[1]])
cluster_VCapAndrogen_pm_Rcade_illumina <- sapply(cluster_indeces, function(i) VCapAndrogen_pm_Rcade_illumina[i])

#Get the correspondence of IluminaID and EnsemblID with the location range
genes_corr_illumina_position_indeces <- sapply(cluster_VCapAndrogen_pm_Rcade_illumina, function(i) 
                                                which(genes_corr_illumina_position$IlluminaID%in%i))

ARcombined_allinfo_cluster <- lapply(genes_corr_illumina_position_indeces, function(i) genes_corr_illumina_position[i,])
ARcombined_cluster_locations <- lapply(genes_corr_illumina_position_indeces, function(i) genes_corr_position[i,])

# As I need to explore the neighbourhood of both AR and the gene under consideration, I take now the locations of AR
# found in ARmotif_analysis using FIMO. 
# For each gene I have multiple AR location. What I can do is take the min and max of all locations (they should all 
# be next to each other) and then get the neighbourhood.


ARcombined_allnames_cluster <- lapply(ARcombined_allinfo_cluster, function(i) paste(i$EnsemblReannotated,
                                                                                    i$IlluminaID,
                                                                                    i$GenomicLocation,sep="."))
corresponding_ARlocations <- lapply(ARcombined_allnames_cluster, function(i) 
  corresponding_all_location[which(corresponding_all_location$sequence.name%in%i),])

for(i in 1:length(corresponding_ARlocations)){
  temp_gene <- NULL
  temp_gene <- as.character(corresponding_ARlocations[[i]]$gene_strand)
  temp_gene[which(temp_gene=="1")] <- "+"
  temp_gene[which(temp_gene=="-1")] <- "-"
  corresponding_ARlocations[[i]]$gene_strand <- temp_gene
}

# AR NEIGHBOURHOOD ----
ar_start <- lapply(corresponding_ARlocations, function(i) i$gene_start-25000+i$ARstart)
ar_end <- lapply(corresponding_ARlocations, function(i) i$gene_start-25000+i$ARstop)#it counts from the beginning

ar_neighbourhood <- lapply(corresponding_ARlocations, function(i)
                            as.data.frame(cbind(as.character(i$sequence.name),
                                                as.character(i$gene_chr),
                                                i$gene_start-25000+i$ARstart,
                                                i$gene_start-25000+i$ARstop,
                                                as.character(i$gene_strand))))
                                                # as.character(i$ARstrand))))
# I use the same strand as the gene because I am concatenating the two dna-strings later on.
rm(ar_start);rm(ar_end)

ARcombined_ar_neighbourhood <- lapply(ar_neighbourhood, function(i) rename(i,c("V1"="SequenceName","V2"="chr","V3"="start","V4"="end","V5"="str")))

# GENES NEIGHBOURHOOD ----
ARcombined_genes_neighbourhood <- lapply(corresponding_ARlocations, function(i) i[,c("sequence.name","gene_chr","gene_start","gene_end","gene_strand")])
ARcombined_genes_neighbourhood <- lapply(ARcombined_genes_neighbourhood, function(i) 
                                          rename(i,c("sequence.name"="SequenceName",
                                                     "gene_chr"="chr",
                                                     "gene_start"="start",
                                                     "gene_end"="end",
                                                     "gene_strand"="str")))


#Concatenate the gene and AR location strings ----
ARcombined_genes_sequences <- ARcombined_ar_sequences <- vector("list", length(ARcombined_genes_neighbourhood))
for(i in 1:length(ARcombined_ar_neighbourhood)){
  ARcombined_genes_sequences[[i]] <- get_sequence(ARcombined_genes_neighbourhood[[i]], 
                                  1500,"NO") 
                                  # ,paste("ARcombined_genes_sequences_cluster",i,".fasta",sep=""))
  ARcombined_ar_sequences[[i]] <- get_sequence(ARcombined_ar_neighbourhood[[i]], 
                               1500,"NO")
                               #paste("ARcombined_ar_sequences_cluster",i,".fasta",sep=""))
  print(i)
}


#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞
# Combine the two sets of sequences ----
#  The tf can bind either around AR peak or around the TSS of the gene under consideration
ARcombined_sequences <- vector("list", length(ARcombined_genes_sequences))
for(i in 1:length(ARcombined_genes_sequences)){
  ARcombined_sequences[[i]] <- vector("list", length(ARcombined_genes_sequences[[i]]))
  for(j in 1:length(ARcombined_genes_sequences[[i]])){
    ARcombined_sequences[[i]][[j]] <- DNAString(paste0(as.character(ARcombined_genes_sequences[[i]][[j]]),
                                                       as.character(ARcombined_ar_sequences[[i]][[j]])))
    
    write.table(paste(">",ARcombined_genes_neighbourhood[[i]]$SequenceName[j],sep=""), 
                file=paste("ARcombined_sequences_cluster",i,".fasta",sep=""), 
                append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
    write.table(toString(ARcombined_sequences[[i]][[j]]), 
                file=paste("ARcombined_sequences_cluster",i,".fasta",sep=""),
                append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
    
  }
}
  
# ARcombined_sequences_set <- lapply(ARcombined_sequences, function(i) DNAStringSet(i))


  























#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞---FUNCTIONS---∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞

# Get the sequences from the UCSC ----

get_sequence <- function(Data,Range, TitleFasta){
  Data <- ordered_tfs_position
  Range<-1500
  TitleFasta <- "/Volumes/Valeria/CRUK/AR_downstream_analysis/VCAP/background_sequences/background_sequences.out"
  Data_range <- Data
  Data_range$chr <- paste("chr",Data_range$chr, sep="")
  Data_range$start <- as.numeric(as.character(Data$start))-Range
  Data_range$end <- as.numeric(as.character(Data$end))+Range
  Data_range$str <- as.character(Data_range$str)
  # Load the database
  library("BSgenome.Hsapiens.UCSC.hg38")
  ls(1)
  sequence_results <- vector("list", nrow(Data_range))
  for(i in 1:nrow(Data_range)){
    sequence_results[[i]] <- getSeq(Hsapiens,
                                    Data_range$chr[i],
                                    Data_range$start[i], Data_range$end[i],
                                    width=NA,
                                    Data_range$str[i])
    # Write the sequences in a FASTA format file 
    if(TitleFasta!="NO"){
      write.table(paste(">",Data$SequenceName[i],sep=""),
                  file=TitleFasta, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
      write.table(toString(sequence_results[[i]]), file=TitleFasta,
                  append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
      
      print(paste(i, nchar(toString(sequence_results[[i]])),sep="-"))
    }
  }
  return(sequence_results)
}

#------------------------------------------------------------------------------------------------------------------

# Percentage matrix for KML ----

percentage_matrix <- function(Matrix){
  perc_matrix <- matrix(0,nrow(Matrix), ncol(Matrix))
  for(i in 1:nrow(Matrix)){
    perc_matrix[i,] <- as.matrix((Matrix[i,]-min(Matrix[i,]))/(max(Matrix[i,])-min(Matrix[i,])))
  }
  
  rownames(perc_matrix) <- rownames(Matrix)
  colnames(perc_matrix) <- colnames(Matrix)
  return(perc_matrix)
}

#------------------------------------------------------------------------------------------------------------------


