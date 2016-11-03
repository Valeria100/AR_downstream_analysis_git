#LNCAP ----


# AUTOCORRELATION --------------------------------------------------------------------------------------------------------------

# Select the genes that have a behavioural change that differs from natural fluctuation.
# Autocorrelation higher than 0 (a threshold) indicate a change in behaviour possibly caused
# by AR 

#I calculate autocorrelation on each dataset and each repetition and then select the genes that pass the autocorrelation filter 
# in both repetitions

calculate_acf <- function(Data,ci){
  acf_data <- NULL
  for(i in 1:nrow(Data)){
    if(all(is.na(Data[i,]))){
      acf_data[i] <- NA
    }else{
      acf_data[i] <- acf(as.numeric(Data[i,]),type="correlation",lag.max=nrow(Data), plot=FALSE)$acf[2]
    }
  }
  # Only the positive acf gets selected. when acf<0 is because the signal goes up and down.
  acf_index <- which(acf_data>=ci) #which(abs(acf_data)>=ci) #
  acf_value <- acf_data[acf_index]
  acf_Data <- cbind(acf_index,acf_value)
  return(acf_Data)
}
# the qnorm... is how acf calculates the confidence intervals you can see in the plots.
# I use it as threshold.
LNA1 <- calculate_acf(LNCapAndrogen1,qnorm((1+0.95)/2)/sqrt(ncol(LNCapAndrogen1)))
LNA2 <- calculate_acf(LNCapAndrogen2,qnorm((1+0.95)/2)/sqrt(ncol(LNCapAndrogen2)))

LNC1 <- calculate_acf(LNCapControl1,qnorm((1+0.95)/2)/sqrt(ncol(LNCapControl1)))
LNC2 <- calculate_acf(LNCapControl2,qnorm((1+0.95)/2)/sqrt(ncol(LNCapControl2)))

#This needs to be less stringent as LNCap variates much less
index_LNCapAndrogen <- sort(union(LNA1[,1],LNA2[,1]))

# probes_acf_LNCap <- probes[index_LNCapAndrogen,]

# VARIANCE AND F-TEST - COMPARE VARIANCES OF TWO POPULATIONS -------------------------------------------------------------------

# Calculate variance for both androgens and controls. 
apply_ftest <- function(a1, a2, c1, c2, ind){
  var11 <- var12 <- var21 <- var22 <- NULL
  for(i in 1:nrow(a1)){
    var11[i] <- var.test(as.numeric(a1[i,]),as.numeric(c1[i,]))$p.value
    var12[i] <- var.test(as.numeric(a1[i,]),as.numeric(c2[i,]))$p.value
    var21[i] <- var.test(as.numeric(a2[i,]),as.numeric(c1[i,]))$p.value
    var22[i] <- var.test(as.numeric(a2[i,]),as.numeric(c2[i,]))$p.value  
  }
  
  #ind allow to refer the indeces to the original dataset
  a <- ind[which(var11<=0.05)]
  b <- ind[which(var12<=0.05)]
  d <- ind[which(var21<=0.05)]
  e <- ind[which(var22<=0.05)]
  
  ind_passed <- sort(unique(c(a,b,d,e)))
  return(ind_passed)
}

#indeces are still based on the original datasets
# ind_acf_ftest_VCap <- apply_ftest(VCapAndrogen1[index_VCap_Androgen_notControl,], VCapAndrogen2[index_VCap_Androgen_notControl,],
#                                   VCapControl1[index_VCap_Androgen_notControl,], VCapControl2[index_VCap_Androgen_notControl,],
#                                   index_VCap_Androgen_notControl)

#Not for LNCap as this doesn't variates enough to measure the difference.

LNCapAndrogen1_genesTF_passed <- LNCapAndrogen1[index_LNCapAndrogen,]
LNCapAndrogen2_genesTF_passed <- LNCapAndrogen2[index_LNCapAndrogen,]
LNCapControl1_genesTF_passed <- LNCapControl1[index_LNCapAndrogen,]
LNCapControl2_genesTF_passed <- LNCapControl2[index_LNCapAndrogen,]

#from the probe file - contains ProbeType, VCapProbID, VCapSymbol, etc..
LNCapProbes_genesTF_passed <- probes[index_LNCapAndrogen,]

#"hgnc_symbol-illumina_probe"
rownames(LNCapAndrogen1_genesTF_passed) <- rownames(LNCapAndrogen2_genesTF_passed) <- probes_LNCap[index_LNCapAndrogen]


#----------------------------------------------------------------------------------------------------------------------------------
# SELECT ONLY THE TFs -------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------


#Check if any of the genes selcted until now are transcription factors
# Script of this is in /media/bo01/     or    /Volumes/Valeria/CRUK/Human_TFs
# tfs <- read.table("/Volumes/Valeria/CRUK/Human_TFs/unique_TFs_symbol.txt", header=TRUE)
load("/Volumes/Valeria/CRUK/Human_TFs/ExpAnnotation_VCap_TFs.RData")
tfs_lncap <- get(load("/Volumes/Valeria/CRUK/Human_TFs/ExpAnnotation_LNCap_TFs.RData"))[,c(1,11,12,13,16)]
#no need for toupper as they already are all upper case.

ind_tfs_LNCap <- which((as.character(probes[index_LNCapAndrogen,5])%in%as.character(tfs_lncap$SymbolReannotated))|
                         (as.character(probes[index_LNCapAndrogen,2])%in%as.character(tfs_lncap$IlluminaID))|
                         (as.character(probes[index_LNCapAndrogen,4])%in%as.character(tfs_lncap$IlluminaID)))

#These are the indeces based on the original datasets.
ind_acf_tfs_LNCap <- index_LNCapAndrogen[ind_tfs_LNCap]

# ind_common <- ind_acf_tfs_LNCap[ind_acf_tfs_LNCap%in%ind_acf_tfs_VCap]

LNCapAndrogen1_passed <- LNCapAndrogen1[ind_acf_tfs_LNCap,]
LNCapAndrogen2_passed <- LNCapAndrogen2[ind_acf_tfs_LNCap,]
LNCapControl1_passed <- LNCapControl1[ind_acf_tfs_LNCap,]
LNCapControl2_passed <- LNCapControl2[ind_acf_tfs_LNCap,]

#from the probe file - contains ProbeType, VCapProbID, VCapSymbol, etc..
LNCapProbes_passed <- probes[ind_acf_tfs_LNCap,]

#"hgnc_symbol-illumina_probe"
rownames(LNCapAndrogen1_passed) <- rownames(LNCapAndrogen2_passed) <- probes_LNCap[ind_acf_tfs_LNCap]

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Ensembl_id_correspondence

tfs_lncap_in <- tfs_lncap[tfs_lncap$IlluminaID%in%LNCapProbes_passed$LNCap.Probe_id,]
tfs_lncap_ord <- tfs_lncap_in[order(match(tfs_lncap_in$IlluminaID,LNCapProbes_passed$LNCap.Probe_id)),]
rm(tfs_lncap_in)

ordered_tfs <- tfs_lncap_ord















#CHECKED UNTIL HERE! THE REST NEED TO BE CHANGED!!!!!!!!!!




























#----------------------------------------------------------------------------------------------------------------------------------
# TIME POINT OF CHANGE ------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------

# Calculate where the behaviour of the gene change


#The repetitions have different time points therefore I cannot calculate the mean.
# I need to operate on each matrix separately.
# Calculate slope (y2-y1)/(x2-x1)

subtract_columns_matrix <- function(Matrix){  
  Matrix_sub <- matrix(0, nrow(Matrix), (ncol(Matrix)-1))  
  for(i in 1:(ncol(Matrix)-1)){  
    Matrix_sub[,i] <- Matrix[,i+1]-Matrix[,i]
  }
  return(Matrix_sub)  
}

subtract_elements_array <- function(Array){  
  Array_sub <- NULL  
  for(i in 1:(length(Array)-1)){  
    Array_sub[i] <- Array[i+1]-Array[i]
  }
  return(Array_sub)
}

divide_row_array <- function(Matrix, Array){
  Matrix_div <- matrix(0, nrow(Matrix), ncol(Matrix))    
  for(i in 1:nrow(Matrix)){
    Matrix_div[i,] <- round(Matrix[i,]/Array,digits=1)
  }
  return(Matrix_div)
}

mean_each2 <- function(Array){
  Array_mean <- NULL
  for(i in 1:(length(Array)-1)){
    Array_mean[i] <- round(mean(c(Array[i], Array[i+1])),digits=2)
  }
  return(Array_mean)
}

time_Androgen <- time_VCapAndrogen1#[which(time_VCapAndrogen1%in%time_VCapAndrogen2)]
time_Control <- intersect(time_VCapControl1,time_VCapControl2)


# MEAN of the repetitions ----

VCapAndrogen_passed_mean <- VCapAndrogen1_passed 
#(VCapAndrogen1_passed[which(time_VCapAndrogen1%in%time_Androgen)]+VCapAndrogen2_passed[which(time_VCapAndrogen2%in%time_Androgen)])/2
VCapControl_passed_mean <- (VCapControl1_passed[which(time_VCapControl1%in%time_Control)]+
                              VCapControl2_passed[which(time_VCapControl2%in%time_Control)])/2

#LOESS of the mean of the repetitions -------------------------------------------------------------------------------------------
# I extrapolate the points of the controls using the time points of the androgen to not loose all the time information

VCapAndrogen_passed_mean_loess <- VCapControl_passed_mean_loess <- matrix(0, nrow(VCapAndrogen_passed_mean), 
                                                                          ncol(VCapAndrogen_passed_mean))
for(i in 1:nrow(VCapAndrogen_passed_mean)){
  VCapAndrogen_passed_mean_loess[i,] <- predict(loess(as.numeric(VCapAndrogen_passed_mean[i,])~time_Androgen, 
                                                      span=0.5))
  VCapControl_passed_mean_loess[i,] <- predict(loess(as.numeric(VCapControl_passed_mean[i,])~time_Control, 
                                                     span=0.5), time_Androgen)
}

# VCapAndrogen_passed_mean_loess_shifted <- cbind(VCapAndrogen_passed_mean_loess[,-1],
#                                                  VCapAndrogen_passed_mean_loess[,ncol(VCapAndrogen_passed_mean_loess)])
VCapAndrogen_passed_mean_shifted <- cbind(VCapAndrogen_passed_mean[,-1],
                                          VCapAndrogen_passed_mean[,ncol(VCapAndrogen_passed_mean)])

# Calculate the SLOPE
# calculate_slope <- function(Matrix, Time){
#   mat_sub <- subtract_columns_matrix(Matrix)
#   time_sub <- subtract_elements_array(Time)
#   matrix_slope <- divide_row_array(mat_sub, time_sub)
#   return(matrix_slope)
# }
# VCapAndrogen_passed_mean_slope <- calculate_slope(VCapAndrogen_passed_mean, time_Androgen)


#-------------------------------------------------------------------------------------------------------------------------------
library(TTR)

androgen_point_change <- point_change <- NULL
androgen_range_change <- range_change <- matrix(0,nrow(VCapAndrogen_passed_mean_loess),2)

thr_ccf <- round(qnorm((1+0.95)/2)/sqrt(ncol(VCapAndrogen_passed_mean_loess)), digits=3)

cc <- ac <- vector("list", nrow(VCapAndrogen_passed_mean_loess))
for(i in 1:nrow(VCapAndrogen_passed_mean_loess)){
  
  cc[[i]] <- ccf(as.numeric(VCapAndrogen_passed_mean_loess[i,]),
                 as.numeric(VCapControl_passed_mean_loess[i,]),
                 lag.max=length(as.numeric(VCapAndrogen_passed_mean_loess[i,])), 
                 type="correlation", plot=FALSE)
  
  # ac[[i]] <- acf(as.numeric(VCapAndrogen_passed_mean_loess[i,]),
  #           lag.max=length(as.numeric(VCapAndrogen_passed_mean_loess[i,])),  
  #           type="correlation", plot=FALSE)
  
  crosscorrelation <- round(cc[[i]]$acf,digits=3)
  
  # autocorrelation <- round(ac[[i]]$acf, digits=3)
  
  Lag <- cc[[i]]$lag
  
  above_thr <- which(abs(crosscorrelation) > 0)
  point_change_ccf <- abs(Lag[above_thr])
  
  
  ##########
  
  thr_ma <- which(abs(VCapAndrogen_passed_mean[i,]-VCapControl_passed_mean_loess[i,])>0.4)
  point_change_ma <- time_Androgen[thr_ma]
  
  maA <- round(SMA(as.numeric(VCapAndrogen_passed_mean_loess[i,]),5),2)
  # maC <- round(SMA(as.numeric(VCapControl_passed_mean_loess[i,]),5),2)
  
  androgen_point_change[i] <- which(round(abs(maA - as.numeric(VCapAndrogen_passed_mean_loess[i,])),digits=2) > 0.1)[1]
  androgen_range_change[i,1] <- androgen_point_change[i]-1
  androgen_range_change[i,2] <- androgen_point_change[i]+1
  
  ##########
  #NA if no time points intersect
  point_change[i] <- sort(intersect(point_change_ccf, point_change_ma-2))[1]
  if(!is.na(point_change[i])&(point_change[i]==0)){
    range_change[i,1] <- point_change[i]
    range_change[i,2] <- point_change[i]+2
  }else{
    range_change[i,1] <- point_change[i]-1
    range_change[i,2] <- point_change[i]+1
  }
  
}

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


# Point change calculated using both androgen and control ----
time_change <- sort(unique(point_change))
# Indeces of the TFs that change at each time point
ind_change <- sapply(time_change, function(i) which(point_change==i), simplify=FALSE)
# Symbol_ProbeID of the TFs that change at each time point
tfs_time <- lapply(ind_change, function(i) rownames(VCapAndrogen_passed_mean)[i])

#Write the list of genes that change in a text file
for(i in 1:length(tfs_time)){
  write.table(tfs_time[[i]], file=paste("TFs_change_", time_change[i], ".txt", sep=""),
              quote=FALSE, col.names=FALSE, row.names=FALSE)
}


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


#Compare the Point change of the two approaches ----

difference_point_change <- point_change - androgen_point_change

hist(difference_point_change)
table(difference_point_change)


#----------------------------------------------------------------------------------------------------------------------------------
# Integrate CHIPSEQ DATA - apply_rcade --------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------
# apply_rcade <- function(Data, Probes_passed, Probes_passed_names){
library(Rcade)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1- Calculate differential expression, fold change (microarray data) ----

library(limma)

Time <- c(time_VCapAndrogen1, time_VCapAndrogen2, time_VCapControl1, time_VCapControl2)
Group <- factor(c(rep("Androgen", length(c(time_VCapAndrogen1, time_VCapAndrogen2))),
                  rep("Control", length(c(time_VCapControl1, time_VCapControl2)))))
FileName <- c(paste("VA1_",time_VCapAndrogen1,sep=""),
              paste("VA2_",time_VCapAndrogen2,sep=""),
              paste("VC1_",time_VCapControl1,sep=""),
              paste("VC2_",time_VCapControl2,sep=""))


merged_VCap <- cbind(VCapAndrogen1_passed, VCapAndrogen2_passed, VCapControl1_passed, VCapControl2_passed)
colnames(merged_VCap) <- FileName

# both TFs and Genes ---

# merged_VCap <- cbind(VCapAndrogen1_genesTF_passed, VCapAndrogen2_genesTF_passed, VCapControl1_genesTF_passed, VCapControl2_genesTF_passed)
# colnames(merged_VCap) <- FileName

library(splines)
X <- ns(Time, df=5)

design <- model.matrix(~Group*X)

fit <- lmFit(merged_VCap,design)

fit <- eBayes(fit)

TopTable <- topTable(fit, adjust="BH", n=nrow(merged_VCap))

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


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 3- Create the DE Matrix and the DE lookup table ----

DE <- cbind(ordered_tfs$EnsemblReannotated, TopTable)
colnames(DE) <- c("Ensembl.Gene.ID","logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
DElookup <- list(GeneID="Ensembl.Gene.ID", logFC="logFC", B="B")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 4- Targets information - A matrix containing information about the .Bam les to be used in the analysis ----

samples_file <- read.csv("/media/bo01/Valeria/CRUK/AR_final_code/Dora_AR/samples.csv", 
                         header=TRUE, colClasses = c(NA,"NULL","NULL","NULL","NULL","NULL",NA,NA,NA,"NULL","NULL",NA,NA)) 

fileID <- as.character(samples_file$SRA)
sampleID <- samples_file$Experiment
factorID <- samples_file$Sample.Input
levels(factorID) <- c("Input", "S")

# 4a- Order the files' name in the same order as in the table ----
bam_bai <- dir("/media/bo01/Valeria/CRUK/AR_final_code/Dora_AR/bam", pattern = ".bam")
file_dir <- bam_bai[-c(grep("bai", bam_bai),grep("samstat", bam_bai))]

file_dir_name <- sapply(strsplit(file_dir,"_"), function(i) i[1])

#x nell'ordine di y --> x[order(match(x,y))]
filepath <- file_dir[order(match(file_dir_name,fileID))]

targets <- data.frame(fileID,sampleID,factorID,filepath)

# 4b- Use only the samples you are interested in analysing ----
interest <- c("SRR039769", "SRR039774", "SRR039773", "SRR039775")
targets_int <- targets[which(targets$fileID%in%interest),]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 5- Annotation information ----

# ordered_tfs_position <- ordered_tfs[,c(1,4,5,6,7)]
Split <- strsplit(ordered_tfs$GenomicLocation,":")
chr <- sapply(Split, function(i) i[[1]])
chr <- gsub("chr", "", chr)
start <- sapply(Split, function(i) i[[2]])
end <- sapply(Split, function(i) i[[3]])
str <- sapply(Split, function(i) i[[4]])
ordered_tfs_position <- cbind(ordered_tfs$EnsemblReannotated, chr, start, end, str)
colnames(ordered_tfs_position) <- c("Ensembl.Gene.ID","chr","start", "end", "str")
ChIPannoZones <- defineBins(ordered_tfs_position, zone=c(-25000, 25000), geneID="Ensembl.Gene.ID")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 6- Prior specification ----
DE.prior <- 0.01
prior.mode <- "keepChIP"
prior <- c("D|C"=0.05, "D|notC"=0.005)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 7- Analysis function ----
Dir <- "/media/bo01/Valeria/CRUK/AR_final_code/Dora_AR/bam"
cl <- NULL
Rcade <- RcadeAnalysis(DE, ChIPannoZones, annoZoneGeneidName = "Ensembl.Gene.ID", ChIPtargets = targets_int, ChIPfileDir = Dir, cl=cl, 
                       DE.prior=DE.prior, prior.mode=prior.mode, prior=prior, DElookup=DElookup)
dir.create("Rcade_files_B")
exportRcade(Rcade, directory = "Rcade_files_B",  cutoffMode = "B", cutoffArg = -5, removeDuplicates = "beforeCutoff")

tiff(paste(getwd(),"/Rcade_files_B/plotMM.tiff",sep=""), width=2000, height=1000)
plotMM(Rcade)
dev.off()
tiff(paste(getwd(),"/Rcade_files_B/plotPCA.tiff",sep=""), width=2000, height=1000)
plotPCA(Rcade)
dev.off()
tiff(paste(getwd(),"/Rcade_files_B/plotBB.tiff",sep=""), width=2000, height=1000)
plotBBB(Rcade)
dev.off()
library(rgl) ##required for plotBBB
tiff(paste(getwd(),"/Rcade_files_B/plotBBB.tiff",sep=""), width=2000, height=1000)
plotBBB(Rcade)
dev.off()
#-------------------------------------------------------------------------------------------------------------------

dir.create("Rcade_files_FDR")
exportRcade(Rcade, directory = "Rcade_files_FDR",  cutoffMode = "FDR", removeDuplicates = "beforeCutoff")

tiff(paste(getwd(),"/Rcade_files_FDR/plotMM.tiff",sep=""), width=2000, height=1000)
plotMM(Rcade)
dev.off()
tiff(paste(getwd(),"/Rcade_files_FDR/plotPCA.tiff",sep=""), width=2000, height=1000)
plotPCA(Rcade)
dev.off()
tiff(paste(getwd(),"/Rcade_files_FDR/plotBB.tiff",sep=""), width=2000, height=1000)
plotBBB(Rcade)
dev.off()
library(rgl) ##required for plotBBB
tiff(paste(getwd(),"/Rcade_files_FDR/plotBBB.tiff",sep=""), width=2000, height=1000)
plotBBB(Rcade)
dev.off()


rcade_output <- getRcade(Rcade)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read the output files
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

rcade_results <- read.csv("Rcade_files_B/DEandChIP.csv", header=TRUE, colClasses = c(NA,rep("NULL",11),NA)) 

#Assignt he time point of change to each gene.
probes_timechange <- cbind(ordered_tfs, point_change)#cbind(ordered_tfs[,c(1,2)], point_change)#

probes_timechange_in <- probes_timechange[probes_timechange[,1]%in%rcade_results[,1],]

rcade_results_time <- probes_timechange_in[order(match(probes_timechange_in$ensembl_gene_id, rcade_results$geneID)),]
rcade_results_time <- unique(rcade_results_time)

write.csv(unique(rcade_results_time[,-3]), file="genes_AR_activated.csv", quote=FALSE, row.names=FALSE)

#BEFORE AR TIME POINT OF CHANGE
pc_ar <- unique(probes_timechange[which(probes_timechange[,2]=="AR"),"point_change"])

rcade_before_ar <- rcade_results_time[which(rcade_results_time[,"point_change"]<=pc_ar),]
dim(rcade_before_ar)
probes_timechange_before_ar <- ordered_tfs[which(ordered_tfs[,1]%in%rcade_before_ar[,1]),]


#AFTER AR TIME POINT OF CHANGE
#I include the time_point=NA in this one
rcade_after_ar <- rcade_results_time[-which(rcade_results_time[,"point_change"]<=pc_ar),]
dim(rcade_after_ar)
probes_timechange_after_ar <- ordered_tfs[which(ordered_tfs[,1]%in%rcade_after_ar[,1]),]


#ACTIVATED BY RCADE
ind_ar_activated <- which(ordered_tfs$ensembl_gene_id%in%rcade_results$geneID)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1- Select the TFs that Rcade identified as AR activated ----
genes_AR <- ordered_tfs[ind_ar_activated,]
VCapProbes_ar_activated <- VCapProbes_passed[ind_ar_activated, 2:3]
VCapAndrogen1_ar_activated <- VCapAndrogen1_passed[ind_ar_activated,]
VCapAndrogen2_ar_activated <- VCapAndrogen2_passed[ind_ar_activated,]
VCapControl1_ar_activated <- VCapControl1_passed[ind_ar_activated,]
VCapControl2_ar_activated <- VCapControl2_passed[ind_ar_activated,]
point_change_ar <- point_change[ind_ar_activated]


#Select the TFs that Rcade identified as NON AR activated ----
genes_not_AR <- ordered_tfs[-ind_ar_activated,]
VCapProbes_not_ar_activated <- VCapProbes_passed[-ind_ar_activated, 2:3]
VCapAndrogen1_not_ar_activated <- VCapAndrogen1_passed[-ind_ar_activated,]
VCapAndrogen2_not_ar_activated <- VCapAndrogen2_passed[-ind_ar_activated,]
VCapControl1_not_ar_activated <- VCapControl1_passed[-ind_ar_activated,]
VCapControl2_not_ar_activated <- VCapControl2_passed[-ind_ar_activated,]

VCapAndrogen_mean_not_ar_activated <- VCapAndrogen_passed_mean[-ind_ar_activated,]
point_change_not_ar <- point_change[-ind_ar_activated]
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 2- Cluster the genes not activated by AR using KML ----
percentage_matrix <- function(Matrix){
  perc_matrix <- matrix(0,nrow(Matrix), ncol(Matrix))
  for(i in 1:nrow(Matrix)){
    perc_matrix[i,] <- as.matrix((Matrix[i,]-min(Matrix[i,]))/(max(Matrix[i,])-min(Matrix[i,])))
  }
  
  rownames(perc_matrix) <- rownames(Matrix)
  colnames(perc_matrix) <- colnames(Matrix)
  return(perc_matrix)
}

VCapAndrogen_percentage_not_ar_activated <- percentage_matrix(VCapAndrogen_mean_not_ar_activated)

library(kml)
kml_VCap_not_ar_activated <- cld(VCapAndrogen_percentage_not_ar_activated)
parWithEuclidean <- parALGO(distanceName="euclidean")
kml(kml_VCap_not_ar_activated,nbClusters=5:12,nbRedrawing=10,toPlot='none',parAlgo=parWithEuclidean)

#commands to choose the partition that best suits me.
# It will return popup window which we can use to select the options we like.
X11(type="Xlib")
choice(kml_VCap_not_ar_activated, typeGraph = "bmp")
# use the arrows to see the partitions - see the menu for other changes
# use the space bar to select the partition I want
# use 'm' to export al the result in the current directory. cvs and bmp files


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 3- Plot all genes in all clusters ----

cc <- 9 #seems to be the best fit
dir.create(paste("KML_not_ar_activated_",cc,sep=""))
cluster_file<-read.table(paste(getwd(),"/KML_not_ar_activated_",cc,"/kml_VCap_not_ar_activated-C",cc,"-1-Clusters.csv",sep=""),header=TRUE, sep=";")

cluster_indeces <- lapply(1:cc, function(i) which(cluster_file[,2]==i))

#Plot all genes in all clusters.
tiff(paste(getwd(),"/KML_not_ar_activated_",cc,"/all_",cc,"_clusters.tiff",sep=""), width=3000, height=2000)
Col <- ceiling(length(cluster_indeces)/4)
par(mfrow=c(4,Col))
for(cl in 1:length(cluster_indeces)){
  plot(time_Androgen,VCapAndrogen_percentage_not_ar_activated[cluster_indeces[[cl]][1],],type="l", 
       main=paste("Cluster", cl, sep=" "),
       cex.main=3, xaxt='n', xlab='', ylab='', cex.axis=2, ylim=c(0,1))
  for(i in cluster_indeces[[cl]]){  
    lines(time_Androgen,VCapAndrogen_percentage_not_ar_activated[i,], col=i)
  }
  axis(1, at=time_Androgen, lab=time_Androgen,cex.axis=2)    
}
dev.off()

#Which genes are in each cluster and what is the point-change for each ----

DF_GeneName_PointChange <- lapply(cluster_indeces, function(i) cbind(rownames(VCapAndrogen_percentage_not_ar_activated[i,]), point_change_not_ar[i]))

mean_point_change_cluster <- sapply(DF_GeneName_PointChange, function(i) mean(as.numeric(i[,2]), na.rm=TRUE))
sd_point_change_cluster <- sapply(DF_GeneName_PointChange, function(i) sd(as.numeric(i[,2]), na.rm=TRUE))

tiff(paste(getwd(),"/KML_not_ar_activated_",cc,"/all_",cc,"_clusters_point_change.tiff",sep=""), width=3000, height=2000)
par(mfrow=c(4,Col))
for(i in 1:length(DF_GeneName_PointChange)){
  hist(as.numeric(DF_GeneName_PointChange[[i]][,2]),
       breaks = length(na.omit(as.numeric(DF_GeneName_PointChange[[i]][,2]))),
       cex.axis=2)
}
dev.off()

#Cluster point-change together with the KML clusters (KML+point_change) ----

genes_not_ar_time <- lapply(cluster_indeces, function(i) cbind(genes_not_AR[i,], point_change_not_ar[i]))

time_cluster<-NULL
time_cluster[[1]] <- c(0:5); time_cluster[[2]] <- c(6:12); time_cluster[[3]] <- c(13:24)

cluster_time_not_ar <- lapply(genes_not_ar_time, function(i)  
  lapply(time_cluster, function(j)
    i[which(i[,9]%in%j),]))


VDR_cluster <- which(unlist(lapply(genes_not_ar_time, function(i) "VDR"%in%i[,2])))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PSCAN ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
chrs <- c(seq(1:22), "X", "Y")

get_range_bedtools <- function(Data, Range){
  start <- Data[,3]-Range
  end <- Data[,4]+Range
  strand <- Data[,5]
  strand[strand==1] <- "+"
  strand[strand==-1] <- "-"
  
  Data_range <- data.frame(Data[,2], start, end, strand)
  colnames(Data_range) <- c("chr", "start", "end", "strand")
  
  return(Data_range)
}


# BACKGROUND ------------------------
dir.create("background_range")
dir.create("background_sequences")

background_range <- get_range_bedtools(ordered_tfs[,c(1,4,5,6,7)],1500)

sapply(chrs, function(i) 
  write.table(background_range[background_range[,1]==i,], 
              file=paste(getwd(),"/background_range/background_range",i,".bed", sep=""), 
              col.names=FALSE, quote=FALSE, sep="\t", row.names=FALSE))

#Run bedtools to get the sequences ----
system("bash getfasta_background_bedtools.sh")
#Merge the sequences in one single file
system("cat /media/bo01/Valeria/CRUK/VCap_analysis/background_sequences/background_sequences*.out > /media/bo01/Valeria/CRUK/VCap_analysis/background_sequences/all_sequences.fasta")
system("cat /media/bo01/Valeria/CRUK/VCap_analysis/background_sequences/background_sequences*.out > /home/cri.camres.org/bo01/Downloads/pscan/all_sequences.fasta")

#FOREGROUND ------------------------

# Only KML ----

genes_not_AR_range <- get_range_bedtools(genes_not_AR[,c(1,4,5,6,7)],1500)

cluster_genes_not_AR_range <- lapply(cluster_indeces, function(i) genes_not_AR_range[i,]) #unique

for(j in 1:length(cluster_genes_not_AR_range)){
  dir.create(paste(getwd(),"/foreground_range/cluster", j, sep=""))
  dir.create(paste(getwd(),"/foreground_sequences/cluster", j, sep=""))
  
  for(i in 1:length(chrs)){
    write.table(cluster_genes_not_AR_range[[j]][which(cluster_genes_not_AR_range[[j]][,1]==chrs[i]),],
                file=paste(getwd(),"/foreground_range/cluster", j,"/foreground_range_cluster",j,"_chr",chrs[i],".bed", sep=""),
                col.names=FALSE, quote=FALSE, sep="\t", row.names=FALSE)
  }
}

#Run bedtools to get the sequences ----
system("bash getfasta_foreground_cluster.sh")


#Merge the sequences in each clusters
for(i in 1:cc){
  system(paste("cat ",getwd(),"/foreground_sequences/cluster",i,"/foreground_sequences_cluster",i,"_chr*.fasta > ",getwd(),"/foreground_sequences/cluster",i,"/foreground_sequences_cluster",i,".multifasta", sep=""))
  system(paste("cat ",getwd(),"/foreground_sequences/cluster",i,"/foreground_sequences_cluster",i,"_chr*.fasta > /home/cri.camres.org/bo01/Downloads/pscan/foreground_sequences_cluster",i,".multifasta", sep=""))
}


# KML+Time ----

#write file with the all clusters elements
for(i in 1:length(cluster_time_not_ar)){
  for(j in 1:length(cluster_time_not_ar[[i]])){
    write.table(paste("Cluster ", i, " - Time ", j, sep=""), file="cluster_time_not_ar.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(cluster_time_not_ar[[i]][[j]], file="cluster_time_not_ar.txt", append=TRUE, quote=FALSE, row.names=FALSE, sep="\t")
    write.table("\\", file="cluster_time_not_ar.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}


genes_not_AR_range <- get_range_bedtools(genes_not_AR[,c(1,4,5,6,7)],1500)

cluster_time_genes_not_AR_range <- lapply(cluster_time_not_ar, function(i) 
  lapply(i, function(j) unique(genes_not_AR_range[which(genes_not_AR[,1]%in%j[,1]),])))

lapply(cluster_time_not_ar, function(i) lapply(i, function(j) which(j[,2]=="VDR")))

dir.create(paste(getwd(),"/foreground_range", sep=""))
dir.create(paste(getwd(),"/foreground_sequences", sep=""))

for(j in 1:length(cluster_time_genes_not_AR_range)){
  for(k in 1:length(cluster_time_genes_not_AR_range[[j]])){
    dir.create(paste(getwd(),"/foreground_range/cluster",j,"_time", k, sep=""))
    dir.create(paste(getwd(),"/foreground_sequences/cluster",j,"_time", k, sep=""))
    for(i in 1:length(chrs)){
      write.table(cluster_time_genes_not_AR_range[[j]][[k]][which(cluster_time_genes_not_AR_range[[j]][[k]][,1]==chrs[i]),],
                  file=paste(getwd(),"/foreground_range/cluster",j,"_time", k,
                             "/foreground_range_cluster",j,"_time", k,"_chr",chrs[i],".bed", sep=""),
                  col.names=FALSE, quote=FALSE, sep="\t", row.names=FALSE)    
    }
  }
}

#Run bedtools to get the sequences ----
#Check the cluster number!!!!
system("bash getfasta_foreground_cluster_time.sh")


#Merge the sequences in each clusters
for(i in 1:cc){
  for(j in 1:length(time_cluster)){
    system(paste("cat /media/bo01/Valeria/CRUK/VCap_analysis/foreground_sequences/cluster",i,"_time",j,"/foreground_sequences_cluster",i,"_time",j,"_chr*.fasta > /media/bo01/Valeria/CRUK/VCap_analysis/foreground_sequences/cluster",i,"_time",j,"/foreground_sequences_cluster",i,"_time",j,".multifasta", sep=""))
    system(paste("cat /media/bo01/Valeria/CRUK/VCap_analysis/foreground_sequences/cluster",i,"_time",j,"/foreground_sequences_cluster",i,"_time",j,"_chr*.fasta > /home/cri.camres.org/bo01/Downloads/pscan/foreground_sequences_cluster",i,"_time",j,".multifasta", sep=""))
  }
}

# POSITION WEIGHTS MATRICES ----


#Read the HOMER file with the PFMs (Frequences) --------------------------------------------------------------------------------

homer_pfm <- readLines(paste(getwd(), "homer_custom_motifs.txt", sep="/"))
homer_names1 <- toupper(sapply(strsplit(names(which(sapply(homer_pfm, substr,1,1)==">")),"\t"), function(i) i[2]))
homer_names <- sapply(strsplit(homer_names1, "\\("), function(i) i[1])
ind_name_start_homer <- which(substr(homer_pfm,1,1)==">")

homer_pfm_list <- vector("list", length(ind_name_start_homer))
for(i in 1:length(ind_name_start_homer)){
  if(i==length(ind_name_start_homer)){
    mat <- homer_pfm[(ind_name_start_homer[i]+1):length(homer_pfm)]
  }else{
    mat <- homer_pfm[(ind_name_start_homer[i]+1):(ind_name_start_homer[i+1]-1)]
  }
  mat2 <- matrix(unlist(strsplit(mat, "\t"), use.names=FALSE), ncol = 4, byrow = TRUE)
  
  homer_pfm_list[[i]] <- apply(mat2, 1,as.numeric)
  rownames(homer_pfm_list[[i]]) <- c("A","C","G","T")
}

names(homer_pfm_list) <- homer_names

library(Biostrings)
homer_pwm_all <- lapply(homer_pfm_list, function(i) PWMEnrich::toPWM(i)$pwm)
# Select the homer motifs that apppear in the possible target list
ind_homer_motifs <- which(homer_names%in%ordered_tfs[,2])

homer_pfm_possibletarget <- homer_pfm_list[ind_homer_motifs]
homer_pwm_possibletarget <- homer_pwm_all[ind_homer_motifs]#lapply(homer_pfm_possibletarget, function(i) toPWM(i)$pwm)

" Homer has 20 TFs in the possible_target_ordered. 19 of which are also in jaspar.
I run both because the pfms are a bit different."

length(names(homer_pfm_possibletarget))
34
length(unique(names(homer_pfm_possibletarget)))
32

# for(i in 1:length(homer_pfm_possibletarget)){
#   write.table(paste(">",names(homer_pfm_possibletarget[i]),sep=""), file="homer_pfm.txt", append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
#   write.table(homer_pfm_possibletarget[[i]], file="homer_pfm.txt", append=TRUE, col.names=FALSE, row.names=FALSE) 
# }


#Run MotifDB --------------------------------------------------------------------------------------

query_motif <- function(gene_symbol){
  library(MotifDb)
  # library(MotIV)
  library(seqLogo)
  
  gene.organism.rows <- NULL
  gene_corr <- NULL
  organism.rows <- grep('Hsapiens', values(MotifDb)$organism, ignore.case=TRUE)
  #Some gene symbols are of the form G1::G2. which means they are combined. 
  # I split the gene symbols where there is the :: and put them in a list. 
  # This way the index is the same as in the array.
  gene_symbol_list <- strsplit(values(MotifDb)$geneSymbol,"::")
  
  for(i in 1:length(gene_symbol)){
    #Check if the gene I am interested in is in any position in my list
    gene_symbol_compared <- sapply(gene_symbol_list, function(j) gene_symbol[i]%in%j)
    #Return the indeces
    geneSymbol.rows <- which(gene_symbol_compared==TRUE)
    
    if(length(geneSymbol.rows)>0){gene_corr <- c(gene_corr,gene_symbol[i])}
    
    gene.organism.rows <- c(gene.organism.rows,intersect(geneSymbol.rows,organism.rows))    
    #     matrix_motifs <- c(matrix_motifs,as.list(gene_motif))
  }
  matrix_motif <- MotifDb[gene.organism.rows]
  
  return(matrix_motif)
}


# possible_target_ordered_noempty <- unique(VCapProbes_passed[-which(VCapProbes_passed[,3]==""),3])

motifdb_pfm <- query_motif(unique(as.character(ordered_tfs[,2])))

# motifdb_motif <- export(motifdb_pfm, "motifdb_motif_memeformat", 'meme')

motifdb_pfm_list <- as.list(motifdb_pfm)

motifdb_names <- values(motifdb_pfm)$geneSymbol#sapply(strsplit(names(motifdb_pfm_list), "-"),"[[",3)
names(motifdb_pfm_list) <- motifdb_names

length(motifdb_names)
261
length(unique(motifdb_names))
111

#number of sequences which contributed to the motif. Known only for some.
# values(motifdb_pfm)$sequenceCount

library(PWMEnrich)
motifdb_pwm_possibletarget <- lapply(motifdb_pfm_list, function(i) PWMEnrich::toPWM(i)$pwm)
# class(motifdb_pwm_possibletarget[[1]])  -> matrix
motifdb_pwm_possibletarget_toPWM <- PWMEnrich::toPWM(motifdb_pfm_list)
#class(motifdb_pwm_possibletarget_toPWM[[1]]) -> PWM

#Run JASPAR ---------------------------------------------------------------------------------------

# library(JASPAR2016)
# library(TFBSTools)
# 
# opts <- list()
# # opts[["species"]] <- 9606 #human 10090 #mouse
# opts[["tax_group"]] <- "vertebrates"
# opts[["collection"]] <- "CORE"
# opts[["matrixtype"]] <- "PWM"
# PFMatrixList <- getMatrixSet(JASPAR2016, opts)
# 
# pwm <- Matrix(PFMatrixList)
# names(pwm) <- name(PFMatrixList)
# 
# ind_possible_target <- which(toupper(names(pwm))%in%VCapProbes_passed[,3])
# pwm_possible <- pwm[ind_possible_target]
# in_jaspar_library <- sort(names(pwm_possible))
# 
# 
# length(unique(toupper(c(unique(names(homer_pfm_possibletarget)),unique(motifdb_names),  unique(motifdb_names)))))
# 115
# aa<-unique(toupper(c(unique(names(homer_pfm_possibletarget)),unique(motifdb_names),  unique(motifdb_names))))
# 
# length(aa[-grep("_",aa)])
# 82

# Prepare for PSCAN -----------------------------------------------------------------------------------------------------------------

library(JASPAR2016)
library(TFBSTools)

opts <- list()
# opts[["species"]] <- 9606 #human 10090 #mouse
opts[["tax_group"]] <- "vertebrates"
opts[["collection"]] <- "CORE"
opts[["matrixtype"]] <- "ICM"
PFMatrixList_pfm <- getMatrixSet(JASPAR2016, opts)

pfm_jaspar <- Matrix(PFMatrixList_pfm)
names(pfm_jaspar) <- toupper(name(PFMatrixList_pfm))

# Merge all the pfms from the different databases -----------------------------------------------------------------------------------
list_motifs <- c(pfm_jaspar, homer_pfm_possibletarget, motifdb_pfm_list)

motifs_names <- make.names(names(list_motifs), unique=TRUE)

motifs_names1 <- gsub("\\.", "", as.character(motifs_names))
other_col <- rep(0,length(motifs_names1))
motifs_names2 <- gsub("\\..", "::", motifs_names1)
df_motif <- data.frame(motifs_names1, other_col, motifs_names2)

write.table(df_motif, file=paste(getwd(),"/matrix_list.txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(df_motif, file="/home/cri.camres.org/bo01/Downloads/pscan/matrix_list.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


# PWM - Write each motif in a single file -------------------------------------------------------------------------------------------------
names(list_motifs) <- motifs_names
dir.create("pscan_files")
for(i in 1:length(list_motifs)){
  write.table(list_motifs[[i]], 
              file=paste(getwd(),"/PWM/",df_motif[i,1],".pfm",sep=""), 
              quote=FALSE, col.names=FALSE, row.names=FALSE)#, append=TRUE)
  write.table(list_motifs[[i]], 
              file=paste("/home/cri.camres.org/bo01/Downloads/pscan/",df_motif[i,1],".pfm",sep=""), 
              quote=FALSE, col.names=FALSE, row.names=FALSE)#, append=TRUE)
}


# Run PSCAN ----

setwd("/home/cri.camres.org/bo01/Downloads/pscan")
system("./pscan -p all_sequences.fasta -g")

# only clusters
for(i in 1:cc){
  system(paste("./pscan -q foreground_sequences_cluster",i,".multifasta -m all_sequences.short_matrix ", sep=""))
}

# clusters and tim
for(i in 1:cc){
  for(j in 1:length(time_cluster)){
    system(paste("./pscan -q foreground_sequences_cluster",i,"_time",j,".multifasta -m all_sequences.short_matrix ", sep=""))
    # when it says aborted.. don't worry it's because the datasets is empty for that time point.
  }
}
#copy all the results files into the working directory
system("cp *.res /media/bo01/Valeria/CRUK/VCap_analysis/foreground_sequences")
system("cp *_time*.res /media/bo01/Valeria/CRUK/VCap_analysis/foreground_sequences")

setwd("/media/bo01/Valeria/CRUK/VCap_analysis")


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PSCAN CORRESPONDENCE ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#Read the results for each cluster

foreground_sequences_cluster <- vector("list", cc)
for(i in 1:cc){
  foreground_sequences_cluster[[i]] <- read.table(paste(getwd(),"/foreground_sequences/foreground_sequences_cluster",i,".multifasta.res",sep=""),
                                                  header=TRUE,nrows=814)
}

foreground_sequences_cluster_time <- vector("list", cc)
for(i in 1:cc){
  foreground_sequences_cluster_time[[i]] <- vector("list", length(time_cluster))
  for(j in 1:length(time_cluster)){
    if(file.exists(paste(getwd(),"/foreground_sequences/foreground_sequences_cluster",i,"_time",j,".multifasta.res",sep=""))){
      foreground_sequences_cluster_time[[i]][[j]] <- read.table(paste(getwd(),"/foreground_sequences/foreground_sequences_cluster",i,"_time",j,".multifasta.res",sep=""),
                                                                header=TRUE,nrows=814)
    }
  }
}

#-----------------------------------------------------------------------------------------------------------------------------------
#only cluster
# head(foreground_sequences_cluster_time[[8]][[1]])
# 
# foreground_sequences_cluster_pvalue <- lapply(foreground_sequences_cluster, function(i) i[which(i$P_VALUE<=0.05),])
# 
# foreground_cluster_not_AR <- lapply(foreground_sequences_cluster, function(i) i[which(i$TF_NAME%in%genes_not_AR$hgnc_symbol),])
# foreground_cluster_not_AR_pvalue <- lapply(foreground_cluster_not_AR, function(i) i[which(i$P_VALUE<=0.05),])
# 
# 
# foreground_cluster_AR <- lapply(foreground_sequences_cluster, function(i) i[which(i$TF_NAME%in%genes_AR$hgnc_symbol),])
# foreground_cluster_AR_pvalue <- lapply(foreground_cluster_AR, function(i) i[which(i$P_VALUE<=0.05),])

#-----------------------------------------------------------------------------------------------------------------------------------
#cluster+time
foreground_sequences_cluster_time_pvalue <- lapply(foreground_sequences_cluster_time, function(i) 
  lapply(i, function(j) j[which(j$P_VALUE<=0.05),]))

foreground_cluster_time_not_AR <- lapply(foreground_sequences_cluster_time, function(i) 
  lapply(i, function(j) j[which(j$TF_NAME%in%genes_not_AR$hgnc_symbol),]))
foreground_cluster_time_not_AR_pvalue <- lapply(foreground_cluster_time_not_AR, function(i) 
  lapply(i, function(j) j[which(j$P_VALUE<=0.05),]))


foreground_cluster_time_AR <- lapply(foreground_sequences_cluster_time, function(i) 
  lapply(i, function(j) j[which(j$TF_NAME%in%genes_AR$hgnc_symbol),]))
foreground_cluster_time_AR_pvalue <- lapply(foreground_cluster_time_AR, function(i) 
  lapply(i, function(j) j[which(j$P_VALUE<=0.05),]))
#The genes selected  at the beginning are 412, but the TFs PWM are 814. So the sum of Ar and nor AR is not the total. Some are from
# genes not selected.

#write files for foreground_cluster_time_not_AR_pvalue - NOT AR

for(i in 1:length(foreground_cluster_time_not_AR_pvalue)){
  for(j in 1:length(foreground_cluster_time_not_AR_pvalue[[i]])){
    write.table(paste("Cluster ", i, " - Time ", j, sep=""), file="foreground_cluster_time_not_AR_pvalue.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(foreground_cluster_time_not_AR_pvalue[[i]][[j]], file="foreground_cluster_time_not_AR_pvalue.txt", append=TRUE, quote=FALSE, row.names=FALSE, sep="\t")
    write.table("\\", file="foreground_cluster_time_not_AR_pvalue.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}

#write files for foreground_cluster_time_AR_pvalue - IN AR

for(i in 1:length(foreground_cluster_time_AR_pvalue)){
  for(j in 1:length(foreground_cluster_time_AR_pvalue[[i]])){
    write.table(paste("Cluster ", i, " - Time ", j, sep=""), file="foreground_cluster_time_AR_pvalue.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(foreground_cluster_time_AR_pvalue[[i]][[j]], file="foreground_cluster_time_AR_pvalue.txt", append=TRUE, quote=FALSE, row.names=FALSE, sep="\t")
    write.table("\\", file="foreground_cluster_time_AR_pvalue.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}


# VDR
ensembl_cluster5_1_AR <- ordered_tfs[which(ordered_tfs[,2]%in%foreground_cluster_time_AR_pvalue[[5]][[1]][,1]),]

rcade_results[which(rcade_results[,1]%in%ensembl_cluster5_1_AR[,1]),]



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Genes activated by Ar + another TF ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ar_bs <- read.table()




