# I want to investigate the AR motifs in time.
# Does the AR motif change when it activates a gene at 5h compare to one at 10 hours?

genes_corr_illumina_position <- androgen_probes_timechange[androgen_probes_timechange$EnsemblReannotated%in%rcade_results$geneID,]
dim(genes_corr_illumina_position)
genes_corr_position <- ordered_tfs_position[ordered_tfs_position$Ensembl.Gene.ID%in%rcade_results$geneID,]
dim(genes_corr_position)

#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞
# Get the sequences from the UCSC ----

genes_location_range <- genes_corr_position
genes_location_range$chr <- paste("chr",genes_location_range$chr, sep="")
genes_location_range$start <- genes_corr_position$start-25000
genes_location_range$end <- genes_corr_position$end+25000
genes_location_range$str[which(genes_location_range$str==+1)] <- "+"
genes_location_range$str[which(genes_location_range$str==-1)] <- "-"

genes_location_range_opposite <- genes_location_range
genes_location_range_opposite$str[which(genes_location_range_opposite$str=="+")] <- "-"
genes_location_range_opposite$str[which(genes_location_range_opposite$str==-"-")] <- "+"

# Load the database
library("BSgenome.Hsapiens.UCSC.hg38")
ls(1)
sequence_genes <- sequence_genes_opposite <- vector("list", nrow(genes_location_range))
for(i in 1:nrow(genes_location_range)){
  sequence_genes[[i]] <- getSeq(Hsapiens,
                              genes_location_range$chr[i],
                              genes_location_range$start[i], genes_location_range$end[i],
                              width=NA,
                              genes_location_range$str[i])
  sequence_genes_opposite[[i]] <- getSeq(Hsapiens,
                                       genes_location_range_opposite$chr[i],
                                       genes_location_range_opposite$start[i], genes_location_range_opposite$end[i],
                                       width=NA,
                                       genes_location_range_opposite$str[i])
  print(i)
} 

# Write the sequences in a FASTA format file 
for(i in 1:length(sequence_genes)){
  write.table(paste(">",genes_corr_illumina_position$EnsemblReannotated[i],".",genes_corr_illumina_position$IlluminaID[i],
                    ".", genes_corr_illumina_position$GenomicLocation[i],
                    sep=""), file="ARactivated_sequences.fasta", append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
  write.table(toString(sequence_genes[[i]]), file="ARactivated_sequences.fasta",
              append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  print(paste(i, nchar(toString(sequence_genes[[i]])),sep="-"))
  
}

#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞

#Get the motifs in MEME format ----
"Got to meme suite  -> download and install -> download databases -> 
jaspar core 201 (non redundant) -> MA009.3 Ar"
"Copy and paste it in a new txt file and save it as ARmotif.meme"

#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞

# Run FIMO from MEME suite 4.11.2 in home ----
# Use the following command in home (bo01)
# fimo /Volumes/Valeria/CRUK/AR_downstream_analysis/VCAP/ARmotif.meme /Volumes/Valeria/CRUK/AR_downstream_analysis/VCAP/ARactivated_sequences.fasta

# Read the fimo results ----
# delete the # fromt he fimo.txt file
fimo_out <- read.table("fimo_out/fimo.txt", header=TRUE, sep="\t")
head(fimo_out)

fimo_split <- strsplit(as.character(fimo_out$sequence.name), "[.]")

ordered_fimo_genes_corr <- ordered_fimo_genes_corr_location <- NULL
for(i in 1:nrow(fimo_out)){
  ind_gene_corr <- which((genes_corr_illumina_position$IlluminaID==fimo_split[[i]][2])&
                           (genes_corr_illumina_position$EnsemblReannotated==fimo_split[[i]][1])&
                           (genes_corr_illumina_position$GenomicLocation==fimo_split[[i]][3]))
  ordered_fimo_genes_corr <- rbind(ordered_fimo_genes_corr, genes_corr_illumina_position[ind_gene_corr,])
  ordered_fimo_genes_corr_location <- rbind(ordered_fimo_genes_corr_location, genes_corr_position[ind_gene_corr,])
 
}

corresponding_all <- cbind(fimo_out,ordered_fimo_genes_corr)
corresponding_all_location <- cbind(fimo_out,
                                    ordered_fimo_genes_corr_location,
                                    corresponding_all$ordered_point_change,
                                    corresponding_all$androgen_ordered_point_change)

hist(table(corresponding_all$matched.sequence))

motifs <- corresponding_all$matched.sequence
levels(motifs) <- 1:length(levels(motifs))

plot(as.numeric(as.character(motifs)), corresponding_all$ordered_point_change)
#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞

#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞
# Difference in TIME ----
#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞

#Set the time windows I want to explore 
win <- 12
time_range <- NULL
for(i in 0:(max(corresponding_all$ordered_point_change)/3)){
  time_range <- c(time_range, pc_ar+(win*i))
}

#GEt the indeces of the genes that change within that time range (t1<x<=t2) == (x<=t2 - x<=t1). I use setdiff!
data_time_ind <- vector("list", length(time_range))
for(i in 1:length(time_range)){
  data_time_ind[[i]] <- setdiff(which(corresponding_all$ordered_point_change<=time_range[i]), unlist(data_time_ind))
  print(i)
}

corresponding_all_time <- lapply(data_time_ind, function(i) corresponding_all[i,])

# Check what are the motifs that appear most in each time window

# Create a Percentage matrix position for each time window ----
sequence_window <- lapply(corresponding_all_time, function(i) DNAStringSet(as.character(i$matched.sequence)))

# Create a fake PWM ----
pcm_sequences_time <- lapply(sequence_window, function(i) if(length(i)!=0){consensusMatrix(i,as.prob=TRUE)[1:4,]})

# Use the e' bioconductor package seqLogo
library(seqLogo)
pdf(paste("Figures/motif_logo_at_time",win,".pdf",sep=""))
par(mfrow=c(2,2))
lapply(pcm_sequences_time, function(i) seqLogo(makePWM(i)))
dev.off()

#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞
#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞
# Difference in TIME ANDROGEN ----
#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞

#Set the time windows I want to explore 
for(win in 3:12){
  time_range <- NULL
  for(i in 0:(max(corresponding_all$androgen_ordered_point_change, na.rm=TRUE)/3)){
    time_range <- c(time_range, pc_ar+(win*i))
  }
  
  #GEt the indeces of the genes that change within that time range (t1<x<=t2) == (x<=t2 - x<=t1). I use setdiff!
  data_time_ind <- vector("list", length(time_range))
  for(i in 1:length(time_range)){
    data_time_ind[[i]] <- setdiff(which(corresponding_all$androgen_ordered_point_change<=time_range[i]), unlist(data_time_ind))
    print(i)
  }
  
  corresponding_all_time <- lapply(data_time_ind, function(i) corresponding_all[i,])
  
  # Check what are the motifs that appear most in each time window
  
  # Create a Percentage matrix position for each time window ----
  sequence_window <- lapply(corresponding_all_time, function(i) DNAStringSet(as.character(i$matched.sequence)))
  
  # Create a fake PWM ----
  pcm_sequences_time <- lapply(sequence_window, function(i) if(length(i)!=0){consensusMatrix(i,as.prob=TRUE)[1:4,]})
  
  # Use the e' bioconductor package seqLogo
  library(seqLogo)
  pdf(paste("Figures/androgen_motif_logo_at_time",win,".pdf",sep=""))
  for(ll in 1:length(pcm_sequences_time)){
    if(length(pcm_sequences_time[[ll]])>0){
      seqLogo(pcm_sequences_time[[ll]])
    }else{
      plot(1)
    }
  }
  dev.off()
}

#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞

#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞
# Difference in SPACE
#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞

corresponding_all_location$start-25000



#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞

#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞
# Difference in TIME and SPACE
#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞



#∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞

















