#!/usr/bin/env Rscript

# this script is intended to be used as a guide for Paired-End RNAseq raw data preprocessing. Quality control is not included. Labels are exclusive for our experimental design, therefore, they should be renamed for a general usage #
# giving the path to the folder that contains the raw data. The script will generate a new folder called "trimm_data" in which trimmed libraries will be stored #

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (path to the folder that contains raw_data).", call.=FALSE)
}

path = args[1]

base=path
raw_data=paste0(base,"/raw_data")
trimmed_data=paste0(base,"/trimmed_data")

### TRIMMING ###

print('Starts trimming')
setwd(base)
# create a directory for storing the output trimmed files, with a bash call inside the R script #
system("mkdir trimmed_data")
setwd(raw_data)
lista=list.files()
# create a numerical vector with the total number of the raw data files in a sequential way #
index = seq(1,length(lista),2)
# start the trimming loop for all the files in /raw_data # 
for (i in index){
# store forward and reverse labels for the same sample in a variable # 
    R1 = lista[i]
    R2 = lista[i+1]  
# generate the new output labels #
    R1_split=strsplit(as.character(R1), ".fastq")
    R2_split=strsplit(as.character(R2), ".fastq") 
    output1 = paste0(R1_split[[1]][1], '_trimm.fq.gz')
    output2 = paste0(R2_split[[1]][1], '_trimm.fq.gz')     
# store the fastp command as a string #
    command = paste0("fastp -w 16 --n_base_limit 5 cut_front_window_size 1 cut_front_mean_quality 30 --cut_front cut_tail_window_size 1 cut_tail_mean_quality 30 --cut_tail -l 20 -i ", R1, " -I ", R2, " -o ", output1, " -O ", output2)
# execute the trimming #
    system(command)
# move both output files to the trimmed_data directory # 
    system(paste0("mv ", output1, " ../trimmed_data"))
    system(paste0("mv ", output2, " ../trimmed_data"))
}

### GENOME INDEXING ###

print('Starts indexing')
setwd(base)
system("mkdir STAR_genome")
# the folder must contain the reference genome in fasta format #
command = paste0("STAR --runThreadN 8 --runMode genomeGenerate --genomeDir STAR_genome --genomeFastaFiles VCost_genome.fasta")
system(command)

### ALIGNMENT ###

print('Starts alignment')
setwd(trimmed_data)
# prepare the output labels for the BAM files #
for (i in index){
      R1_split=strsplit(as.character(output1), ".fq")
      output_STAR = paste0("sorted_", R1_split[[1]][1], ".bam")
# aligning paired-end libraries and setting the output into BAM format sorted by position #
      command =paste0("STAR --runThreadN 16 --runMode alignReads --genomeDir /home/genomes --readFilesCommand zcat --readFilesIn ", output1, " ", output2, " --outSAMtype BAM SortedByCoordinate")
      system(command)
# removing the trimmed fastq files #
      command = paste0("rm ", output1)
      system(command)
      command = paste0("rm ", output2)
      system(command)
# rename the output BAM files #
      command = paste0("mv Aligned.sortedByCoord.out.bam ", output_STAR)
     system(command)
# removing the output logs #
      system("rm Log*")
      system("rm SJ.out.tab")
      
### SUMMARIZING ###
 
print('Starts summarizing')
for (i in index){
# prepare labels for the output matrix #
    matrix = paste0(R1_split[[1]][1], "_matrix_counts.txt")
    command = paste0('featureCounts -p -T 16 -f -C -a /home/annotations/VCost.v3_25_final.gtf', ' -o ', matrix, ' ', output_STAR)
    system(command)
# we remove the BAM file #
    system(paste0("rm ", output_STAR))
}
      

