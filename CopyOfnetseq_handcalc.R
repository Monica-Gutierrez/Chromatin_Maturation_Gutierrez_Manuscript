
#genes <-  read.csv("~/MANUSCRIPT/temp_figs/corr_complete_UCSCgeneTable_corrected_fullLengthCORR_4NucsMdptsandacf_10minPulse_ALL_PTMS_4_20_18.csv", stringsAsFactors = F)
genes <- read.csv("~/MANUSCRIPT/temp_figs/corr_complete_UCSCgeneTable_corrected_fullLengthCORR_4NucsMdptsandacf_DensFixed_July2018.csv", stringsAsFactors = F)

na <- which(is.na(genes$gene) == T)
genes <- genes[-na, ]

file_name = c("netseq_merged") #"DM538_sacCer3_2017-11-12-02-45") #

for(i in 1:nrow(genes)){
  if(genes$strand[i] == "-"){
  genes$net_start[i] <- genes$end[i] - 400
  genes$net_end[i] <- genes$end[i] + 100
  } else {
    genes$net_start[i] <- genes$start[i] - 100
    genes$net_end[i] <- genes$start[i] + 400
  }
}

file_of_interest.df <- genes 
f=1
#data files
file_name.bam = (paste("/data/home/mpg22/weiner_rando_2015/remaining_marks/",file_name[f],".bam", sep=''))
file_name.bam.bai = paste("/data/home/mpg22/weiner_rando_2015/remaining_marks/",file_name[f],".bam.bai",sep='')

#file_name.bam = (paste("/data/home/mpg22/NET_seq/",file_name[f],".bam", sep=''))
#file_name.bam.bai = paste("/data/home/mpg22/NET_seq/",file_name[f],".bam.bai",sep='')

for (m in 1:nrow(file_of_interest.df)){
  chr = (file_of_interest.df[m,"chr"])
  if(is.na(chr) == T){
    next
  }
  new_start= (file_of_interest.df[m, "net_start"])
  new_end= (file_of_interest.df[m, "net_end"])
  
  chr.gr = GRanges(seqnames= chr, ranges = IRanges(start =new_start , end = new_end ))
  
  p = ScanBamParam(what = c("rname", "strand", "pos", "qwidth"),which = chr.gr)
  
  A_reads.l = scanBam(file = file_name.bam,
                      index = file_name.bam.bai,
                      param = p)
  
  #All the information from the range is in the first entry of the output_reads.l list (this can be       accessed by ou$
  # str(output_reads.l[[1]]) to see list structure
  
  #create a new GenomicRanges object for the reads from this list:
  A_reads.gr = GRanges(seqnames = A_reads.l[[1]]$rname,
                       ranges = IRanges(start = A_reads.l[[1]]$pos, 
                                        width = A_reads.l[[1]]$qwidth),
                       strand = A_reads.l[[1]]$strand)
  
  #changing the genomic ranges so that it uses the midpoint instead of the length of the read.                          
  ranges(A_reads.gr) = IRanges(start=mid(ranges(A_reads.gr)), width=1)
  
  ss_data.df <- as.data.frame(A_reads.gr)
  
  C <-  (nrow(ss_data.df)) #Number of reads mapped to a gene
  #N <- 37576990#Total mapped reads in the experiment
  #L <- genes$net_end[m] - genes$net_start[m] #length in base-pairs for a gene
  
  genes$net_seq_500[m] <-  C/500 #(10^9 * C)/(N * L)
  
  cat("\t done with counts in gene",m ,"\n")
}

write.csv(genes, "~/MANUSCRIPT/temp_figs/corr_complete_UCSCgeneTable_corrected_fullLengthCORR_4NucsMdptsandacf_DensFixed_500bpsNETseq_July2018.csv", row.names = F)


#samtools idxstats /data/home/mpg22/weiner_rando_2015/htz1_sacCer3.bam | awk 'NR<17 {s+=$3} END {print s}'
#h2az_bed[grep("N99733", h2az_bed$id),]
#h2az_bed[grep("N99575", h2az_bed$id),]
#samtools idxstats /data/home/mpg22/weiner_rando_2015/remaining_marks/DM538_sacCer3_2017-11-12-02-45.bam| awk 'NR<17 {s+=$3} END {print s}'
