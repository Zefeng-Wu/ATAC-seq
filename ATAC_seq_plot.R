# make plot of ATAC-seq
require(stringr)
require(tidyverse)
require(ggpubr)

##################################
###### 1.1 raw reads summary
##################################

Reads_Count <- function(fastq_file){
  system("date")
  message("Calculate the reads number of sample: ",fastq_file)
  reads_num_in_file <- as.numeric(str_split_fixed(system(paste(c("wc -l ",fastq_file),collapse =  ""),intern = TRUE),pattern = " ",n = 2)[,1])
  return(reads_num_in_file)
}

file_list<-dir("/data/wzf/Medicago_ATAC_public/1Fastq",full.names = TRUE)
sample_info <-data.frame(Sample_file_name=file_list,stringsAsFactors = FALSE)
sample_info$reads_num <-apply(sample_info,1,function(x)Reads_Count(x[1])/4) # Caution!!!
sample_info$Sample_id <- str_split_fixed(basename(sample_info$Sample_file_name),pattern = "_[0-9].fastq",n = 2)[,1]
#plot
reads_summary<-sample_info%>%group_by(Sample_id)%>%summarise(Reads_number=sum(reads_num))
ggtexttable(reads_summary, rows = NULL,theme = ttheme("mOrange")) #theme = ttheme("blank")

##################################################################
######### 1.2 raw mapped reads summary (based on the filenames in the 1.1)
##################################################################

Alignment_Count <-function(bam_file){
  system("date")
  message("Calculate the reads number of sample: ",bam_file)
  reads_num_in_file <- as.numeric(str_split_fixed(system(paste(c("samtools view -c -F 260 -@ 20 ",bam_file),collapse =  ""),intern = TRUE),pattern = " ",n = 2)[,1])
  return(reads_num_in_file)
}

file_list<-dir("/data/wzf/Medicago_ATAC_public/2raw_alignement_results",full.names = TRUE)
reads_summary$Aligned_reads <-apply(reads_summary,1,function(x) Alignment_Count(file_list[grep(x[1],file_list)]))
ggtexttable(reads_summary, rows = NULL,theme = ttheme("mOrange"))

####################################################################
######### 1.3 uniq mapped reads summary (based on the filenames in the 1.1)
###################################################################

file_list <- grep(dir("/data/wzf/Medicago_ATAC_public/5remove_duplicated_results/",full.names = TRUE),invert = T,value = T,pattern = ".bai")
reads_summary$Uniq_mappped_reads <-apply(reads_summary,1,function(x) Alignment_Count(file_list[grep(x[1],file_list)]))
ggtexttable(reads_summary, rows = NULL,theme = ttheme("mOrange"))


####################################################################
######### 1.4 peak number summary (based on the filenames in the 1.1)
###################################################################
library(ChIPpeakAnno)
library(rtracklayer)
Peak_number<-function(peak_file){
  message(peak_file)
  peak_gr<- import.bed(peak_file)
  return(length(peak_gr))
}
reads_summary<-reads_summary[reads_summary$Sample_id!="MT-gDNA",]
file_list <-grep(dir("../Medicago_ATAC_public/8b_peak_calling_result_homer",full.names = TRUE),value = T,pattern = ".bed")
reads_summary$Peak_number <-apply(reads_summary,1,function(x) Peak_number(file_list[grep(x[1],file_list)]))
ggtexttable(reads_summary, rows = NULL,theme = ttheme("mOrange"))


####################################################################
######### 1.5 peak distribution using merged peak at least two replicates
###################################################################
library(ChIPseeker)
library(GenomicFeatures)
library(ggsci)

tx <- makeTxDbFromGFF("~/genome/Medicago_truncatula/Medicago_truncatula.MedtrA17_4.0.31.gtf")
options(ChIPseeker.downstreamDistance = 1000)

plot_list <-list()
file_list <-dir("../Medicago_ATAC_public/9merged_replicated_peaks",full.names = TRUE,pattern = ".bed")
for (peak in file_list){
  if (!endsWith(x = peak,suffix = "resize.bed")){
    message(peak)
    peak_name <- basename(peak)
    peak_gr<-import.bed(peak)
    #covplot(peak_gr)
     #plotAvgProf2(peak_gr,TxDb = tx,upstream = 2000,downstream = 2000)
    peakAnno <- annotatePeak(peak_gr, tssRegion=c(-2000, 0),TxDb=tx)
    #plotAnnoPie(peakAnno)
    data <- peakAnno@annoStat
    data$ymax <- cumsum(data$Frequency) ## Compute the cumulative percentages (top of each rectangle)
    data$ymin <- c(0, head(data$ymax, n=-1)) ## Compute the bottom of each rectangle
    data[,c(2,3,4)]<-round(data[,c(2,3,4)],digits = 2)
    p<-ggplot(data = data, aes(ymax=ymax, ymin=ymin,xmax=4,xmin=3, fill = Feature))+
      geom_rect()+
      coord_polar(theta = "y") +
      geom_label(aes(label=paste(Frequency,"%"),x=3.5,y=(ymax+ymin)/2),inherit.aes = TRUE, show.legend = FALSE) +
      theme_void() +
      scale_fill_npg()+
      xlim(1,4)+
      theme(legend.position=c(.5, .5))+
      ggtitle(peak_name)
    plot_list[[peak_name]]<-p
  }
}

library(ggpubr)
ggarrange(plotlist = plot_list,common.legend = TRUE,ncol = 3)+theme(plot.margin = unit(c(2,2,2,2), "lines"),plot.title = element_text(hjust = 0.5))



####################################################################
######### 1.6 diffbind
###################################################################
#---> step1 prepare sample info datafile

sample_info <-data.frame(SmapleID=c("Control","Control","Control","Submerge","Submerge","Submerge","RootTip","RootTip"),
                         Tissue = "Rootip",
                         Factor="Root",
                         Condition=c("Control","Control","Control","Submerge","Submerge","Submerge","RootTip","RootTip"),
                         Treatment="Soil",
                         Replicate=c(1,2,3,1,2,3,1,2),
                         bamReads=c("7downsampled_replicates_results/Mt_ATAC_C1.uniq.sorted.dedup.bam",
                                    "7downsampled_replicates_results/Mt_ATAC_C2.uniq.sorted.dedup.bam",
                                    "7downsampled_replicates_results/Mt_ATAC_C3.uniq.sorted.dedup.bam",
                                    "7downsampled_replicates_results/Mt_ATAC_S1.uniq.sorted.dedup.bam",
                                    "7downsampled_replicates_results/Mt_ATAC_S2.uniq.sorted.dedup.bam",
                                    "7downsampled_replicates_results/Mt_ATAC_S3.uniq.sorted.dedup.bam",
                                    "7downsampled_replicates_results/Mt_root.tip.ATAC_rep1.uniq.sorted.dedup.bam",
                                    "7downsampled_replicates_results/Mt_root.tip.ATAC_rep2.uniq.sorted.dedup.bam"),
                         controlID=c("Contorl","Control","Control","Submerge","Submerge","Submerge","RooTip","RootTip"),
                         bamControl=NA,
                         Peaks =  c("8b_peak_calling_result_homer/Mt_ATAC_C1.uniq.sorted.dedup.peak.bed",
                                   "8b_peak_calling_result_homer/Mt_ATAC_C2.uniq.sorted.dedup.peak.bed",
                                   "8b_peak_calling_result_homer/Mt_ATAC_C3.uniq.sorted.dedup.peak.bed",
                                   "8b_peak_calling_result_homer/Mt_ATAC_S1.uniq.sorted.dedup.peak.bed",
                                   "8b_peak_calling_result_homer/Mt_ATAC_S2.uniq.sorted.dedup.peak.bed",
                                   "8b_peak_calling_result_homer/Mt_ATAC_S3.uniq.sorted.dedup.peak.bed",
                                   "8b_peak_calling_result_homer/Mt_root.tip.ATAC_rep1.uniq.sorted.dedup.peak.bed",
                                   "8b_peak_calling_result_homer/Mt_root.tip.ATAC_rep2.uniq.sorted.dedup.peak.bed"),
                         PeakCaller="bed")
#### save and reload
write.csv(sample_info,file = "result_picture/sample_info.csv",row.names = FALSE)
tamoxifen <- dba(sampleSheet="result_picture/sample_info.csv")

## plot sample similarity
plot(tamoxifen)
tamoxifen <- dba.count(tamoxifen, summits=250) #https://www.biostars.org/p/401188/

### differental
tamoxifen <- dba.contrast(tamoxifen, categories=DBA_CONDITION)
plot(tamoxifen)
tamoxifen <- dba.analyze(tamoxifen,bFullLibrarySize = TRUE)
plot(tamoxifen, contrast=1)
tamoxifen.DB <- dba.report(tamoxifen)

## pca plot
dba.plotPCA(tamoxifen,DBA_TISSUE,label=DBA_CONDITION)

dba.plotPCA(tamoxifen, contrast=1,label=DBA_TISSUE)
dba.plotVolcano(tamoxifen)


####################################################################
######### 1.7 fragment size
###################################################################
library(GenomicAlignments)
library(magrittr)
library(dplyr)
library(ggplot2)

atacReads <- readGAlignmentPairs("7downsampled_replicates_results/Mt_ATAC_C1.uniq.sorted.dedup.bam", param = ScanBamParam(mapqFilter = 2, 
                                            flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE),
                                            what = c("qname", "mapq", "isize")))
atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)

fragLenPlot <- table(insertSizes)%>%data.frame(stringsAsFactors = FALSE) 
colnames(fragLenPlot)<-c("InsertSize","Count")
fragLenPlot$InsertSize<-as.numeric(fragLenPlot$InsertSize)
ggplot(fragLenPlot,aes(x = InsertSize, y = Count)) + 
  geom_line(color="steelblue")+
  theme_bw()+
  theme(text=element_text(size = 20))+
  geom_vline(xintercept = c(180, 247), colour = "red") + geom_vline(xintercept = c(315,  437), colour = "darkblue") + 
  geom_vline(xintercept = c(100), colour = "darkgreen")

####################################################################
######### 1.8 soggi
###################################################################
require(soGGi)

genes <-genes(tx)
Get_marks_from_filenames<-function(peak_file_name){
  basename<-basename(peak_file_name)
  mark_name<-unlist(strsplit(basename,"\\."))[1]
  return(mark_name)
}

df<-data.frame(xIndex=NA,Sample=NA,Score=NA)

plot_data<-function(histone_file){
  chipExample1 <- regionPlot(histone_file,samplename = Get_marks_from_filenames(histone_file),testRanges = genes,format = "bam",distanceUp = 2000,distanceDown = 2000,style = "percentOfRegion",method = "spline")
  df<-rbind(df,plotRegion(chipExample1)$data)
  return(df)
}

bam_signal_list <-list()
## run functuon for each bam
bam_file_list <- grep(dir("../Medicago_ATAC_public/7c_merge_replicates_bams",full.names = TRUE),invert = T,value = T,pattern = ".bai")
for (bam in bam_file_list){
  bam_signal_list[basename(bam)]<-plot_data(bam)
}

df_final<-rbind(dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8,dd9,dd10,dd11) # change top 5
df_final$Level<-stringr::str_split_fixed(df_final$Sample,pattern = "_",2)[,2]
df_final$Sample<-stringr::str_split_fixed(df_final$Sample,pattern = "_",2)[,1]
df_final$Sample<-factor(df_final$Sample,levels=c("H3K4me3","H3K36me3","H3K9ac","H3K27ac","H3K4me2","DNaseI","H3K9me2","H3K4me1","H3K27me1","H3K27me3","Control"))
p<-ggplot(df_final,aes(x=xIndex,y=Score,colour=Level))+geom_line()+scale_x_discrete(limits=seq(1,400,100),labels=c("-1kb","TSS","TTS","1kb"))+theme(text = element_text(size=18),plot.margin = unit(c(1,1,1,1), "cm"))+ylab("Mean RPM")+xlab("Region")+facet_wrap(~Sample,nrow=3,scales = "free")


##################################
#### Gviz; genome track viwer/ IGV/IGB/trackViwer 
#################################
library(Gviz)
library(GenomicFeatures)

options(ucscChromosomeNames=FALSE)

focus_gene <- "MTR_8g494210"

mt.tx <-makeTxDbFromGFF("~/genome/Medicago_truncatula/Medicago_truncatula.MedtrA17_4.0.31.gtf")
genes <-genes(mt.tx)

region.start <- start(genes[genes$gene_id==focus_gene])-10000
region.end <- end(genes[genes$gene_id==focus_gene])+10000


## make different track
axis_track = GenomeAxisTrack()
gene_track = GeneRegionTrack(tr,fill="blue",lwd=2,lex=1,col.line="red",col="transparent",lty=1,transcriptAnnotation=TRUE)

plot_track_list = list(axis_track,gene_track)
                  
plotTracks(plot_track_list,    
           from = region.start,
           to = region.end,
           showID = T,   
           sizes = c(1,1),
           chromosome = region.chr_name,
           genome = "mtr17")

##################
### Option2
library(rtracklayer)
options(ucscChromosomeNames=FALSE)
focus_gene <- "MTR_8g494210"

mt.gtf <-import("~/genome/Medicago_truncatula/Medicago_truncatula.MedtrA17_4.0.31.gtf")
#mt.gtf<-mt.gtf[mt.gtf$type!="transctipt"]
mt.gtf.new <- GRanges(mt.gtf,
                      source = mt.gtf$source,
                      type = mt.gtf$type,
                      score = mt.gtf$score,
                      phase = mt.gtf$phase,
                      gene = mt.gtf$gene_id,
                      transcript = mt.gtf$transcript_id,
                      symbol = sprintf("%s:%s",mt.gtf$gene_id,mt.gtf$transcript_id))

RefSeq.gtf.select = mt.gtf.new[mt.gtf.new$gene == focus_gene]     # 先把所有名字叫TP53的基因选出来

region.start = min(start(RefSeq.gtf.select)) - 10000  #基因序列开始的位置
region.end = max(end(RefSeq.gtf.select)) + 10000  #基因序列结束的位置
region.chr_name = as.character(seqnames(RefSeq.gtf.select))[1]

#下面我们要把这个区域里面所有和这个区域有关的基因全部都选择出来
select_region.GR_obj = GRanges(seqnames = region.chr_name,
                               ranges = IRanges(start = region.start,end = region.end))
select_region.GR_obj

region.start <- start(genes[genes$gene_id==focus_gene])-10000
region.end <- end(genes[genes$gene_id==focus_gene])+10000

# select region all genes
region.Refseq.select = mt.gtf.new[mt.gtf.new %over% select_region.GR_obj]
region.Refseq.select<-region.Refseq.select[region.Refseq.select$type!="gene"&
                                          region.Refseq.select$type!="transcript"]
## make different track
##1.make axis and gene model track
axis_track = GenomeAxisTrack(range = IRanges(start = region.start,end = region.end))
gene_track = GeneRegionTrack(region.Refseq.select,transcriptAnnotation = TRUE,
                             chromosome = region.chr_name,
                             name = "Gene model",
                             fill="blue4",
                             lwd=2,
                             lex=1,
                             col.line="red3",
                             col="transparent",
                             lty=1)
                             #background.panel = "#FFFEDB")


## make peak track
library(ChIPseeker)
peak1 <- readPeakFile("../Medicago_ATAC_public/9merged_replicated_peaks/Mt_ATAC_C.bed")
peak2 <- readPeakFile("../Medicago_ATAC_public/9merged_replicated_peaks/Mt_ATAC_S.bed")

peak1_track <- AnnotationTrack(peak1,name = "MT.control")
peak2_track <- AnnotationTrack(peak1,name = "MT.submerge")

### make alignment peak (it is prefer to use normalzied bigwig or bam after merge all replicate)
align_track1 <- AlignmentsTrack("../Medicago_ATAC_public/7downsampled_replicates_results/Mt_ATAC_C1.uniq.sorted.dedup.bam",
                               type = "coverage",
                               col.coverage="Orange",
                               fill.coverage="Orange",
                               name="MT.control")
align_track2 <- AlignmentsTrack("../Medicago_ATAC_public/7downsampled_replicates_results/Mt_ATAC_S1.uniq.sorted.dedup.bam",
                                type = "coverage",
                                col.coverage="green4",
                                fill.coverage="green4",
                                name="MT.control")


plot_track_list = list(axis_track,peak1_track,peak2_track,align_track1,align_track2,gene_track)

plotTracks(plot_track_list,    
           from = region.start,
           to = region.end,
           showID = T,   
           sizes = c(0.5,0.5,0.5,1,1,1),
           chromosome = region.chr_name,
           col = NULL,
           rotation.title=0,
           cex.title=1,               # title size
           title.width=2,             # title panel width
           #extend.left=10,
           
           background.title = "gray", # left panel colot
           col.border.title="white",  # left panel border color
           col.title="black")         # title name color



### ggtern 可视化三种处理下基因的表达和atac-seq


