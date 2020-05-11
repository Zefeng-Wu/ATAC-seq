#!/bin/bash
export LANG="en_US.UTF-8"

### FastQC inspect
fastq_files_dir="/data/wzf/Medicago_ATAC_public/1Fastq"
cd $fastq_files_dir/../ # as the working directory
mkdir 1a_fastqc
for m in $(ls $fastq_files_dir/*.fastq); do echo $m; fastqc $m  -o 1a_fastqc; done
multiqc 1a_fastqc

### Fastq reads adapter detecting

mkdir 1b_adapters
adapter_detector="/home/WZF/softwares/detect_adapter.py"

for m in $(ls $fastq_files_dir/*.fastq); do echo $m; python3 $adapter_detector $m > 1b_adapters/$(basename $m).adpt;done

### FASTQ read adaptor trimming

Adapter_error_rate=​0.1
mkdir 1c_trimed_fastq

for fq in $(ls $fastq_files_dir/*.fastq | sed 's/_[1-2].fastq//g'| uniq)
do echo $fq
	prefix=$(basename $fq)
	a_adapter=$(cat 1b_adapters/${prefix}_1.fastq.adpt)
	A_adapter=$(cat 1b_adapters/${prefix}_2.fastq.adpt)

	if ([ -z $(grep '[^[:space:]]' 1b_adapters/${prefix}_1.fastq.adpt) ] && [ -z $(grep '[^[:space:]]' 1b_adapters/${prefix}_2.fastq.adpt) ]) ;then    
			echo "No adapter needed to remove, just copy the orignal raw reads into the tirmed fastq directory !"
			cp ${fq}_1.fastq 1c_trimed_fastq/${prefix}_1.fastq
			cp ${fq}_2.fastq 1c_trimed_fastq/${prefix}_2.fastq
	fi
		
	if ([ ! -z $(grep '[^[:space:]]' 1b_adapters/${prefix}_1.fastq.adpt) ] && [ ! -z $(grep '[^[:space:]]' 1b_adapters/${prefix}_2.fastq.adpt) ]) ;then
		echo "Both ends have adapters needed to remove !"
		cutadapt -m 5 -e $Adapter_error_rate -a $a_adapter -A $A_adapter -o 1c_trimed_fastq/${prefix}_1.fastq -p 1c_trimed_fastq/${prefix}_2.fastq  ${fq}_1.fastq ${fq}_2.fastq >>cutadapt.logs
	fi

	if ([ ! -z $(grep '[^[:space:]]' 1b_adapters/${prefix}_1.fastq.adpt) ] && [ -z $(grep '[^[:space:]]' 1b_adapters/${prefix}_2.fastq.adpt) ]) ;then
		echo "Left reads have adapter needed to remove !"
		cp 1b_adapters/${prefix}_1.fastq.adpt 1b_adapters/${prefix}_2.fastq.adpt
		cutadapt -m 5 -e $Adapter_error_rate -a $a_adapter -A $A_adapter -o 1c_trimed_fastq/${prefix}_1.fastq -p 1c_trimed_fastq/${prefix}_2.fastq  ${fq}_1.fastq ${fq}_2.fastq >>cutadapt.logs
	fi

	if ([ -z $(grep '[^[:space:]]' 1b_adapters/${prefix}_1.fastq.adpt) ] && [ ! -z $(grep '[^[:space:]]' 1b_adapters/${prefix}_2.fastq.adpt) ]) ;then
		echo "Right reads have adapter nedded to remove!"
		cp 1b_adapters/${prefix}_2.fastq.adpt 1b_adapters/${prefix}_1.fastq.adpt
		cutadapt -m 5 -e $Adapter_error_rate -a $a_adapter -A $A_adapter -o 1c_trimed_fastq/${prefix}_1.fastq -p 1c_trimed_fastq/${prefix}_2.fastq  ${fq}_1.fastq ${fq}_2.fastq >>cutadapt.logs
	fi
done

### Read alignment (Bowtie2 aligner)

mkdir 2raw_alignement_results
mkdir 3uniq_mapped_results
mkdir 4uniq_mapped_results_sort
mkdir 5remove_duplicated_results
mkdir 6remove_Mt_Pt_results

bwt2_idx=/home/WZF/genome/Medicago_truncatula/Medicago_truncatula
mapping_threads=15

for fq in $(ls 1c_trimed_fastq/*.fastq | sed 's/_[1-2].fastq//g'| uniq); do echo $fq; prefix=$(basename $fq); bowtie2 -p $mapping_threads -x $bwt2_idx -1 ${fq}_1.fastq -2 ${fq}_2.fastq  --no-unal --no-mixed --no-discordant -S 2raw_alignement_results/${prefix}.sam >>bt2_mapping_log; done  

### get uniq mapping alignment (or mapping quality>=2) # will lead to un-paired reads
for m in $(ls 2raw_alignement_results/*.sam); do samtools view -Sh $m | grep -e "^@" -e "XM:i:[012][^0-9]" | grep -v "XS:i:"  >3uniq_mapped_results/$(basename ${m%.sam}).uniq.sam;done # some strict
# nohup sh -c 'for m in $(ls 2raw_alignement_results/*.sam); do samtools view -Sh $m | grep -e "^@" -e "XM:i:[012][^0-9]" | grep -v "XS:i:"  >3uniq_mapped_results/$(basename ${m%.sam}).uniq.sam;done' &
### convert sam to sorted bam
for m in $(ls 3uniq_mapped_results/*.sam); do echo $m; picard SortSam  INPUT=$m OUTPUT=4uniq_mapped_results_sort/$(basename ${m%.sam}).sorted.bam SO=coordinate; done

### Remove PCR duplicated reads (or Picard-tools)
for m in $(ls 4uniq_mapped_results_sort/*.bam); do echo $m; samtools rmdup $m 5remove_duplicated_results/$(basename ${m%.bam}).dedup.bam; done
#for m in $(ls 4uniq_mapped_results_sort/*.bam);do picard MarkDuplicates INPUT=$m OUTPUT=5remove_duplicated_results/$m METRICS_FILE=5remove_duplicated_results/$m.mark_dup_metrix.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false


### Remove the reads of chloroplast or mitochondrial genomes
for m in $(ls 5remove_duplicated_results/*.bam); do samtools index $m;done
for m in $(ls 5remove_duplicated_results/*.bam); do samtools idxstats $m | cut -f 1 | grep -v 'Mt|Pt' | xargs samtools view -b $m > 6remove_Mt_Pt_results/$(basename $m);done

### Random sampling to make replicates with eaqul number of reads
mkdir 7downsampled_replicates_results
samples=$(ls 6remove_Mt_Pt_results/*.bam | xargs -n 1 basename | sed 's/\.uniq\.sorted\.dedup\.bam//g' | sed -r 's/[0-9]{1,}$//g'| sed s'/*\///g' | sort| uniq)

for sample in $samples
do 	#echo $sample #Mt_ATAC_C
	
	sample_repeats_files=$(ls 6remove_Mt_Pt_results/$sample*)
	sample_repeats_num=$(ls 6remove_Mt_Pt_results/$sample*|wc -l) #3
	
	if [ $sample_repeats_num == 1 ]
	then 
		echo $sample ": has no replicates, just copy the bam file!"	
		cp 6remove_Mt_Pt_results/$sample* 7downsampled_replicates_results/
	
	else
		temp_min_reads_number=()

		for sample_repeat_file in $sample_repeats_files
		do  
			count=$(samtools view -c $sample_repeat_file)
			temp_min_reads_number+=($count)
		#echo ${temp_min_reads_number[@]}
		done
		
		{ read min; max=$(tail -n1); } < <(printf "%s\n" "${temp_min_reads_number[@]}" | sort -n)
		min_reads_number=$min
				
		for sample_repeat_file in $sample_repeats_files
		do 
			count=$(samtools view -c $sample_repeat_file)
			if [ $count == $min_reads_number ]
			then
				echo $sample": has replicates, and the minumum reads number is: " ,$min_reads_number,"and this the replcate with minimal reads:" $(basename $sample_repeat_file)
				cp $sample_repeat_file 7downsampled_replicates_results
			else
				pick_up_ratio=`echo "scale=4; $min_reads_number/$count"|bc`
				echo $sample": has replicates, and the minumum reads number is: " ,$min_reads_number, "; the present pick up ratio is: ", $pick_up_ratio
				picard DownsampleSam I=$sample_repeat_file O=7downsampled_replicates_results/$(basename $sample_repeat_file) P=$pick_up_ratio
			fi
		done
		
	fi
done


### fragment shift bam (deeptools:may lose some information,ATACseqQC,or bedpe)
mkdir 7a_fragment_shift_result
#for m in $(ls 7downsampled_replicates_results/*.bam); do alignmentSieve --numberOfProcessors 25 --ATACshift --bam $m -o 7a_fragment_shift_result/$(basename $m).tmp.bam;done
#for m in $(ls 7a_fragment_shift_result/*.tmp.bam); do samtools sort -@ 20 -O bam -o 7a_fragment_shift_result/$(basename ${m%.tmp.bam}) $m; samtools index -@ 20 7a_fragment_shift_result/$(basename ${m%.tmp.bam}); rm $m; done
R script ATACseqQC.R 

### insepct fragment size distribution
mkdir 7b_fragemnt_size_distribution
R script ATACseqQC.R 

## peak calling (homer,macs,hopspot:FDR < 0.01 with 8 reads,HMMRATAC)

mkdir 8a_peak_calling_result_hmmratac #(erro with no tag )

for m in $(ls 7a_fragment_shift_result/*.bam); do echo $m; samtools view -H $m| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > 8a_peak_calling_result_hmmratac/genome.info; java -jar ~/softwares/HMMRATAC_V1.2.10_exe.jar -b $m -i $m.bai -g 8a_peak_calling_result_hmmratac/genome.info -o 8a_peak_calling_result_hmmratac/$(basename $m) ;done
#---------------------------------------------------------------------------------------------------------------
mkdir 8b_peak_calling_result_homer
for m in $(ls 7a_fragment_shift_result/*.bam);do makeTagDirectory 8b_peak_calling_result_homer/$(basename $m) $bam -format bam -mapq 2; done 
for m in $(ls 7a_fragment_shift_result/*.bam);do Findpeaks 8b_peak_calling_result_homer/$(basename $m)/ -region -minDist 150 -o 8b_peak_calling_result_homer/$(basename $m);done  #(pc,2018,2018 pj)

### peak merger between replciates by bedtools to get union THS sites (or pooled data, and two replicated peak)
mkdir 9merged_replicated_peaks

bedtools intersect -a a.bed -b b.bed  -wa -wb -f 0.5 -F 0.5 -e > intersect.bed 
bedtools merge intersect.bed -d 150 



#### visualization 
mkdir 10deep_tools_result 
for down_sample in $(ls 7downsampled_replicates_results/*.bam); do deeptools -noralizeUsingRPKM -extendReads 150 wig ;done

### peak annotations (PeakAnnotator,annotatePeaks.pl,Peakannotation)
mkdir 11peak_annotation_results
for peak in $(merged_replciated_peak/); do R script $peak -o 11peak_annotation_result/; done  

#### motif anlysis
mkdir 12Motif_enrich

extend peak to 300bp or 500bp
MEME-ChIP pipeline fasta in peaks compared to known motif
findMotifsGenome.pl
 
## ATAC-Seq Footprinting
mkdir 13Motif_footprint
pyDNase dnase_average_profile.py -A

## Defining High-Confidence Target Sites for Transcription Factors
fimo 

## TF-target gene netowrk





#########################################
## RNA-Seq pipeline
#########################################
hisat v2.0.4 with parameter -dta-cufflinks

### visualization
Deeptools2 -noralizeUsingRPKM


## references
# Chromatin Accessibility Landscape in Human Early Embryos and Its Association with Evolution
# Evolutionary flexibility in flooding response circuitry in angiosperms
# Profiling of Accessible Chromatin Regions across Multiple Plant Species and Cell Types Reveals Common Gene Regulatory Principles and New Control Modules



## SAM file description
# AS:i:-1:比对得分，-0.6-0.6L (l为序列长度)
# XS: 对于多个比对，最优比对的得分值(过滤多个比对位置的reads) e.g. cd 2raw_alignement_results; samtools view Mt_ATAC_C1.sam | grep "SRR8763747.10000024"
# XN:i: 不确定碱基的数目
# XM:i: 错配数目(用于过滤没有比对上的reads)
# XO:i: gap的数目
# XG:i:<N>: gap延伸数目
# NM:i:<N> 最小编辑距离
# YF:Z:<S> 为什么被过滤掉
# YT:Z:<S> UU:不属于配对;CP: pair aligned concordantly; DP indicates the read was part of a pair and the pair aligned discordantly. Value of UP indicates the read was part of a pair but the pair failed to aligned either concordantly or discordantly.
#MD:Z:<S> 错配的reference



















