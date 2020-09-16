## ATAC-pipeline
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
    # for m in $(ls 2raw_alignement_results/*.sam); do echo $m;  samtools view -bSq 2 $m  >3uniq_mapped_results/$(basename ${m%.sam}).q2.bam;done

### convert sam to sorted bam
    for m in $(ls 3uniq_mapped_results/*.sam); do echo $m; picard SortSam  INPUT=$m OUTPUT=4uniq_mapped_results_sort/$(basename ${m%.sam}).sorted.bam SO=coordinate; done
    # for m in $(ls 3uniq_mapped_results/*.bam); do samtools sort $m 4uniq_mapped_results_sort/$(basename ${m%.q2.bam}).sorted; done


### Remove PCR duplicated reads (or Picard-tools)
    for m in $(ls 4uniq_mapped_results_sort/*.bam); do echo $m; samtools rmdup $m 5remove_duplicated_results/$(basename ${m%.bam}).dedup.bam; done
    #for m in $(ls 4uniq_mapped_results_sort/*.bam);do picard MarkDuplicates INPUT=$m OUTPUT=5remove_duplicated_results/$m     METRICS_FILE=5remove_duplicated_results/$m.mark_dup_metrix.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false


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


### fragment shift bam and fragement size distribution (deeptools:lose some information and the resulted bam can not be used in other softs,ATACseqQC,or bedpe)
    mkdir 7a_fragment_shift_result
    mkdir 7b_fragemnt_size_distribution

    #for m in $(ls 7downsampled_replicates_results/*.bam); do alignmentSieve --numberOfProcessors 25 --ATACshift --bam $m -o 7a_fragment_shift_result/$(basename $m).tmp.bam;done
    #for m in $(ls 7a_fragment_shift_result/*.tmp.bam); do samtools sort -@ 20 -O bam -o 7a_fragment_shift_result/$(basename ${m%.tmp.bam}) $m; samtools index -@ 20    7a_fragment_shift_result/$(basename ${m%.tmp.bam}); rm $m; done
    #samtools index -@ 8 sample1.shifted.bam

    for m in $(ls 7downsampled_replicates_results/*.bam)
    do
	echo $m
	if [ -e 7a_fragment_shift_result/$(basename $m) ]
	then	
		echo $m, ": reads shifting is existed, and skip!"
	 
	else
		Rscript Shift_bam.R -b $m -o 7a_fragment_shift_result/$(basename $m) -f 7b_fragemnt_size_distribution/ # the resulted bam are paired bam fies,memory cost ~10 x bam_size, and tiem costing!
	fi
    done

### merge replicates bams #(from original replcated bams or sampled replcated bams, used for footprint) 
    mkdir 7c_merge_replicates_bams 

    samples=$(ls 7a_fragment_shift_result/*.bam | xargs -n 1 basename | sed 's/\.uniq\.sorted\.dedup\.bam//g' | sed -r 's/[0-9]{1,}$//g'| sed s'/*\///g' | sort| uniq)

    for sample in $samples
    do 	#echo $sample #Mt_ATAC_C
	
	sample_repeats_files=$(ls 7a_fragment_shift_result/$sample*.bam)
	sample_repeats_num=$(ls 7a_fragment_shift_result/$sample*.bam |wc -l) #3
	
	if [ $sample_repeats_num == 1 ]
	then 
		echo $sample ": has no replicates, just copy the bam file!"	
		cp 7a_fragment_shift_result/$sample*.bam 7c_merge_replicates_bams/
	else
		samtools merge -@ 10  7c_merge_replicates_bams/$sample.bam 7a_fragment_shift_result/$sample*.bam
	fi
    done


### peak calling (homer,macs(macs2 callpeak -f BAMPE -q 0.001 --broad- cutoff 0.01 -g 17e9 --bw 300),hopspot:FDR < 0.01 with 8 reads,HMMRATAC)

## MACS2
    #make macs2
    #macs2 callpeak -f BAMPE -g hs --keep-dup all --cutoff-analysis -n sample1  -t sample1.shifted.bam --outdir macs2/sample1 2> macs2.log


    mkdir 8a_peak_calling_result_hmmratac #(erro with no tag )

    for m in $(ls 7a_fragment_shift_result/*.bam); do echo $m; samtools view -H $m| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > 8a_peak_calling_result_hmmratac/genome.info; java -jar ~/softwares/HMMRATAC_V1.2.10_exe.jar -b $m -i $m.bai -g 8a_peak_calling_result_hmmratac/genome.info -o 8a_peak_calling_result_hmmratac/$(basename ${m%.bam}) ;done
#---------------------------------------------------------------------------------------------------------------
    mkdir 8b_peak_calling_result_homer # need sam file

    for m in $(ls 7a_fragment_shift_result/*.bam); do samtools view -@ 5 -h -o ${m%.bam}.sam $m; done
    for m in $(ls 7a_fragment_shift_result/*.sam);do makeTagDirectory 8b_peak_calling_result_homer/$(basename ${m%.sam}) $m -single -format sam -mapq 2; done  # -single: chromsome number >100
    for m in $(ls 7a_fragment_shift_result/*.sam);do findPeaks 8b_peak_calling_result_homer/$(basename ${m%.sam})/ -region -minDist 150 -o 8b_peak_calling_result_homer/$(basename ${m%.sam});done  #(pc,2018,2018 pj)
    for m in $(ls 8b_peak_calling_result_homer/*.peak); do grep -v "#" $m  | awk 'OFS="\t"{print $2,$3,$4,$1,$8,$5}' > $m.bed;done

### Homer peak merger between replciates by bedtools to get union THS sites (or pooled data, and two replicated peak,remove scaffold)
    mkdir 9merged_replicated_peaks

    samples=$(ls 8b_peak_calling_result_homer/*.bed | xargs -n 1 basename | sed 's/\.uniq\.sorted\.dedup\.peak\.bed//g' | sed -r 's/[0-9]{1,}$//g'| sed s'/*\///g' | sort| uniq)

    for sample in $samples
    do 	#echo $sample #Mt_ATAC_C
	
	sample_repeats_files=$(ls 8b_peak_calling_result_homer/$sample*.bed)
	sample_repeats_num=$(ls 8b_peak_calling_result_homer/$sample*.bed|wc -l) #3
	
	if [ $sample_repeats_num == 1 ]
	then 
		echo $sample ": has no replicates, just copy the bam file!"	
		cp 8b_peak_calling_result_homer/$sample*.bed 9merged_replicated_peaks/
		grep -v "scaffold" 9merged_replicated_peaks/$sample*.bed > 9merged_replicated_peaks/$sample*.bed
	fi
	
	if [ $sample_repeats_num == 2 ]
	then
		repeat1_peak=$(ls 8b_peak_calling_result_homer/$sample*.bed| head -n 1)
		repeat2_peak=$(ls 8b_peak_calling_result_homer/$sample*.bed| tail -n 1)
		
		bedtools intersect -a $repeat1_peak -b $repeat2_peak  -wa -wb -f 0.5 -F 0.5 > intersect_rep1_rep2.bed
		cut -f1,2,3,4,5,6 intersect_rep1_rep2.bed > inter_rep1.bed
		cut -f7,8,9,10,11,12 intersect_rep1_rep2.bed > inter_rep2.bed
		cat inter_rep1.bed inter_rep2.bed | bedtools sort| bedtools merge -d 150 | grep -v "scaffold" > 9merged_replicated_peaks/$sample.bed
		rm intersect_rep1_rep2.bed inter_rep*.bed
	fi

	if [ $sample_repeats_num == 3 ]
	then
		repeat1_peak=$(ls 8b_peak_calling_result_homer/$sample*.bed | head -n 1)
		repeat2_peak=$(ls 8b_peak_calling_result_homer/$sample*.bed | sed -n 2p)
		repeat3_peak=$(ls 8b_peak_calling_result_homer/$sample*.bed | tail -n 1)
		
		bedtools intersect -a $repeat1_peak -b $repeat2_peak  -wa -wb -f 0.5 -F 0.5 > intersect_rep1_rep2.bed
		cut -f1,2,3,4,5,6 intersect_rep1_rep2.bed > inter_rep1_2.bed
		cut -f7,8,9,10,11,12 intersect_rep1_rep2.bed > inter_rep2_1.bed

		bedtools intersect -a $repeat2_peak -b $repeat3_peak  -wa -wb -f 0.5 -F 0.5 > intersect_rep2_rep3.bed
		cut -f1,2,3,4,5,6 intersect_rep2_rep3.bed > inter_rep2_3.bed
		cut -f7,8,9,10,11,12 intersect_rep2_rep3.bed > inter_rep3_2.bed
		
		bedtools intersect -a $repeat3_peak -b $repeat1_peak  -wa -wb -f 0.5 -F 0.5 > intersect_rep3_rep1.bed
		cut -f1,2,3,4,5,6 intersect_rep3_rep1.bed > inter_rep3_1.bed
		cut -f7,8,9,10,11,12 intersect_rep3_rep1.bed> inter_rep1_3.bed

		cat inter_rep1_2.bed inter_rep2_1.bed inter_rep2_3.bed inter_rep3_2.bed inter_rep3_1.bed  inter_rep1_3.bed | bedtools sort| bedtools merge -d 150 | grep -v "scaffold" > 9merged_replicated_peaks/$sample.bed
		rm intersect_rep1_rep2.bed intersect_rep2_rep3.bed intersect_rep3_rep1.bed inter_rep*.bed
	fi
	done

### MACS peak merge
    mkdir 9merged_replicated_peaks

    samples=$(ls 8macs/*.narrowPeak | xargs -n 1 basename | sed 's/\.narrowPeak//g' | sed -r 's/[0-9]{1,}$//g'| sed s'/-///g' | sort| uniq)

    for sample in $samples
    do 	#echo $sample #Mt_ATAC_ C
	
	sample_repeats_files=$(ls 8macs/$sample*.narrowPeak)
	sample_repeats_num=$(ls 8macs/$sample*.narrowPeak | wc -l) #count sample replicate number
	
	if [ $sample_repeats_num == 1 ]
	then 
		echo $sample ": has no replicates, just copy the bam file!"	
		cp 8macs/$sample*.narrowPeak 9merged_replicated_peaks/
		#grep -v "scaffold" 9merged_replicated_peaks/$sample*.narrowPeak > 9merged_replicated_peaks/$sample*.narrowPeak  # exclude scaffold peak
	fi
	
	if [ $sample_repeats_num == 2 ]
	then
		repeat1_peak=$(ls 8macs/$sample*.narrowPeak| head -n 1)
		repeat2_peak=$(ls 8macs/$sample*.narrowPeak| tail -n 1)
		
		bedtools intersect -a $repeat1_peak -b $repeat2_peak  -wa -wb -f 0.5 -F 0.5 > intersect_rep1_rep2.bed   # there may be a problem when a longer peak and a shorter peak are overlaped
		cut -f1,2,3,4,5,6,7,8,9,10 intersect_rep1_rep2.bed > inter_rep1.bed
		cut -f11,12,13,14,15,16,17,18,19,20 intersect_rep1_rep2.bed > inter_rep2.bed
		cat inter_rep1.bed inter_rep2.bed | bedtools sort| bedtools merge -d 150 | grep -v "scaffold" > 9merged_replicated_peaks/$sample.narrowPeak
		rm intersect_rep1_rep2.bed inter_rep*.bed
	fi

	if [ $sample_repeats_num == 3 ]
	then
		repeat1_peak=$(ls 8macs/$sample*.narrowPeak | head -n 1)
		repeat2_peak=$(ls 8macs/$sample*.narrowPeak | sed -n 2p)
		repeat3_peak=$(ls 8macs/$sample*.narrowPeak | tail -n 1)
		
		bedtools intersect -a $repeat1_peak -b $repeat2_peak  -wa -wb -f 0.5 -F 0.5 > intersect_rep1_rep2.bed
		cut -f1,2,3,4,5,6 intersect_rep1_rep2.bed > inter_rep1_2.bed
		cut -f7,8,9,10,11,12 intersect_rep1_rep2.bed > inter_rep2_1.bed

		bedtools intersect -a $repeat2_peak -b $repeat3_peak  -wa -wb -f 0.5 -F 0.5 > intersect_rep2_rep3.bed
		cut -f1,2,3,4,5,6 intersect_rep2_rep3.bed > inter_rep2_3.bed
		cut -f7,8,9,10,11,12 intersect_rep2_rep3.bed > inter_rep3_2.bed
		
		bedtools intersect -a $repeat3_peak -b $repeat1_peak  -wa -wb -f 0.5 -F 0.5 > intersect_rep3_rep1.bed
		cut -f1,2,3,4,5,6 intersect_rep3_rep1.bed > inter_rep3_1.bed
		cut -f7,8,9,10,11,12 intersect_rep3_rep1.bed> inter_rep1_3.bed

		cat inter_rep1_2.bed inter_rep2_1.bed inter_rep2_3.bed inter_rep3_2.bed inter_rep3_1.bed  inter_rep1_3.bed | bedtools sort| bedtools merge -d 150 | grep -v "scaffold" > 9merged_replicated_peaks/$sample.bed
		rm intersect_rep1_rep2.bed intersect_rep2_rep3.bed intersect_rep3_rep1.bed inter_rep*.bed
	fi
    done
### IDR merge peaks between samples
    for m in $(ls 8macs/*-1_peaks.narrowPeak); do prefix=$(basename ${m%-1_peaks.narrowPeak}); idr --samples 8macs/$prefix-1_peaks.narrowPeak 8macs/$prefix-3_peaks.narrowPeak --input-file-type narrowPeak --idr-threshold 0.05 --output-file 9.3merged_replicated_peaks_idr/$prefix.idr.peak.txt --plot; done

#### visualization 
    mkdir 10deep_tools_result 
    for down_sample in $(ls 7downsampled_replicates_results/*.bam); do deeptools -noralizeUsingRPKM -extendReads 150 wig ;done

### peak annotations (PeakAnnotator,annotatePeaks.pl,Peakannotation)
    mkdir 11peak_annotation_results
    for peak in $(merged_replciated_peak/); do R script $peak -o 11peak_annotation_result/; done  

### motif enrichment  anlysis (meme-chip pipeline,need resize peak from homer, et. al)
    mkdir 12Motif_enrich 
    TARGET_LENGTH=300

    for m in $(ls 9merged_replicated_peaks/*.bed); do echo $m; awk -vF=${TARGET_LENGTH} 'BEGIN{ OFS="\t"; }{ len=$3-$2; diff=F-len; flank=int(diff/2); upflank=downflank=flank; if (diff%2==1) { downflank++; }; print $1, $2-upflank, $3+downflank; }' $m | bedtools sort > ${m%.bed}.resize.bed;done # resize peak length to constent length

    for m in $(ls 9merged_replicated_peaks/*resize.bed);do bedtools getfasta -fi ~/genome/Medicago_truncatula/Medicago_truncatula.MedtrA17_4.0.31.dna.toplevel.fa -fo ${m%.resize.bed}.fa  -bed $m; done # get fasta sequences

    for m in $(ls 9merged_replicated_peaks/*.fa);do meme-chip -oc 12Motif_enrich/$(basename ${m%.fa}) -db ~/genome/Medicago_truncatula/Mtr_TF_binding_motifs.meme -meme-p 20 $m; done

#findMotifsGenome.pl
 
### ATAC-Seq Footprinting in peak (pyDNase (peak>100bp,100 MB reads,need merged bam and merged peak from replciates ?); HINT-ATAC(2019,Genome Biology))
    mkdir 13Motif_footprint
    for m in $(ls 7c_fragment_shift_result/*.bam); do wellington_footprints.py 9merged_replicated_peaks/Mt_ATAC_C.bed $m 13Motif_footprint/$(basename m{%.bam}) -o $(basename m{%.bam}) -p 10 -A; done

 
#rgt-hint footprinting --help

## Defining High-Confidence Target Sites for Transcription Factors
    fimo 

## TF-target gene netowrk
    bedtoos gene.bed > gene_up2kb_down_2.5kb.bed
    bedtools intersect gene_up2kb_down_2.5kb.bed footprint.bed > focus.footprint.txt

    motif in footprint with -2k,gene body 

