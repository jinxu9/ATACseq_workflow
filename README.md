##	work flow for ATACseq data analysis 
Version : 1 
Description: This is not a pipeline. Just notes for people who doesn't familar with ATACseq data analysis to catch up. 

Data 2017-08-14 
Author : jinxu9@stanford.edu 


Before going to each step, it's better if you already know basic command lines and the most frequentely used bioinformation tools , such as :

1. awk 

2. sed 

3. samtools 

4. bedtools 




## Step 1 : Make peak list 
	1. Merge peaks for different samples / experimental condition/ cell types to get a union set of peaks. 
```
mergeZINBAPeaks EC_peaks_FDR3En5  Brain_peaks.narrowPeak.FDR3En5 Liver_peaks.narrowPeak.FDR3En5  Lung_peaks.narrowPeak.FDR3En5 Skin_peaks.narrowPeak.FDR3En5 >EC_peaks_FDR3En5.bed  

```
	2. Merge samples to get more representative peaks for cell type / condition etc. 
```
	# skin
	samtools merge   Skin_EC_merge.sort.rmdup.bam ../12-skin/Mapping/12-skin.pe.q10.sort.rmdup.bam ../8-skin/Mapping/8-skin.pe.q10.sort.rmdup.bam ../skin/skin_CD31_43K.pe.q10.sort.bam.rmdup.bam & 

	# Liver
	samtools merge  Liver_EC_merge.sort.rmdup.bam ../15-liver/Mapping/15-liver.pe.q10.sort.rmdup.bam ../16-liver/Mapping/16-liver.pe.q10.sort.rmdup.bam ../7-liver/Mapping/7-liver.pe.q10.sort.rmdup.bam & 
	#Brain
	samtools merge  Brain_EC_merge.sort.rmdup.bam ../5-brain/Mapping/5-brain.pe.q10.sort.rmdup.bam ../9-brain/Mapping/9-brain.pe.q10.sort.rmdup.bam ../Brain_Ad2_10/Mapping/Brain_Ad2_10.pe.q10.sort.rmdup.bam ../Brain_Ad2_5/Mapping/Brain_Ad2_5.pe.q10.sort.rmdup.bam &  
	# Lung
	samtools merge  Lung_EC_merge.sort.rmdup.bam ../6-lung/Mapping/6-lung.pe.q10.sort.rmdup.bam ../Lung_Ad2_11/Mapping/Lung_Ad2_11.pe.q10.sort.rmdup.bam ../Lung_Ad2_6/Mapping/Lung_Ad2_6.pe.q10.sort.rmdup.bam  &

```
	Peak calling again :
```
	for file in `ls *_merge.sort.rmdup.bam`
	do
	name=`echo $file |awk -F_ '{print $1}'`
	macs2  callpeak -t  $file -f BAM  -g mm --outdir ./  -q 0.01 -n $name   --nomodel  --shift 0  &  
	done
```
## Step 2: Count peak signal 
	Count peaks for each sample and make a final matrix 
```
	bed="ECs_merge_peaks.bed"
	bedtools multicov  -bams  sample1.bam sample2.bam .sample3.bam  -bed  \$bed  > \$bed.AddPeakCounts 

```	

## Step 3: Global profile from peak signal  (Read more from the manual of DESeq2, an R package) 
	1. Normalization by sequencing depth / quartile / enrichment score etc .
	2. Overview with PCA and clustering 
```
	See code/DESeq2_PCA.R 
```

## Step 4: Differential peak call (Read more from the manual of DESeq2, an R package)
	Detect the differential peaks between different condition/ cell type etc. 
```
	See code/DESeq2_DEP.R 
```

## Step 5: Visualization on UCSC genome browser (Read more from https://genome.ucsc.edu/goldenpath/help/customTrack.html )
	1. Convert the alignment files into bigwig format, which can be uploaded and visualized from UCSC genome browser. 
```
	/seq/ATAC-seq/Code/bam2bed_shift.pl $file 
	genomeCoverageBed -bg -split -i \$pref.bed -g  /home/jinxu/DB/mmu9/mm9_UCSC_genome/mm9_all.chrsize > \$pref.bedGraph
	norm_bedGraph.pl  \$pref.bedGraph  \$pref.norm.bedGraph 1>log 2>err 
	bedGraphToBigWig  \$pref.norm.bedGraph  /home/jinxu/DB/mmu9/mm9_UCSC_genome/mm9_all.chrsize  \$pref.norm.bw

```
	2. Upload to any  locations can be access using URL
	
	3. Configure and save in UCSC browser. 
```
	https://genome.ucsc.edu/goldenpath/help/customTrack.html
```

## Step 6: Annotation of [differetial] peaks (Add later) 
	1. Annotation by GREAT
	2. Motif enrichment test 
	3. TF foot printing test
