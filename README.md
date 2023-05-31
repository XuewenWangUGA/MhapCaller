# MhapCaller
Macrohaplotype caller for STR, SNP and InDel from NGS long read sequencing data, especially for PacBio HiFi reads
The MhapCaller calls targeted STRs/SSRs, SNPs and InDels simutanunously in a row from each single NGS read and clusters the variants into phased haplotype string. MhapCaller is the best tool to analysis all haplotypic variants from genetically inherited DNA. It suits for diploid, polyploid, and DNA mixtures from many individuals, e.g. DNA forensics. MhapCaller is programmed in Java with parallele computing enabled so it can run in any computing platforms. 

The MhapCaller is also integrated into a pipeline with Python3. The pipeline takes the fastq reads as input, align reads to the genome, sort and index alignment, and then run MhapCaller. Additional tools are also available as util tools.

The whole package is also available in a Singularity container.

 MhapCaller runs fast. It takes ~ 2 mins for a merged data from four full runs of PacBio smrtcell HiFi reads on a Linux machine on Cluster with 12 threads,and maximumm allowed 120G RAM. 

## Latest version
V0.3

## Usage
for help: 

`java -jar MHapCaller.jar`

## Command
MHapCaller v0.3

usage: `java -jar -Xmx100G MHapCaller.jar [options]`

 **-a,--strAnchorFile** <arg>  &nbsp; configure file of anchors for STRs in
                                 tabular plain text
 
**-c,--MinReadCount** &nbsp; integer, minimum read count to report a
                                 macrohaplotype, default [2]
 
 **-d,--inDelBedFile** <arg> &nbsp; configure file of InDels in BED format in
                                 tabular plain text
 
 **-i,--input** <arg> &nbsp; input BAM file with a path
 
 **-l,--homopolymerLen** &nbsp; integer, homopolymer length threshold to
                                 skip call variant,default [4]
 
 **-m,--MaxMismatches**  &nbsp; integer, maximum distance of edit
                                 distance in anchor match, default [1]
 
 **-n,--snpPanelPosFile** <arg> &nbsp; configure file of SNPs in tabular plain
                                 text
 
 **-o,--output** <arg> &nbsp; output file name for macrohaplotype
                                 result
 
 **-p,--MinReadProportion** &nbsp; double, minimum read proportion of all
                                 reads for each locus to report a macrohaplotype, default
                                 [0.001]
 
 **-q,--LowQualThreshold** &nbsp; integer, phred-scale quality threshold
                                 for variant bases, default [15]
 
 **-r,--refGenomeFastaFile** <arg> &nbsp; Genome reference sequence file in .fasta
                                 format
 
 **-t,--ThreadNumber** &nbsp; integer, the number of computing threads,
                                 default [12]

 
 ## Example 
 
    java -jar -Xmx120G MHapCaller0.3.jar -i /mnt/data0/xuewen/macrohaplotype/hg002/Q40hifi/hg002.8kampl.Q40.fastq.gz_GRCh38.bam -o hg002.8kampl.Q40_GRCh38.tsv -a /mnt/data0/xuewen/macrohaplotype/scripts/MHvar_v0.3/CODISSTR_anchor.XW.config_v0.2.txt -d /mnt/data0/xuewen/macrohaplotype/scripts/MHvar_v0.3/MHindels_v0.2.bed -n /mnt/data0/xuewen/macrohaplotype/scripts/MHvar_v0.3/MHsnps.pos.txt -r /mnt/data0/xuewen/hg38giab/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta -l 2 -m 1 -q 15 -c 100 -p 0.01 -t 12

## Preprocess pipeline for MhapCaller
 This pipeline will preprocess the orginal PacBio hifi reads in BAM, quality control, umi-analysis (optional), statistica summary, read-reference alignment, sort, and index. The results are read-reference alignment in BAM format and BAM index.
 
 `MH_umi_dedup_map_pipe.sh`
 
 

## Input of MhapCaller
 1. a read-reference alignment file in [BAM format](https://en.wikipedia.org/wiki/SAM_(file_format))

 2. Configure files of tartgeted variant positional information for each of STRs, SNPs, and InDels
 
 ## Output of MhapCaller

Macrohaplotypes and supported reads count in tabular text file .tsv. 
E.g., the two macrohaplotypes around FBI's CODIS loci D2S441 in 8 kb PacBio HiFi reads of benchmark reference hg002. Locus details are at [FBI](https://www.fbi.gov/how-we-can-help-you/dna-fingerprint-act-of-2005-expungement-policy/codis-and-ndis-fact-sheet)
     
The first line shows the total reads for this locus D2441. The position lines present the coordinate of each targeted variant site on chr2 in human genome Hg38. Then the macrohaplotypes, including three parts of SNPs, InDels and STRs separated by ";".
     
    #Total hapVar: 	D2S441	27694	
    #Markername	Counts	HapVarLen	hapVar(s)
    #			Position:68011922,68011990,68012066,68012109,68012198,68012335,68012415,68012463,68012596,68012616,68012754,68012776,68012877,68012888,68013086,68013124,68013388,68013405,68013406,68013664,68013692,68013770,68013988,68014065,68014213,68014245,68014581,68014700,68015087,68015100,68015130,68015131,68015330,68015388,68015498,68015787,68015976,68016033,68016045,68016108,68016576,68016583,68017126,68017147,68017242,68017364,68017504,68017573,68017732,68017739,68018169,68018176,68018322,68018330,68012435,68012478,68012608,68012694,68012722,68012762,68012781,68011917
    
    D2S441	8646	136	
    G,A,A,A,G,C,T,A,T,C,C,T,G,T,T,G,C,C,G,C,A,A,G,T,A,T,A,G,A,A,C,G,T,G,C,G,G,G,T,C,C,G,C,G,G,A,A,C,C,C,T,A,A,G;TA,CACATATATA,CA,CA,CACATATATA,CATATATA,ATGT;TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA

    D2S441	7467	152	
    G,A,A,A,G,C,T,A,T,T,C,T,A,T,C,G,C,C,T,C,A,G,G,T,G,T,A,G,A,A,C,G,C,G,C,G,G,G,T,C,C,G,C,G,A,A,A,C,C,C,T,G,A,G;T|TAC,CACATATAT|C,CA,CA,CACATATATA,CATATAT,T;TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATTTATCTATCTA
    
   ## Visualization 
   The macrohaplotypes can be visualized or edited in a text processor, Microsoft excels, and spread sheet.
   
   The graphically view and compare the macrohaplotypes, please use our graphic tool: [USAT](https://github.com/XuewenWangUGA/USAT)
   
   ## Citation
   coming soon





