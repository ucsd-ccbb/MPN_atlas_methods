# Germline Variants Analysis
* **Analysis by Guorong Xu, Ph.D. and Adam Mark, M.S.**
* **Methods document by Amanda Birmingham**
* **Center for Computational Biology & Bioinformatics, University of California, San Diego**
 
Each sample is represented by one pair of forward- and reverse-read fastq files, and the following steps are run once for each sample.  The `$file_name` variable is set to the base part of the fastq pair's file names; for example, for a pair of fastq files named `mysample1_R1.fq` and `mysample1_R2.fq`, `$file_name` = `mysample1`.  `$rootdir` is the directory on the system where the software tools are installed, `$resource` is the directory where tool resources such as indices and genome files are stored, and `$workspace` is the working directory for the analysis. These steps are bash scripts designed to be run from the command line on an Amazon Linux (AMI 2017.09.1) instance.  All resource files referenced below were downloaded from ftp://ftp.broadinstitute.org/bundle/b37/ .  

1. Align two paired-end fastq files with BWA against the human genome, convert output to a stream with samblaster, and write to a bam file:

		groupID="20170718"
		phenotype=$file_name
		numThreads=16
		bwa=$rootdir/bwa/bwa-0.7.12/bwa
		samblaster=$rootdir/samblaster/samblaster # version 0.1.21
		samtools=$rootdir/samtools/samtools-1.1/samtools
		genomeSeq=$resource/bwa/human_g1k_v37.fasta
		$bwa mem -M -t $numThreads -R '@RG\tID:1\tPL:illumina\tPU:'$groupID'\tSM:' \
			$phenotype -v 1 $genomeSeq $workspace/$file_name'_R1.fq' \
			$workspace/$file_name'_R2.fq' |$samblaster | 
			$samtools view -Sb - > $workspace/$file_name.bam
			
2. Sort and index the bam file with sambamba:

		numThreads=16
		sambamba=$rootdir/sambamba/0.4.7/bin/sambamba
		
		tmpDir=$workspace/temp
		
		echo "Sort is processing: $workspace/$file_name.bam"
		
		### Sort BAM File ####
		$sambamba sort -t $numThreads -m 10G --tmpdir $tmpDir \
			-o $workspace/$file_name.sort.bam $workspace/$file_name.bam
		
		$sambamba index -t $numThreads $workspace/$file_name.sort.bam \
			$workspace/$file_name.sort.bam.bai
		
2. Mark duplicates in the sorted bam with Picard and index it with sambamba:

		echo "running picard..."
		## Picard
		java -jar -Djava.io.tmpdir=$workspace/temp -Xms250m -Xmx20g \
			$rootdir/picard-1.96/MarkDuplicates.jar \
			INPUT=$workspace/$file_name.sort.bam OUTPUT=$workspace/$file_name.dedup.bam \
			METRICS_FILE=$workspace/$file_name.metrics.txt AS=true \
			VALIDATION_STRINGENCY=LENIENT
		
		echo "$(date): running sambamba..."
		## Sambamba
		$rootdir/sambamba/0.4.7/bin/sambamba index -t 16 \
			$workspace/$file_name.dedup.bam $workspace/$file_name.dedup.bam.bai


3. Run the following steps once for each chromosome (the autosomes, X, Y and MT), setting the `$chromosome` variable to to the current chromsome number or letter each time.

	1. Split the bam file up by chromosome using samtools and index it with sambamba:

			numThreads=1
			samtools=$rootdir/samtools/samtools-1.1/samtools
			sambamba=$rootdir/sambamba/0.4.7/bin/sambamba
			
			outputBAM=$file_name.$chromosome.bam
			
			### Split BAM####
			$samtools view -b $workspace/$file_name.dedup.bam \
				$chromosome > $workspace/$outputBAM
			$sambamba index -t $numThreads $workspace/$outputBAM \
				$workspace/$outputBAM.bai

	1. Perform local realignment of indels and correct systematic errors in base quality scores using bedtools 2.26.0 and GATK's sequence data processing tools (RealigerTargetCreator, IndelRealigner, BaseRecalibrator, and PrintReads): 

			echo "running bedtools..."
			$rootdir/bedtools2/bin/bedtools genomecov -split -ibam \
				$file_name.$chromosome.bam -bga \
				-g $resource/bwa/human_g1k_v37.fasta.fai \
				-max 70001 > $file_name.$chromosome.bed
			
			echo "$running RealignerTargetCreator..."
			java -Djava.io.tmpdir=$workspace/temp -Xmx20g \
				-jar $rootdir/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar \
				-T RealignerTargetCreator \
				-R $resource/bwa/human_g1k_v37.fasta \
				-I $file_name.$chromosome.bam \
				-known $resource/1000G_phase1.indels.b37.vcf \
				-known $resource/Mills_and_1000G_gold_standard.indels.b37.vcf \
				-o $workspace/intervals.$chromosome.interval_list -L $chromosome
			
			echo "$running IndelRealigner..."
			java -Djava.io.tmpdir=$workspace/temp -Xmx20g \
				-jar $rootdir/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar \
				-T IndelRealigner -R $resource/bwa/human_g1k_v37.fasta \
				-I $workspace/$file_name.$chromosome.bam \
				-known $resource/1000G_phase1.indels.b37.vcf \
				-known $resource/Mills_and_1000G_gold_standard.indels.b37.vcf \
				-targetIntervals $workspace/intervals.$chromosome.interval_list \
				-o $workspace/$file_name.realign.$chromosome.bam -L $chromosome
			
			echo "running BaseRecalibrator..."
			java -Djava.io.tmpdir=$workspace/temp -Xmx20g \
				-jar $rootdir/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar \
				-T BaseRecalibrator -R $resource/bwa/human_g1k_v37.fasta \
				-I $workspace/$file_name.$chromosome.bam \
				-knownSites $resource/dbsnp_138.b37.vcf \
				-knownSites $resource/hapmap_3.3.b37.vcf \
				-o $workspace/$file_name.$chromosome.grp \
				-dcov 2 -L $workspace/$file_name.$chromosome.bed
			
			echo "running PrintReads..."
			java -Djava.io.tmpdir=$workspace/temp -Xmx20g \
				-jar $rootdir/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar \
				-T PrintReads -R $resource/bwa/human_g1k_v37.fasta \
				-I $workspace/$file_name.realign.$chromosome.bam \
				-BQSR $workspace/$file_name.$chromosome.grp \
				-o $workspace/$file_name.final.$chromosome.bam -rf BadCigar \
				-L $chromosome --no_pg_tag
		
   1. Run GATK's HaplotypeCaller to identify variants

			echo "running gatk..."
			java -Xms454m -Xmx20g -XX:+UseSerialGC -Djava.io.tmpdir=$workspace/temp \
				-jar $rootdir/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar \
				-T HaplotypeCaller \
				-R $resource/bwa/human_g1k_v37.fasta \
				-I $workspace/$file_name.final.$chromosome.bam \
				-L $workspace/$file_name.$chromosome.bed \
				--out $workspace/$file_name.$chromosome.vcf.gz \
				--annotation BaseQualityRankSumTest \
				--annotation FisherStrand \
				--annotation GCContent \
				--annotation HaplotypeScore \
				--annotation HomopolymerRun \
				--annotation MappingQualityRankSumTest \
				--annotation MappingQualityZero \
				--annotation QualByDepth \
				--annotation ReadPosRankSumTest \
				--annotation RMSMappingQuality \
				--annotation DepthPerAlleleBySample \
				--annotation Coverage \
				--annotation ClippingRankSumTest \
				--standard_min_confidence_threshold_for_calling 30.0 \
				--standard_min_confidence_threshold_for_emitting 30.0 \
				--dbsnp $resource/dbsnp_132_b37.leftAligned.vcf
				
4. Merge the vcf files for all chromosomes into a single vcf using tabix and vcftools:		

		bgzip=$rootdir/tabix-0.2.6/bgzip
		tabix=$rootdir/tabix-0.2.6/tabix
		vcfConcat=$rootdir/vcftools_0.1.12b/bin/vcf-concat
		vcfSort=$rootdir/vcftools_0.1.12b/bin/vcf-sort
		
		tmpDir=$workspace/tmp
		
		for vcf_file in $workspace/$file_name.*.vcf; 
		do
		    file_list=$file_list$vcf_file" "
		done
		
		echo $file_name " is processing..."
		
		$vcfConcat $file_list | $bgzip -c > $workspace/$file_name.vcf.gz
		$vcfSort -t $tmpDir $workspace/$file_name.vcf.gz > $workspace/$file_name.raw.vcf
		
5. Compress and index the merged vcf file using bgzip and bcftools:	
	
		bcftools=$rootdir/bcftools-1.5/bcftools
		
		$bgzip -c $workspace/$file_name.raw.vcf > $workspace/$file_name.raw.vcf.gz
		$bcftools index $workspace/$file_name.raw.vcf.gz
	
There is now one vcf file for each sample.  Merge these per-sample vcfs into a single vcf that contains information for all samples, using bcftools:

	$bcftools merge [list the raw.vcf.gz files for each sample here]-o merged.vcf

Use ANNOVAR 2017Jun01 to annotate the multi-sample vcf, with `$inputfolder` set to the directory in which ANNOVAR inputs should be placed:

	inputfile="merged"
	
	perl $rootdir/ANNOVAR/convert2annovar.pl -format vcf4 -allsample -withfreq \
		-includeinfo $inputfolder"/"$inputfile".vcf" > $inputfolder"/"$inputfile".avinput"
	
	perl $rootdir/ANNOVAR/table_annovar.pl $inputfolder"/"$inputfile".avinput" \
		$rootdir/ANNOVAR/humandb/ -buildver hg19 -otherinfo \
		-out $inputfolder"/"$inputfile -remove \
		-protocol refGene,knownGene,esp6500siv2_all.1000g2015aug_all,snp138,\
		ljb26_all,cg46,cg69,popfreq_all_20150413,clinvar_20170130,cosmic70,nci60 \
		-operation g,g,f,f,f,f,f,f,f,f,f,f -nastring . -csvout
	
	perl $rootdir/ANNOVAR/table_annovar.pl --vcfinput $inputfolder"/"$inputfile".vcf" \
		$rootdir/ANNOVAR/humandb/ -buildver hg19 -otherinfo \
		-out $inputfolder"/"$inputfile -remove \
		-protocol refGene,knownGene,esp6500siv2_all,1000g2015aug_all,snp138,\
		ljb26_all,cg46,cg69,popfreq_all_20150413,clinvar_20170130,cosmic70,nci60 \
		-operation g,g,f,f,f,f,f,f,f,f,f,f -nastring . 
		
The ANNOVAR results are filtered in R version 3.4.1 using a script that takes six positional arguments.  In order, these are:

1. A path to the input multi-sample vcf file created by the above code.  Used only if the sixth parameter (load) is FALSE.
2. The name of the output csv file to create.  Also, if the sixth argument (load) is true, this is treated as the base name of a `.Rdata` file containing pre-processed data to load.  
3. A string indicating which genome to use; set to "b37" for the described analyses.  Used only if the sixth argument (load) is FALSE.
4. A TRUE or FALSE indicating whether to filter the annotations, retaining only those categorized as downstream, upstream, intergenic, intronic, or ncRNA\_intronic, and categorized as synonymous single nucleotide variants, and having ExAC\_ALL frequency < 0.05 and 1000G\_ALL frequency < 0.05.
5. A string indicating the variant caller used, either "gatk", "vardict", or "mutect".
6. A TRUE or FALSE indicating whether to load the annotations from a pre-existing `.Rdata` file.

The script is shown below:

	library(VariantAnnotation) # version 1.24.5
	library(Rsamtools) # version 1.30.0
	library(pbapply) # version 1.3-4
	library(data.table) # version 1.10.4-3
	library(magrittr) # version 1.5
	
	args <- commandArgs(TRUE)
	vcf_path <- args[1]
	csv_output <- args[2]
	genome <- args[3]
	filter <- args[4]
	var_caller <- args[5]
	load <- args[6]
	
	.collapse <- function (...) {
	  paste(unlist(list(...)), sep = ",", collapse = ",")
	}
	
	#filter vcf object by annotation fields
	filter_vcf_object <- function(annotations, filter=filter, var_caller=var_caller){         
	  if(filter==TRUE){
	    annotations <- annotations[!unlist(info(annotations)$Func.refGene) %in% c(
	    	"downstream", "upstream", "intergenic", "intronic", "ncRNA_intronic")]
	    annotations <- annotations[!unlist(info(annotations)$Func.knownGene) %in% c(
	    	"downstream", "upstream", "intergenic", "intronic", "ncRNA_intronic")]
	    annotations <- annotations[!unlist(info(annotations)$ExonicFunc.knownGene) 
	    	%in% c("synonymous_SNV")]
	    annotations <- annotations[unlist(info(annotations)$ExAC_ALL) < 0.05]
	    annotations <- annotations[unlist(info(annotations)$`1000G_ALL`) < 0.05]
	  }
	  ID <- paste0(seqnames(annotations), ":", start(annotations), "_", 
	  	as.character(ref(annotations)), "/", sapply(alt(annotations), as.character))
	  if(var_caller %in% c("gatk", "vardict")){
	    AD <- geno(annotations)$AD
	    AD.df <- data.frame(lapply(1:ncol(AD), function(i) sapply(AD[,i], 
	    	function(j) round(j[2]/(j[1]+j[2]), 2))))
	    AD <- lapply(data.table(AD), sapply, .collapse) %>% data.frame
	    GT <- cbind(AD, AD.df)
	  }
	  if(var_caller=="mutect"){
	    AD <- geno(annotations)$AD
	    AD.df <- geno(annotations)$FA
	    AD <- lapply(data.table(AD), sapply, .collapse) %>% data.frame(
	    	check.names = FALSE)
	    GT <- cbind(AD, AD.df)
	  }
	  colnames(GT) <- c(paste0(colnames(AD), "_AD"), paste0(colnames(AD), "_AF"))
	  annotations <- cbind(ID=ID, info(annotations), data.table(GT))
	  annotations <- lapply(annotations, sapply, .collapse) %>% data.frame
	  names(annotations)[1:12] <- c("ID", "Gene", "Consequence", "AA_Change", 
	  	"Function_UCSC", "Function_refGene", "dbSNP", "CADD", "ExAC",
	                                "1000G", "COSMIC_70", "NCI-60")
	  names(annotations) <- gsub("X", "", names(annotations))
	  data.table(annotations, check.names=FALSE)
	}
	
	#read large annovar vcf into R                          
	yield_annovar <- function(vcf_path, genome, rdata_file=rdata_file, save_rdata){
	  param <- ScanVcfParam(geno=c("AD", "FA"),
	                        info=c("Gene.knownGene",
	                               "ExonicFunc.knownGene",
	                               "AAChange.refGene",
	                               "Func.knownGene",
	                               "Func.refGene",
	                               "snp138",
	                               "CADD_phred",
	                               "ExAC_ALL",
	                               "1000G_ALL",
	                               "cosmic70",
	                               "nci60"))
	  tab <- TabixFile(vcf_path, yieldSize=100000)
	  open(tab)
	  file_list <- list()
	  i=1
	  while (nrow(file_chunk <- readVcf(tab, genome=genome, param=param))){
	    file_list[[i]] <- file_chunk
	    cat("vcf chunk:", i, "\n")
	    i=i+1}
	  close(tab)
	  save(file_list, file=rdata_file)
	  file_list
	}    
	
	write.table_with_header <- function(x, file, header, ...){
	  cat(header, '\n',  file = file)
	  write.table(x, file, append = T, sep="\t", quote=FALSE, row.names=FALSE)
	}
	
	header=
	"Field Descriptions
	ID: Chrom:Position_ReferenceAllele/AlternateAllele
	Gene: Gene symbol
	Consequence: Variant consequence
	AA_Change: Amino acid change, if applicable
	Function_UCSC: location of variant, from UCSC gene annotation
	Function_refGene: location of variant, from RefGene annotation
	CADD: CADD functional impact score
	ExAC: ExAC population allele frequency
	1000G: 1000 Genomes allele frequency
	COSMIC_70: COSMIC (Catalogue Of Somatic Mutations In Cancer) ID if 
	present in database
	NCI-60: Variant present in NCI-60 human tumor cell lines
	Column data for samples contains counts for the reference and 
	alternate alleles, and alternate allele frequency (ref, alt, AF), 
	in the order listed"
	
	if(load==TRUE){
	  cat("Loading RData file\n")
	  load(gsub("txt", "RData", csv_output))
	  cat("Filtering Variants\n")
	  filt_list <- pblapply(file_list, filter_vcf_object, filter=filter, 
	  	var_caller=var_caller)
	} else {
	  annovar_list <- yield_annovar(vcf_path, genome=genome, 
	  	rdata_file=gsub("txt", "RData", csv_output), save_rdata=FALSE)
	  cat("Filtering Variants\n")
	  filt_list <- pblapply(annovar_list, filter_vcf_object, filter=filter, 
	  	var_caller=var_caller)}
	
	cat("Concatenating Variants\n")
	vcf <- rbindlist(filt_list)
	names(vcf) <- gsub("X", "", names(vcf))
	write.table_with_header(vcf, csv_output, header, sep="\t", quote=FALSE, 
		row.names=FALSE)

This script produces a tab-separated output file of filtered annotations.


## Software References

* **bwa**: Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. (2013) arXiv:1303.3997v2 [q-bio.GN].
* **samblaster**: Faust GG, Hall IM. SAMBLASTER: fast duplicate marking and structural variant read extraction. Bioinformatics. 2014 Sep 1;30(17):2503-5.
* **samtools**: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9.
* **sambamba**: Tarasov A, Vilella AJ, Cuppen E, Nijman IJ, Prins P. Sambamba: fast processing of NGS alignment formats. Bioinformatics. 2015 Jun 15;31(12):2032-4. 
* **picard**: http://broadinstitute.github.io/picard/ .
* **bedtools**: Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2.
* **GATK**: DePristo MA, Banks E, Poplin R, Garimella KV, Maguire JR, Hartl C, Philippakis AA, del Angel G, Rivas MA, Hanna M, McKenna A, Fennell TJ, Kernytsky AM, Sivachenko AY, Cibulskis K, Gabriel SB, Altshuler D, Daly MJ. A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nat Genet. 2011 May;43(5):491-8.
* **Tabix**: Li H. Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics. 2011 Mar 1;27(5):718-9.
* **vcftools**: Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE, Lunter G, Marth GT, Sherry ST, McVean G, Durbin R; 1000 Genomes Project Analysis Group. The variant call format and VCFtools. Bioinformatics. 2011 Aug 1;27(15):2156-8.
* **bcftools**: https://github.com/samtools/bcftools .
* **ANNOVAR**: Wang K, Li M, Hakonarson H. ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. Nucleic Acids Res. 2010 Sep;38(16):e164.
* **R**: R Core Team . R: A language and environment for statistical computing. (2016) R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/ .
* **VariantAnnotation**: Obenchain V, Lawrence M, Carey V, Gogarten S, Shannon P, Morgan M. VariantAnnotation: a Bioconductor package for exploration and annotation of genetic variants. Bioinformatics. 2014 Jul 15;30(14):2076-8.
* **Rsamtools**: Morgan M, Pagès H, Obenchain V and Hayden N. Rsamtools: Binary alignment (BAM), FASTA, variant call (BCF), and tabix file import. (2017) R package version 1.30.0, http://bioconductor.org/packages/release/bioc/html/Rsamtools.html .
* **pbapply**: https://cran.r-project.org/web/packages/pbapply/index.html .
* **data.table**: https://github.com/Rdatatable/data.table .
* **magrittr**: https://cran.r-project.org/web/packages/magrittr/index.html .