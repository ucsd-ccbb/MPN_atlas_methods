# RNA Editing Analysis
* **Analysis by Adam Mark, M.S.**
* **Methods document by Amanda Birmingham**
* **Center for Computational Biology & Bioinformatics, University of California, San Diego**
 
Each sample is represented by one pair of forward- and reverse-read fastq files.  `$rootdir` is the directory on the system where the software tools are installed, `$resource` is the directory where tool resources such as indices and genome files are stored, and `$workspace` is the working directory for the analysis. The below steps are bash scripts designed to be run from the command line on Ubuntu Linux 14.04 or higher.

Begin by setting up the variables:
	
	# Environment
	cpu="$(cat /proc/cpuinfo | grep processor | wc -l)
	genome=$resource/sequences/Hsapiens/hg19/ucsc.hg19.fasta
	# genome file downloaded from 
	# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
	
	# star
	genomeDir=$resource/star_index/hg19
	gtf=$resource/star_index/hg19/Homo_sapiens.GRCh37.87.gtf
	# gtf file downloaded from
	# ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/
	#   Homo_sapiens.GRCh37.87.gtf.gz
	ergid=1234
	rgpl=ILLUMINA
	rglb=TrueSeq
	
	# RNA variant calling
	dbsnp=$resource/variation/dbsnp_138.hg19.vcf.gz
	# dbsnp file downloaded from
	# ftp://ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz
	mills=$resource/variation/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
	# mills file downloaded from 
	# ftp://ftp.broadinstitute.org/bundle/hg19/
	#	Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
	g1000=$resource/variation/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
	# g1000 file downloaded from 
	# ftp://ftp.broadinstitute.org/bundle/hg19/
	# 	1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
	
	# SNPiR
	SNPiR=$rootdir/SNPiR
	RepeatMasker=$resource/annotation/RepeatMasker.hg19.bed
	# RepeatMasker and gene_annotation files produced according to 
	# instructions from step 1 at
	# https://github.com/rpiskol/SNPiR/blob/4ccef45038942ad37adf2693774b99408c05b1a6/README
	AluRegions=$resource/annotation/Alu.hg19.bed
	# AluRegions file is a subset of the RepeatMasker file produced by running
	# grep Alu $resource/annotation/RepeatMasker.hg19.bed > $resource/annotation/Alu.hg19.bed
	gene_annotation=$resource/SNPiR/db/SNPiR_annotation.hg19
	radar=$resource/db/Human_AG_all_hg19_v2.txt
	# radar file downloaded from 
	# http://lilab.stanford.edu/GokulR/database/Human_AG_all_hg19_v2.txt
	darned=$resource/db/hg19.txt
	# darned file downloaded from
	# https://darned.ucc.ie/static/downloads/hg19.txt
	rnaEditDB=$resource/db/rnaEditDB.txt
	# rnaEdit file produced by merging and sorting 
	# radar and darned files above
	
	# bedtools
	intersectBed=$rootdir/bedtools2/bin/intersectBed
	
Run the following steps are run once for each sample.  The `$sample` variable is set to the base part of the fastq pair's file names; for example, for a pair of fastq files named `mysample1_R1.fq` and `mysample1_R2.fq`, `$sample` = `mysample1`.  

1. Set the sample-specific variables:

		# star
		dna_sample=$(echo $sample | cut -d'-' -f1)
		fastq1="$sample"_R1.fastq.gz
		fastq2="$sample"_R2.fastq.gz
		bam=$workspace/"$sample"_Aligned.sortedByCoord.out.bam
		dedup_bam=$workspace/"$sample"_dedup.bam
		dedup_sorted_bam=$workspace/"$sample"_dedup_sorted.bam
		split_bam=$workspace/"$sample"_split.bam
		realigned_bam=$workspace/"$sample"_realigned.bam
		bqsr_bam=$workspace/"$sample"_bqsr.bam
		
		# Variant calling
		raw_vcf=$workspace/"$sample".vcf
		gatk_tmp_vcf=$workspace/"$sample"_gatk_tmp.vcf
		gatk_filt_vcf=$workspace/"$sample"_gatk_filt.vcf
			
2. Align the forward and reverse reads using STAR 2.5.2b and index the resulting bam file with samtools 1.3:

	 	# Alignment
	 	echo "Creating RNA BAM"
	 	cd $workspace;
	 	STAR \
	 		--genomeDir $genomeDir \
	 		--readFilesIn $workspace/$fastq1 $workspace/$fastq2 \
	 		--readFilesCommand zcat \
	 		--twopassMode Basic \
	 		--quantMode GeneCounts \
	 		--sjdbGTFfile $gtf \
	 		--alignIntronMax 200000 \
	 		--alignMatesGapMax 200000 \
	 		--outFilterMismatchNoverLmax 0.04 \
	 		--outFilterIntronMotifs RemoveNoncanonical \
	 		--outSAMtype BAM SortedByCoordinate \
	 		--outBAMsortingThreadN 16 \
	 		--outSAMunmapped Within \
	 		--outSAMattrRGline ID:$rgid SM:"$sample" PL:$rgpl LB:$rglb \
	 		--outSAMstrandField intronMotif \
	 		--runThreadN 16 \
	 		--chimSegmentMin 25 \
	 		--chimJunctionOverhangMin 25 \
	 		--outFileNamePrefix "$sample"_;

 		samtools index $bam

 		echo "Done aligning"
		
3. Mark duplicates in the bam and sort it with Picard 2.10.7, then index it with samtools:

	 	echo "Marking Duplicates"
	 	java -jar -XX:ParallelGCThreads=$cpu $rootdir/picard/picard.jar \
	 		MarkDuplicates I=$bam O=$dedup_bam CREATE_INDEX=true \
	 		VALIDATION_STRINGENCY=SILENT \
	 		M=$workspace/"$sample"_dedup_metrics.txt;
	
	 	java -jar -XX:ParallelGCThreads=$cpu $rootdir/picard/picard.jar \
	 		SortSam I=$dedup_bam O=$dedup_sorted_bam SORT_ORDER=coordinate;
	
	 	samtools index $dedup_sorted_bam
	 	
	 	echo "Done marking duplicates"

4. Splits reads that contain Ns in their CIGAR string using GATK 3.8's sequence data processing tools (SplitNCigarReads): 

	 	echo "SplitNCigarReads"
	 	java -jar -Xmx28g $rootdir/gatk/GenomeAnalysisTK.jar \
	 		-T SplitNCigarReads \
	 		-R $genome \
	 		-I $dedup_sorted_bam \
	 		-o $split_bam \
	 		-rf ReassignOneMappingQuality \
	 		-RMQF 255 \
	 		-RMQT 60 \
	 		-U ALLOW_N_CIGAR_READS \
	 		-fixNDN

	 	echo "Completed SplitNCigarReads"
		
5. Perform local realignment of indels and correct systematic errors in base quality scores using GATK's sequence data processing tools (RealigerTargetCreator, IndelRealigner, BaseRecalibrator, and PrintReads):

	 	echo "Performing indel realignment"
	 	java -jar -Xmx28g $rootdir/gatk/GenomeAnalysisTK.jar \
	 		-T RealignerTargetCreator \
	 		-nt $cpu \
	 		-R $genome \
	 		-I $split_bam \
	 		-known $g1000 \
	 		-known $mills \
	 		-o $workspace/"$sample".intervals;
	
	 	java -jar -Xmx28g $rootdir/gatk/GenomeAnalysisTK.jar \
	 		-T IndelRealigner \
	 		-R $genome \
	 		-I $split_bam \
	 		-targetIntervals $workspace/"$sample".intervals \
	 		-known $g1000 \
	 		-known $mills \
	 		-o $realigned_bam;

	 	echo "Done performing indel realignment"
	 	
	 	echo "\nPerforming BQSR\n"
	 	java -jar -Xmx28g $rootdir/gatk/GenomeAnalysisTK.jar \
	 		-T BaseRecalibrator \
	 		-R $genome \
	 		-I $realigned_bam \
	 		-knownSites $dbsnp \
	 		-knownSites $mills \
	 		-knownSites $g1000 \
	 		-o $workspace/"$sample"_recal_data.table \
	 		-nct $cpu;
	
	 	java -jar -Xmx28g $rootdir/gatk/GenomeAnalysisTK.jar \
	 		-T PrintReads \
	 		-R $genome \
	 		-I $realigned_bam \
	 		-BQSR $workspace/"$sample"_recal_data.table \
	 		-o $bqsr_bam
	 		
	 	echo "Done performing BQSR"	 	
				
6. Run GATK's HaplotypeCaller to identify variants in the RNA:

	 	echo "Calling Variants with HaplotypeCaller"
	 	java -jar -Xmx28g $rootdir/gatk/GenomeAnalysisTK.jar \
	 		-T HaplotypeCaller \
	 		-R $genome \
	 		-I $bqsr_bam \
	 		-dontUseSoftClippedBases \
	 		-U ALLOW_N_CIGAR_READS \
	 		-stand_call_conf 20 \
	 		--dbsnp $dbsnp \
	 		-out_mode EMIT_VARIANTS_ONLY \
	 		-rf BadCigar \
	 		-o $raw_vcf;
		
7. Filter out of the vcf file those variants with unacceptable strand bias, quality, or filtered read depth, and those that do not pass GATK's own filter:
	
		java -jar -Xmx28g $rootdir/gatk/GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R $genome \
			-V $raw_vcf \
	 		-window 35 \
	 		-cluster 3 \
	 		-filterName FS -filter "FS > 30.0" \
	 		-filterName QD -filter "QD < 2.0" \
	 		-filterName Quality -filter "QUAL < 20" \
	 		--genotypeFilterExpression "DP < 5" --genotypeFilterName DP \
	 		--setFilteredGtToNocall \
			-o $gatk_tmp_vcf;
	
		java -jar -Xmx28g $rootdir/gatk/GenomeAnalysisTK.jar \
		 -T SelectVariants \
		 -R $genome \
		 -V $gatk_tmp_vcf \
		 --excludeNonVariants \
		 --excludeFiltered \
		 -o $gatk_filt_vcf
		
		# keep only variants that passed GATK's filter
		### filter $gatk_filt_vcf for FILTER=="PASS"
		awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' \
			$gatk_filt_vcf > $gatk_filt_vcf.pass


8. Do additional variant filering using SNPiR (commit 5b571da90628bf01e6b6e2472f6efd3f917c6f41 at https://github.com/rpiskol/SNPiR ) and bedtools 2.22.1:

		echo "Performing SNPir filtering steps, "
		echo "except for RNA editing site removal";
		# convert vcf format into custom SNPiR format and 
		# filter out variants with quality <20
		$SNPiR/convertVCF.sh $gatk_filt_vcf.pass \
			$workspace/"$sample"_raw_variants.txt 20

		# filter out variants having a total depth of less than 5 reads or 
		# fewer than 2 reads with an alternate allele
		awk -F '\t|,' '{if (($3 >= 5) && ($4 >= 2)) {print}}' \
			$workspace/"$sample"_raw_variants.txt > $workspace/"$sample"_raw_variants.mincov.txt
	
		# filter out mismatches at read ends
		# note: add the -illumina option if your reads are in 
		# Illumina 1.3+ quality format
		$SNPiR/filter_mismatch_first6bp.pl \
			-infile $workspace/"$sample"_raw_variants.mincov.txt \
			-outfile $workspace/"$sample"_raw_variants.rmhex.txt \
			-bamfile $bqsr_bam
			
9. Identify variants in known Alu elements:
			
		# identify SNPiR-filtered variants in known Alu elements
		awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' \
			$workspace/"$sample"_raw_variants.rmhex.txt | \
			$intersectBed \
			-wa \
			-header \
			-a stdin \
			-b $AluRegions \
			> $workspace/"$sample"_alu_sites.txt 
			
10. For later use, capture variants whose locations are in known Alu elements as well as those whose locations both are within Alu elements and are known to be rna-editing sites:
		
		# TODO: Ask Adam precisely where this is used; can't find in dnamethod file
		# capture in VCF format ALL SNPiR-filtered variants in known Alu elements 
		# (not only those variants at known RNA-editing sites);
		# save these for later, when they will be used in identifying novel RNA editing sites
		$intersectBed \
			-wa \
			-header \
			-a $gatk_filt_vcf.pass \
			-b $workspace/"$sample"_alu_sites.txt \
			> $workspace/"$sample"_alu_rnaseq_vars.vcf 

		# filter out of the list of variants in known Alu elements those 
		# that are not at known rna-editing sites
		# Note $rnaEditDB is the radar and darned databases concatenated
		echo "Filtering for Alu RNA editing sites with RADAR and DARNED"
		$intersectBed \
			-wa \
			-a $workspace/"$sample"_alu_sites.txt \
			-b $rnaEditDB \
			| uniq \
			> $workspace/"$sample"_alu_rnaedit_sites.txt 
		
		# TODO: Ask Adam where these are used
		# capture in VCF format those SNPiR-filtered variants that are 
		# in known Alu elements AND at known rna-editing sites
		awk '{OFS="\t";print $1,$2,$2}' $workspace/"$sample"_alu_rnaedit_sites.txt | \
			$intersectBed \
			-wa \
			-header \
			-a $gatk_filt_vcf.pass \
			-b stdin \
			| uniq \
			> $workspace/"$sample"_alu_rnaedit_sites.vcf 
		
11. Identify variants NOT in known Alu elements:

		# identify SNPiR-filtered variants NOT in known Alu elements
		awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' \
			$workspace/"$sample"_raw_variants.rmhex.txt | \
			$intersectBed \
			-v \
			-a stdin \
			-b $AluRegions \
			> $workspace/"$sample"_raw_variants.rmhex.nonalu.txt

12. Further filter the variants not in known Alu elements:
						
		# TODO: Verify with Adam that the filters below here are applied 
		# ONLY to variants at non-alu sites
		# filter out variants in simple repeat regions
		awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' \
			$workspace/"$sample"_raw_variants.rmhex.nonalu.txt | \
			$intersectBed \
			-v \
			-a stdin \
			-b $RepeatMasker | \
			cut -f1,3-7 > $workspace/"$sample"_raw_variants.rmhex.nonalu.rmsk.txt
	
		# filter out variants at intronic sites that are 
		# within 4bp of splicing junctions;
		# make sure your gene annotation file is in UCSC 
		# text format and sorted by chromosome and 
		# transcript start position
		$SNPiR/filter_intron_near_splicejuncts.pl \
			-infile $workspace/"$sample"_raw_variants.rmhex.nonalu.rmsk.txt \
			-outfile $workspace/"$sample"_raw_variants.rmhex.nonalu.rmsk.rmintron.txt \
			-genefile $gene_annotation
	
		# filter out variants in homopolymers
		$SNPiR/filter_homopolymer_nucleotides.pl \
			-infile $workspace/"$sample"_raw_variants.rmhex.nonalu.rmsk.rmintron.txt \
			-outfile $workspace/"$sample"_raw_variants.rmhex.nonalu.rmsk.rmintron.rmhom.txt \
			-refgenome $genome
	
		# filter out variants caused by mismapped reads;
		# this may take a while depending on the number of variants 
		# to screen and the size of the reference genome
		# note: add the -illumina option if your reads are in 
		# Illumina 1.3+ quality format
		$SNPiR/BLAT_candidates.pl \
			-infile $workspace/"$sample"_raw_variants.rmhex.nonalu.rmsk.rmintron.rmhom.txt \
			-outfile $workspace/"$sample"_raw_variants.rmhex.nonalu.rmsk.rmintron.rmhom.rmblat.txt \
			-bamfile $bqsr_bam \
			-refgenome $genome
		
13. For later use, capture variants whose locations are not in known Alu elements as well as those whose locations not in Alu elements BUT are known to be rna-editing sites:

		# capture in VCF format ALL additionally filtered variants not in known Alu elements 
		# (not only those variants at known RNA-editing sites);
		awk '{OFS="\t";print $1,$2,$2}' \
			$workspace/"$sample"_raw_variants.rmhex.nonalu.rmsk.rmintron.rmhom.rmblat.txt | \
			$intersectBed \
			-wa \
			-header \
			-a $gatk_filt_vcf.pass \
			-b stdin \
			| uniq \
			> $workspace/"$sample"_nonalu_rnaseq_vars.vcf 

		# TODO: Ask Adam why earlier filtering for alu+rnaedit had uniq filter while this doesn't
		# filter out of the list of additionally filtered variants not in known Alu elements those 
		# that are also not at known rna-editing sites
		# Note $rnaEditDB is the radar and darned databases concatenated
		awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' \
			$workspace/"$sample"_raw_variants.rmhex.nonalu.rmsk.rmintron.rmhom.rmblat.txt | \
			$intersectBed \
			-wa \
			-a stdin \
			-b $rnaEditDB \
			> $workspace/"$sample"_nonalu_rnaedit_sites.txt 
		
		# capture in VCF format those additionally filtered variants that are 
		# not in known Alu elements but are at known rna-editing sites
		echo "Selecting RNA editing sites from VCF"
		$intersectBed \
			-wa \
			-header \
			-a $gatk_filt_vcf.pass \
			-b $workspace/"$sample"_nonalu_rnaedit_sites.txt \
			| uniq \
			> $workspace/"$sample"_nonalu_rnaedit_sites.vcf

14. In the remainder of this pipeline, combine the RNA variant analysis performed above with the DNA variant analysis described elsewhere to further stratify the list of potential editing sites (TODO: Ask Adam whether it is germline or somatic); begin by ensuring you have the necessary files:
	1. The per-chromosome bed files (TODO: Ask Adam precisely what is in these bed files) for this sample for all 24 chromosomes (1-22, X, and Y).  These should be named in the format `"$dna_sample"-PB."$chr".bed` (for example, `mysample1-PB.14.bed` for chromosome 14 of the sample `mysample1`).
	2. The gzipped vcf file containing the variants identified in the DNA for this sample.  This should be named in the format `"$dna_sample"-PB.raw.vcf.gz` (for example, `mysample1-PB.raw.vcf.gz` for the sample `mysample1`).
	
15. Index the sample's bam file:

		samtools index $bqsr_bam
		
16. Filter the per-chromosome bed files to remove any TODO:regions?? with a read depth of less than five reads, and condense the remaining TODO:regions?? into a single file:

		# Only consider sites where DNA coverage is greater than 4 reads 
		for chr in {1..22} X Y; do
			awk '{ if ($4 < 5) { print }}' \
				$workspace/"$dna_sample"-PB."$chr".bed \
				> $workspace/"$dna_sample-PB"."$chr".filt.bed
		done
		
		cat $workspace/"$dna_sample"-PB.*.filt.bed > $workspace/"$dna_sample"-PB.filt.bed

17. Identify variants whose locations both are in known Alu elements and are covered in the DNA sequencing to at least the minimum read depth: TODO: change alu_regions.txt to alu_sites.txt from here down

		# filter out of the list of variants in known Alu elements those 
		# that are not in the list of sites with adequate DNA read depth
		$intersectBed \
				-v \
				-header \
				-a $workspace/"$sample"_alu_regions.txt \
				-b $workspace/"$dna_sample-PB".filt.bed \
				> $workspace/"$sample"_alu.mincov.txt

18. For later use, capture variants whose locations both are in known Alu elements and are covered in the DNA sequencing to at least the minimum read depth, but whose genomic locations were NOT found to contain DNA variants:

		# note: concatenate the string "chr" to the beginning of each 
		# non-comment line.  Since such lines begin with the chromosome
		# identifier, this adds "chr" to the chromosome id in the 
		# DNA variant file to match the format of chromosome identifiers 
		# in the RNA variant files
		zcat $workspace/"$dna_sample"-PB.raw.vcf.gz | \
				awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | \
				$intersectBed \
				-v \
				-header \
				-a $workspace/"$sample"_alu.mincov.txt \
				-b stdin \
				> $workspace/"$sample"_rna_edit.alu.txt

		# capture these variants in VCF format
		awk '{OFS="\t";print $1,$2,$2}' \
			$workspace/"$sample"_rna_edit.alu.txt | \
			$intersectBed \
			-wa \
			-header \
			-a $gatk_filt_vcf.pass \
			-b stdin \
			> $workspace/"$sample"_rna_edit.alu.vcf

		# TODO: Doesn't it seem like these files could use better names?

19. Repeat the filtering commands in step 12 but replace the input `"$sample"_raw_variants.rmhex.nonalu.txt` with the earlier-produced input `"$sample"_raw_variants.rmhex.txt` that does NOT exclude variants in known Alu elements.  To reflect this, elide `nonalu` from the names of all output files produced in this re-run.  Thus, the output of this work will be a file of the name format `"$sample"_raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.txt`.
20. For later use, capture variants whose locations both are in known Alu elements and are covered in the DNA sequencing to at least the minimum read depth, but whose genomic locations were NOT found to contain DNA variants:

		# capture in VCF format all additionally filtered variants
		awk '{OFS="\t";print $1,$2,$2}' \
			$workspace/"$sample"_raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.txt | \
			$intersectBed \
			-wa \
			-header \
			-a $gatk_filt_vcf.pass \
			-b stdin \
			> $workspace/"$sample"_rna_edit.vcf
			
		# filter out variants whose positions are not covered 
		# in the DNA sequencing to at least the minimum read depth
		$intersectBed \
			-v \
			-header \
			-a $workspace/"$sample"_rna_edit.vcf \
			-b $workspace/"$dna_sample-PB".filt.bed \
			> $workspace/"$sample"_rna_edit.mincov.vcf
			
		# filter out variants whose genomic locations were found 
		# to contain DNA variants.
		# note: concatenate the string "chr" to the beginning of each 
		# non-comment line.  Since such lines begin with the chromosome
		# identifier, this adds "chr" to the chromosome id in the 
		# DNA variant file to match the format of chromosome identifiers 
		# in the RNA variant files 
		zcat $workspace/"$dna_sample"-PB.raw.vcf.gz | \
			awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | \
			$intersectBed \
			-v \
			-header \
			-a $workspace/"$sample"_rna_edit.mincov.vcf \
			-b stdin \
			> $workspace/"$sample"_rna_edit.final.vcf
			
		# TODO: Doesn't it seem like these files could use better names?

	
There is now a suite of vcf files for each sample.  

# TODO: Ask Adam what has replaced Annovar!

## Software References

* **star**: Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21.
* **samtools**: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9.
* **picard**: http://broadinstitute.github.io/picard/ .
* **GATK**: DePristo MA, Banks E, Poplin R, Garimella KV, Maguire JR, Hartl C, Philippakis AA, del Angel G, Rivas MA, Hanna M, McKenna A, Fennell TJ, Kernytsky AM, Sivachenko AY, Cibulskis K, Gabriel SB, Altshuler D, Daly MJ. A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nat Genet. 2011 May;43(5):491-8.
* **SNPiR**: https://github.com/rpiskol/SNPiR .
* **bedtools**: Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2.

