#!/bin/bash

bam_dir=/mnt/data1/adam/jamieson/holm/bam
genome=/mnt/data1/gxu/software/bwa_index/bwa/human_g1k_v37.fasta
regions=/mnt/data1/adam/software/annotation/UCSC.knownGene.GRCh37.exons.bed # UCSC exons
lcr=/mnt/data1/adam/software/annotation/LCR-hs37d5.bed # Heng Li low complexity regions
encode_blacklist1=/mnt/data1/adam/software/annotation/wgEncodeDacMapabilityConsensusExcludable.bed # Low mapability regions from Encode
encode_blacklist2=/mnt/data1/adam/software/annotation/wgEncodeDukeMapabilityRegionsExcludable.bed
refFlat=/mnt/data1/adam/software/annotation/refFlat_GRCh37.txt
out_dir=/mnt/data1/adam/jamieson/holm/cnvkit/pooled_normal/noLCR
normal_dir=$out_dir/normals
tumor_dir=$out_dir/tumor
mkdir -p $normal_dir
mkdir -p $tumor_dir


# Create genome accessibility bed file
cnvkit.py access $genome -o $out_dir/access.bed -x $lcr -x $encode_blacklist1 -x $encode_blacklist2
# Intersect these regions from exons bed file
intersectBed -wa -header -a $regions -b $out_dir/access.bed > $out_dir/access_exons.bed

cnvkit.py batch \
		/mnt/data1/adam/jamieson/holm/bam/tumor/*bam \
		--normal /mnt/data1/adam/jamieson/holm/bam/normal/*bam \
		--method wgs \
		--fasta $genome \
		--output-dir $out_dir \
		--annotate $refFlat \
		--access $out_dir/access_exons.bed \
		--drop-low-coverage \
		--processes 14


# The following commands are equal to cnvkit.py batch command
#cnvkit.py autobin $bam_dir/*/*.bam -t $regions -g $out_dir/access.bed --annotate $refFlat -m wgs

# For each tumor sample...
#for sample in $(cat /mnt/data1/adam/jamieson/holm/pb_samples.txt | grep -v "^#"); do
#	cnvkit.py coverage $bam_dir/tumor/$sample.final.ucscExons.bam $regions -o $tumor_dir/$sample.targetcoverage.cnn -p 14
#	cnvkit.py coverage $bam_dir/tumor/$sample.final.ucscExons.bam $regions -o $tumor_dir/$sample.antitargetcoverage.cnn -p 14
#done

# For each control sample...
#for sample in $(cat /mnt/data1/adam/jamieson/holm/pb_normal_samples.txt | grep -v "^#"); do
#        cnvkit.py coverage $bam_dir/normal/$sample.final.ucscExons.bam $regions -o $normal_dir/$sample.targetcoverage.cnn -p 14
#        cnvkit.py coverage $bam_dir/normal/$sample.final.ucscExons.bam $regions -o $normal_dir/$sample.antitargetcoverage.cnn -p 14
#done

# With all normal samples...
#cnvkit.py reference $normal_dir/*.{,anti}targetcoverage.cnn --fasta $genome --no-edge --no-gc -o $normal_dir/my_reference.cnn

# Call CNV for each tumor sample
for sample in $(cat /mnt/data1/adam/jamieson/holm/pb_samples.txt | grep -v "^#"); do
	#cnvkit.py fix $tumor_dir/$sample.targetcoverage.cnn $tumor_dir/$sample.antitargetcoverage.cnn $normal_dir/my_reference.cnn --no-edge --no-gc -o $tumor_dir/$sample.cnr
	cnvkit.py segmetrics --ci -o $tumor_dir/$sample.final.ucscExons.segmetrics -s $tumor_dir/$sample.final.ucscExons.cns $tumor_dir/$sample.final.ucscExons.cnr
	cnvkit.py call --filter ci -o $tumor_dir/$sample.final.ucscExons.segmetrics.call $tumor_dir/$sample.final.ucscExons.segmetrics
	cnvkit.py export seg -o $tumor_dir/$sample.final.ucscExons.seg $tumor_dir/$sample.final.ucscExons.segmetrics.call
	cnvkit.py export vcf -o $tumor_dir/$sample.final.ucscExons.cnvkit.vcf $tumor_dir/$sample.final.ucscExons.segmetrics.call
done
