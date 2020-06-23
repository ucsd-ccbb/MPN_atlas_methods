#!/bin/bash
#$ -N delly_call
#$ -cwd
#$ -pe smp 16
#$ -l h_vmem=60G

workspace=/scratch/workspace/$sample
mkdir -p $workspace
genome=/shared/workspace/software/sequences/Hsapiens/GRCh37/human_g1k_v37.fasta
s3URL=s3://lab-data-bucket/raw-sequence-data/20171102_Frida_Jamieson_DNASeq/wgs_analysis
s3URL=s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis
aws s3 cp $s3URL/$sample/$sample.final.bam $workspace/
aws s3 cp $s3URL/$sample/$sample.final.bai $workspace/


runDellyCall () {
	delly call \
	-o $workspace/"$sample"_delly.bcf \
	-g $genome \
	$workspace/$sample.final.bam
}

runDellyCall $sample
#rm -r $workspace
