#!/bin/bash
#$ -N cnvkit
#$ -cwd
#$ -pe smp 16
#$ -l h_vmem=30G

sample=$sample
workspace=/scratch/workspace/$sample
mkdir -p $workspace

genome=/shared/workspace/software/sequences/Hsapiens/GRCh37/human_g1k_v37.fasta
s3URL=s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/structual_variants
regions=/shared/workspace/software/annotation/Exome-AZ_V2_GRCh37.bed
refFlat=/shared/workspace/software/annotation/refFlat_GRCh37.txt
bam=$workspace/"$sample".final.bam

aws s3 cp $s3URL/$sample/$sample.final.bam $workspace/ --quiet
aws s3 cp $s3URL/$sample/$sample.final.bam.bai $workspace/ --quiet

/home/ubuntu/anaconda3/bin/cnvkit.py target $regions --split -o $workspace/$sample.target.bed --avg-size 500
/home/ubuntu/anaconda3/bin/cnvkit.py antitarget -g $regions $workspace/$sample.target.bed -o $workspace/$sample.antitarget.bed

/home/ubuntu/anaconda3/bin/cnvkit.py coverage -p 16 $workspace/"$sample".final.bam $workspace/$sample.target.bed -o $workspace/$sample.targetcoverage.cnn
/home/ubuntu/anaconda3/bin/cnvkit.py coverage -p 16 $workspace/"$sample".final.bam $workspace/$sample.antitarget.bed -o $workspace/$sample.antitarget.cnn

/home/ubuntu/anaconda3/bin/cnvkit.py reference -f $genome -o $workspace/$sample.targetcoverage-flatbackground.cnn -t $workspace/$sample.target.bed -a $workspace/$sample.antitarget.bed
/home/ubuntu/anaconda3/bin/cnvkit.py fix -o $workspace/$sample.cnvkit.cnr $workspace/$sample.targetcoverage.cnn $workspace/$sample.antitarget.cnn $workspace/$sample.targetcoverage-flatbackground.cnn
/home/ubuntu/anaconda3/bin/cnvkit.py segment -o $workspace/$sample.cnvkit.cns $workspace/$sample.cnvkit.cnr --threshold 0.00001
/home/ubuntu/anaconda3/bin/cnvkit.py segmetrics --ci -o $workspace/$sample.cnvkit.segmetrics -s $workspace/$sample.cnvkit.cns $workspace/$sample.cnvkit.cnr
/home/ubuntu/anaconda3/bin/cnvkit.py call --filter ci -o $workspace/$sample.cnvkit.segmetrics.call $workspace/$sample.cnvkit.segmetrics
/home/ubuntu/anaconda3/bin/cnvkit.py export seg -o $workspace/$sample.cnvkit.seg $workspace/$sample.cnvkit.segmetrics.call

cp $workspace/$sample.cnvkit.cnr /shared/workspace/projects/califano/cnvkit/$sample/
cp $workspace/$sample.cnvkit.cns /shared/workspace/projects/califano/cnvkit/$sample/
cp $workspace/$sample.cnvkit.segmetrics /shared/workspace/projects/califano/cnvkit/$sample/
cp $workspace/$sample.cnvkit.segmetrics.call /shared/workspace/projects/califano/cnvkit/$sample/
cp $workspace/$sample.cnvkit.seg /shared/workspace/projects/califano/cnvkit/$sample/

aws s3 cp $workspace/$sample.cnvkit.cnr $s3URL/$sample/cnvkit/ --quiet
aws s3 cp $workspace/$sample.cnvkit.cns $s3URL/$sample/cnvkit/ --quiet
aws s3 cp $workspace/$sample.cnvkit.segmetrics $s3URL/$sample/cnvkit/ --quiet
aws s3 cp $workspace/$sample.cnvkit.segmetrics.call $s3URL/$sample/cnvkit/ --quiet
aws s3 cp $workspace/$sample.cnvkit.seg $s3URL/$sample/cnvkit/ --quiet

rm -r $workspace
