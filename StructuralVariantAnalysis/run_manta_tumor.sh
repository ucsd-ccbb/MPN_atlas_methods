
#!/bin/bash
#$ -N manta_tumor
#$ -cwd
#$ -pe smp 32
#$ -l h_vmem=60G

# Run LUMPY express

export PATH=$PATH:/shared/workspace/software/anaconda2/bin:/shared/workspace/software/tabix/tabix-0.2.6:/shared/workspace/software/bcftools-1.7:/shared/workspace/software/samtools/samtools-1.5:/shared/workspace/software/bedtools2/bin

sample=$sample
workspace=/scratch/workspace/$sample
mkdir -p $workspace

aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/20171102_Frida_Jamieson_DNASeq/wgs_analysis/$sample/"$sample".final.bam $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/20171102_Frida_Jamieson_DNASeq/wgs_analysis/$sample/"$sample".final.bam.bai $workspace
#aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/$sample/"$sample".final.bam $workspace
#aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/$sample/"$sample".final.bam.bai $workspace

genome=/shared/workspace/software/sequences/Hsapiens/GRCh37/human_g1k_v37.fasta

runMantaTumorOnly () {
        bam=$1
        /shared/workspace/software/anaconda2/bin/configManta.py \
        --tumorBam $bam \
        --referenceFasta $genome \
        --runDir $workspace
        $workspace/runWorkflow.py -m local -j 32

	aws s3 cp $workspace/results/variants/tumorSV.vcf.gz s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/manta/$sample/results/variants/
        aws s3 cp $workspace/results/variants/tumorSV.vcf.gz.tbi s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/manta/$sample/results/variants/
}

bam=$workspace/"$sample".final.bam
runMantaTumorOnly $bam

