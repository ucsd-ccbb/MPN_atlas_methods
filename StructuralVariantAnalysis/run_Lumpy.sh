#!/bin/bash
#$ -N lumpy
#$ -cwd
#$ -pe smp 32
#$ -l h_vmem=60G

# Run LUMPY express

export PATH=/shared/workspace/software/anaconda2/bin:/shared/workspace/software/tabix/tabix-0.2.6:/shared/workspace/software/bcftools-1.7:/shared/workspace/software/samtools/samtools-1.5:/shared/workspace/software/bedtools2/bin:/shared/workspace/software/anaconda2/bin:/shared/workspace/software/vcftools/src/perl:/shared/workspace/software/anaconda2/share/lumpy-sv-0.2.14a-2/scripts:$PATH

sample=$sample
workspace=/scratch/$sample
mkdir -p $workspace
cd $workspace

aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/20171102_Frida_Jamieson_DNASeq/wgs_analysis/$sample/"$sample".final.bam $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/20171102_Frida_Jamieson_DNASeq/wgs_analysis/$sample/"$sample".final.bam.bai $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/$sample/"$sample".final.bam $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/$sample/"$sample".final.bam.bai $workspace

# Extract the discordant paired-end alignments.
samtools view -b -F 1294 $workspace/"$sample".final.bam > $workspace/"$sample".discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h $workspace/"$sample".final.bam \
    | /shared/workspace/software/anaconda2/share/lumpy-sv-0.2.14a-2/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > $workspace/"$sample".splitters.unsorted.bam

# Sort both alignments
samtools sort $workspace/"$sample".discordants.unsorted.bam -o $workspace/"$sample".discordants.bam
samtools sort $workspace/"$sample".splitters.unsorted.bam -o $workspace/"$sample".splitters.bam

lumpyexpress \
    -B $workspace/"$sample".final.bam \
    -S $workspace/"$sample".splitters.bam \
    -D $workspace/"$sample".discordants.bam \
    -o $workspace/"$sample".lumpy.vcf

sed -i 's/myelofibrosis/'"$sample"'.final/g' $workspace/"$sample".lumpy.vcf

aws s3 cp --quiet $workspace/$sample.lumpy.vcf s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/lumpy/$sample/