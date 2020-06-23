 #!/bin/bash
#$ -N sv2
#$ -cwd
#$ -pe smp 32
#$ -l h_vmem=60G

## Merge SV outputs from Manta, Melt, Speedseq with SV2

export PATH=/shared/workspace/software/tabix/tabix-0.2.6/:/shared/workspace/software/bcftools-1.7/:/shared/workspace/software/samtools/samtools-1.5/:/shared/workspace/software/bedtools2/bin/:/shared/workspace/software/anaconda2/bin:/shared/workspace/software/vcftools/src/perl:/shared/workspace/software/anaconda2/share/lumpy-sv-0.2.14a-2/scripts:$PATH
export PERL5LIB=/shared/workspace/software/vcftools/src/perl

sample=$sample
workspace=/scratch/$sample
mkdir -p $workspace
cd $workspace

exclude_bed=/shared/workspace/software/annotation/btu356_LCR-hs37d5.bed/btu356_LCR-hs37d5.bed

#aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/20171102_Frida_Jamieson_DNASeq/wgs_analysis/gatk/$sample/ $workspace --recursive
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/20171102_Frida_Jamieson_DNASeq/wgs_analysis/$sample/$sample.final.bam $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/20171102_Frida_Jamieson_DNASeq/wgs_analysis/$sample/$sample.final.bam.bai $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/$sample/$sample.raw.vcf.gz $workspace
#aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/gatk/$sample/ $workspace --recursive
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/$sample/$sample.final.bam $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/$sample/$sample.final.bam.bai $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/manta/$sample/results/variants/diploidSV.vcf.gz $workspace
# aws s3 cp s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/structual_variants/$sample/speedseq/$sample.sv.vcf.gz $workspace
# aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/melt/$sample/ALU.final_comp.vcf $workspace
# aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/melt/$sample/LINE1.final_comp.vcf $workspace
# aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/melt/$sample/SVA.final_comp.vcf $workspace

rm $workspace/$sample.MT.vcf.gz

### filter $gatk_filt_vcf for FILTER=="PASS"
zcat $workspace/diploidSV.vcf.gz | awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' | bgzip > $workspace/diploidSV.filt.vcf.gz

# Run Lumpy
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

# exclude low complexity regions
lumpyexpress \
    -B $workspace/"$sample".final.bam \
    -S $workspace/"$sample".splitters.bam \
    -D $workspace/"$sample".discordants.bam \
    -x $exclude_bed \
    -o $workspace/"$sample".lumpy.vcf

bgzip -f $workspace/"$sample".lumpy.vcf
tabix -p vcf $workspace/"$sample".lumpy.vcf.gz

# Merge melt output vcf
# vcf-concat $workspace/ALU.final_comp.vcf $workspace/LINE1.final_comp.vcf $workspace/SVA.final_comp.vcf > $workspace/melt.vcf
# vcf-sort $workspace/melt.vcf | grep -E "^#|PASS" > $workspace/melt.sorted.vcf
# sed 's/'"$sample"'.final/'"$sample"'/g' $workspace/melt.sorted.vcf | bgzip > $workspace/melt.sorted.vcf.gz
# tabix -p vcf $workspace/melt.sorted.vcf.gz

# Concatenate SNV vcf
#vcf-concat $workspace/$sample.*.vcf.gz > $workspace/$sample.snv.vcf
#vcf-sort $workspace/$sample.snv.vcf | bgzip > $workspace/$sample.snv.vcf.gz
#tabix -p vcf $workspace/$sample.snv.vcf.gz

# Create PED file
if zcat $workspace/$sample.raw.vcf.gz | grep myelofibrosis; then
    echo -e "0\tmyelofibrosis\t0\t0\t1" > $workspace/$sample.ped
else
    echo -e "0\t$sample\t0\t0\t1" > $workspace/$sample.ped
fi

# Run SV2
sv2 \
        -i $workspace/$sample.final.bam \
        -v $workspace/diploidSV.vcf.gz $workspace/"$sample".lumpy.vcf.gz \
        -snv $workspace/$sample.raw.vcf.gz \
        -ped $workspace/$sample.ped \
        -O $workspace
sed -i 's/myelofibrosis/'"$sample"'/g' $workspace/sv2_genotypes/sv2_genotypes.vcf

# upload results to S3
aws s3 cp --quiet $workspace/$sample.lumpy.vcf.gz s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/lumpy/$sample/
aws s3 cp --quiet $workspace/sv2_genotypes/sv2_genotypes.vcf s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/SV2/$sample/
aws s3 cp --quiet $workspace/sv2_genotypes/sv2_genotypes.txt s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/SV2/$sample/

rm -r $workspace
