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

exclude_bed=/shared/workspace/software/annotation/btu356-suppl_data/btu356_LCR-hs37d5.bed/btu356_LCR-hs37d5.bed

aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/20171102_Frida_Jamieson_DNASeq/wgs_analysis/gatk/$sample/ $workspace --recursive
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/20171102_Frida_Jamieson_DNASeq/wgs_analysis/$sample/$sample.final.bam $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/20171102_Frida_Jamieson_DNASeq/wgs_analysis/$sample/$sample.final.bam.bai $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/gatk/$sample/ $workspace --recursive
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/$sample/$sample.final.bam $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/$sample/$sample.final.bam.bai $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/manta/$sample/results/variants/diploidSV.vcf.gz $workspace
aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/lumpy/$sample/$sample.lumpy.vcf.gz $workspace
# aws s3 cp s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/structual_variants/$sample/speedseq/$sample.sv.vcf.gz $workspace
# aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/melt/$sample/ALU.final_comp.vcf $workspace
# aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/melt/$sample/LINE1.final_comp.vcf $workspace
# aws s3 cp --quiet s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/melt/$sample/SVA.final_comp.vcf $workspace

rm $workspace/$sample.MT.vcf.gz

### filter $gatk_filt_vcf for FILTER=="PASS"
zcat $workspace/diploidSV.vcf.gz | awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' | bgzip > $workspace/diploidSV.filt.vcf.gz

# Concatenate SNV vcf
vcf-concat $workspace/$sample.*.vcf.gz > $workspace/$sample.snv.vcf
vcf-sort $workspace/$sample.snv.vcf | bgzip > $workspace/$sample.snv.vcf.gz
tabix -p vcf $workspace/$sample.snv.vcf.gz

# Create PED file
if zcat $workspace/$sample.snv.vcf.gz | grep myelofibrosis; then
    echo -e "0\tmyelofibrosis\t0\t0\t1" > $workspace/$sample.ped
else
    echo -e "0\t$sample\t0\t0\t1" > $workspace/$sample.ped
fi

# Run SV2
sv2 \
        -i $workspace/$sample.final.bam \
        -v $workspace/diploidSV.filt.vcf.gz $workspace/"$sample".lumpy.vcf.gz \
        -snv $workspace/$sample.snv.vcf.gz \
        -ped $workspace/$sample.ped \
        -O $workspace
sed -i 's/myelofibrosis/'"$sample"'/g' $workspace/sv2_genotypes/sv2_genotypes.vcf

# upload results to S3
aws s3 cp --quiet $workspace/$sample.lumpy.vcf.gz s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/lumpy/$sample/
aws s3 cp --quiet $workspace/sv2_genotypes/sv2_genotypes.vcf s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/SV2/$sample/
aws s3 cp --quiet $workspace/sv2_genotypes/sv2_genotypes.txt s3://lab-data-bucket/raw-sequence-data/ccbb_analysis/wgs_analysis/SV2/$sample/

rm -r $workspace