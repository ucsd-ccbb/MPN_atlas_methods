#!/bin/bash

# Merge Lumpy and Manta SV callsets --> Filter SV from Controls

sample=$1
echo $sample
workspace=/mnt/data1/adam/jamieson/holm/sv


# Separate translocations (interchromosomal) from DELS/DUPS/INS/INV (intrachromosomal)
# Intra
# grep -v "SVTYPE=BND" $workspace/manta/$sample/results/variants/tumorSV.vcf | \
# 	awk -F '\t' '{if($0 ~ /#/) print; else if($7 == "PASS") print}' | grep -v ";IMPRECISE" > $workspace/manta/$sample/results/variants/tumorSV.intra.vcf
# grep -v "SVTYPE=BND|##FORMAT=<ID=BD" $workspace/lumpy/$sample/$sample.lumpy.vcf | grep -v ";IMPRECISE" > $workspace/lumpy/$sample/$sample.lumpy.intra.vcf
# # Inter
# grep -v "SVTYPE=DEL\|SVTYPE=DUP\|SVTYPE=INS" $workspace/manta/$sample/results/variants/tumorSV.vcf | \
# 	awk -F '\t' '{if($0 ~ /#/) print; else if($7 == "PASS") print}' | grep -v ";IMPRECISE" > $workspace/manta/$sample/results/variants/tumorSV.inter.vcf
# grep -v "SVTYPE=DEL\|SVTYPE=DUP\|SVTYPE=INS|##FORMAT=<ID=BD" $workspace/lumpy/$sample/$sample.lumpy.vcf | \
# 	grep -v ";IMPRECISE" > $workspace/lumpy/$sample/$sample.lumpy.inter.vcf

# filterAndMerge () {
# 	mkdir -p $workspace/filtered/$sample 
# 	type=$1
# 	echo $workspace/manta/$sample/results/variants/tumorSV."$type".vcf > $workspace/filtered/$sample/"$sample"."$type"_file_list.txt
# 	echo $workspace/lumpy/$sample/$sample.lumpy."$type".vcf >> $workspace/filtered/$sample/"$sample"."$type"_file_list.txt
# 	SURVIVOR merge \
# 		$workspace/filtered/$sample/"$sample"."$type"_file_list.txt \
# 		500 2 1 1 0 0 \
# 		$workspace/filtered/$sample/"$sample"."$type".filtered.sv.vcf
# }

# filterAndMerge intra
# filterAndMerge inter

filterAndMerge () {
	mkdir -p $workspace/filtered/$sample
	sed 's/myelofibrosis/'$sample'/g' $workspace/manta/$sample/results/variants/tumorSV.vcf | awk -F '\t' '{if($0 ~ /#/) print; else if($7 == "PASS") print}' | grep -v ";IMPRECISE" > $workspace/manta/$sample/results/variants/tumorSV.filt.vcf
	sed 's/myelofibrosis/'$sample'/g' $workspace/lumpy/$sample/$sample.lumpy.vcf | grep -v ";IMPRECISE" > $workspace/lumpy/$sample/$sample.lumpy.filt.vcf

	echo $workspace/manta/$sample/results/variants/tumorSV.filt.vcf > $workspace/filtered/$sample/"$sample"_file_list.txt
	echo $workspace/lumpy/$sample/$sample.lumpy.filt.vcf >> $workspace/filtered/$sample/"$sample"_file_list.txt
	SURVIVOR merge \
		$workspace/filtered/$sample/"$sample"_file_list.txt \
		500 2 1 1 0 0 \
		$workspace/filtered/$sample/"$sample".filtered.sv.vcf
}

filterAndMergeWithImprecice () {
	mkdir -p $workspace/filtered/$sample
	sed 's/myelofibrosis/'$sample'/g' $workspace/manta/$sample/results/variants/tumorSV.vcf | awk -F '\t' '{if($0 ~ /#/) print; else if($7 == "PASS") print}' > $workspace/manta/$sample/results/variants/tumorSV.filt.withimprecice.vcf
	sed 's/myelofibrosis/'$sample'/g' $workspace/lumpy/$sample/$sample.lumpy.vcf > $workspace/lumpy/$sample/$sample.lumpy.filt.withimprecice.vcf

	echo $workspace/manta/$sample/results/variants/tumorSV.filt.withimprecice.vcf > $workspace/filtered/$sample/"$sample"_file_list.withimprecice.txt
	echo $workspace/lumpy/$sample/$sample.lumpy.filt.withimprecice.vcf >> $workspace/filtered/$sample/"$sample"_file_list.withimprecice.txt
	SURVIVOR merge \
		$workspace/filtered/$sample/"$sample"_file_list.withimprecice.txt \
		500 2 1 1 0 0 \
		$workspace/filtered/$sample/"$sample".filtered.sv.withimprecice.vcf
}

#filterAndMerge
filterAndMergeWithImprecice


