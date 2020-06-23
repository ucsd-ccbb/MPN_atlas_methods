#!/bin/bash

# Filter controls from MPN samples
workspace=/mnt/data1/adam/jamieson/holm/sv/filtered
controls_file_list=$workspace/controls/controls_file_list.txt
mpn_file_list=$workspace/mpn/mpn_file_list.txt
all_file_list=$workspace/mpn/all_file_list.txt

# Merge MPN controls 
if [ -f "$controls_file_list" ]; then
	rm $controls_file_list
fi
touch $controls_file_list
for vcf in $(ls $workspace/*PB/*PB.filtered.sv.withimprecice.vcf | grep -E "743|780|792|795"); do echo $vcf >> $controls_file_list ;  done

SURVIVOR merge \
	$controls_file_list \
	500 1 1 1 0 0 \
	$workspace/controls/MPN_controls.sv.withimprecice.vcf

# Merge MPN samples 
if [ -f "$mpn_file_list" ]; then
	rm $mpn_file_list
fi
touch $mpn_file_list
for vcf in $(ls $workspace/*PB/*PB.filtered.sv.withimprecice.vcf | grep -vE "743|780|792|795"); do echo $vcf >> $mpn_file_list ;  done
SURVIVOR merge \
	$mpn_file_list \
	500 1 1 1 0 0 \
	$workspace/mpn/MPN.sv.withimprecice.vcf

# Merge all samples
if [ -f "$all_file_list" ]; then
	rm $all_file_list
fi
touch $all_file_list
for vcf in $(ls $workspace/*PB/*PB.filtered.sv.withimprecice.vcf); do echo $vcf >> $all_file_list ;  done
SURVIVOR merge \
	$all_file_list \
	500 1 1 1 0 0 \
	$workspace/mpn/MPNandControls.sv.withimprecice.vcf


