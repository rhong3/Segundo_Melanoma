#!/bin/bash

location_of_py_file="compare_groups_outliers.py"
location_of_outliers_file="/Users/rh2740/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM600/ola_table.txt"
gene_column_name="Accession"
fdr_cut_off=0.05
blue_or_red="red"
output_qvals="True"
frac_filter="0.3"

output_prefix="/Users/rh2740/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM600/more_"
group2_label="less"
group2_list="/Users/rh2740/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM600/G1.txt"
group1_label="more"
group1_list="/Users/rh2740/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM600/G2.txt"

python3  ${location_of_py_file} \
--outliers_table  ${location_of_outliers_file} \
--gene_column_name ${gene_column_name} \
--fdr_cut_off ${fdr_cut_off} \
--output_prefix ${output_prefix} \
--group1_label ${group1_label} \
--group1_list ${group1_list} \
--group2_label ${group2_label} \
--group2_list ${group2_list} \
--blue_or_red ${blue_or_red} \
--output_qvals ${output_qvals} \
--frac_filter ${frac_filter}
