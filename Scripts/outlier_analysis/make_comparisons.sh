#!/bin/bash

location_of_py_file="compare_groups_outliers.py"
location_of_outliers_file="../../Results/transcriptomics/OLA/ola_table.txt"
gene_column_name="Gene.name"
fdr_cut_off=0.05
blue_or_red="red"
output_qvals="True"
frac_filter="0.3"

output_prefix="../../Results/transcriptomics/OLA/5yr-survival-compare_alive"
group1_label="dead"
group1_list="../../Results/transcriptomics/OLA/G2.txt"
group2_label="alive"
group2_list="../../Results/transcriptomics/OLA/G1.txt"

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
