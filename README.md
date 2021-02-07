# AH-v-MH
A repository to hold code for the manuscript: The cis-regulatory effects of modern human-specific variants



Found here: https://www.biorxiv.org/content/10.1101/2020.10.07.330761v1




Original oligo - barcode association:

Fastqs were mapped to the original oligo library (Oligos_library_joint_noDups.fasta) using map_barcodes_comb.sh. 

Bams were processed to create lists of barcodes using process_barcode_oligo_assoc.py.

Results (oligo-barcode mappings) were converted to fasta with txt_to_fasta.py.


MPRA analysis:

RNA and DNA fastqs were concatenated into single R1, R2 and R3 files per replicate, and their lengths checked, with concatenate_and_check_fastqs.sh.

Reads were aligned to the barcode fasta using alignBowtie_RNADNA.sh.

Bams were processed to barcode counts using RNA_DNA_processing.py.

Barcode counts were formatted for MPRAnalyze using format_for_MPRAnalyze.py.

Differential activity was calculated through MPRAnalyze with MPRAnalyze_full_parallel.R.

