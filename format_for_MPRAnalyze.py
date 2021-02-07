#!/usr/bin/env python
#
#SBATCH --job-name=reformat.py
#SBATCH --output=job_errors/reformat.%j.out
#SBATCH --error=job_errors/reformat.%j.err
#SBATCH --time=2:59:00
#SBATCH -p normal,owners,hns,dpetrov
#SBATCH --mem=20000
#SBATCH -c 16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cweiss19@stanford.edu
#
#
#
#
#

import os
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import itertools 
import pickle
import random

project_dir = '/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/barcode_counts/subsections' 
barcode_mapping_fasta = '/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/alignments/oligos_to_barcodes_comb_30_3.fasta'
original_oligo_list = '/scratch/users/cweiss19/MPRA/DATA/DNA_barcode_associations/Oligos_library_joint_noDups.fasta'
data_df_path = f'{project_dir}/barcode_counts_UMI_1_2.txt'
pickle_path = f'{project_dir}/MPRAnalyze_dict_nocontrols_1_2.pickle'

 


def main():    
    ### prep the data for looping ###
    
    def turn_into_locus(f):
        oligo_split = f.split("_") 
        locus = oligo_split[0] + '_' + oligo_split[1] + '_' + oligo_split[3]   # remove allele and bc_seq
        return locus
    
    def get_allele_source(f):
        oligo_split = f.split("_")    
        allele = oligo_split[2]   # either "Ancestral" or "Derived"
        return allele
    
    def get_bcseq(f):
        oligo_split = f.split("_")    
        bc_seq = oligo_split[-1]   # 15 bp bc seq is the last item always
        return bc_seq
    
    def get_oligo(f):
        oligo_split = f.split("_")
        oligo = oligo_split[0] + '_' + oligo_split[1] + '_' + oligo_split[2] + '_' + oligo_split[3]   # remove the bc_seq
        return oligo
    
    data_df = pd.read_csv(data_df_path, sep='\t', index_col = 'oligo_bc')
    
    #print(len(data_df.index), 'number of barcodes in original df')
    
    #remove ctrls from the original dataframe by creating list without them, and dropping those index rows
    indices = list(data_df.index.values)
    pos_controls = [index for index in indices if 'ctrl' in index]
    neg_controls = [index for index in indices if 'scrambled' in index]
    pos_and_neg_controls = pos_controls + neg_controls
    #print(len(pos_and_neg_controls), 'number of rows to drop')
    data_df.drop(pos_and_neg_controls, inplace=True)
    #data_df.drop(pos_controls, inplace=True)
    #print(len(data_df.index), 'number of barcodes after dropping controls')
    
    #add easily accessible info for locus, allele source and bcsequence in new columns
    data_df['locus'] = data_df.index.map(turn_into_locus)
    data_df['allele_source'] = data_df.index.map(get_allele_source)
    data_df['bc_seq'] = data_df.index.map(get_bcseq)
    data_df['oligo'] = data_df.index.map(get_oligo)
    
    #figure out the oligo with the most barcodes, per allele source
    max_bc_info_dict = {'Ancestral': [], 'Derived': []}
    df_allele_grouped = data_df.groupby('allele_source')
    for each_allele, allele_group in df_allele_grouped:
        max_barcodes = 0
        max_oligo_id = 'test'
        df_oligo_grouped = allele_group.groupby('oligo')
        for each_oligo, barcodes in df_oligo_grouped:
            number_barcodes = len(barcodes.index)
            if number_barcodes > max_barcodes:
                max_barcodes = number_barcodes
                max_oligo_id = each_oligo
        max_bc_info_dict[each_allele] = max_barcodes
        #print('For the', each_allele, 'allele', max_oligo_id, 'had the most barcodes, with', max_barcodes, 'barcodes.')
    
    # make a table using the very first original oligos, create a list of 'loci' that has the ancestral or derived categorization removed
    fasta = pysam.FastaFile(original_oligo_list)
    fasta_reference_names = fasta.references
    loci = []   # generate a list of all the loci (oligos without ancestral or derived info)
    for each_oligo in fasta_reference_names:
        if 'ctrl' in each_oligo:    # remove positive controls
            continue
        if 'scrambled' in each_oligo:      # remove negative controls
            continue
        locus = turn_into_locus(each_oligo)
        loci.append(locus)
    no_dups_loci = list(dict.fromkeys(loci))  # there will be dups because of ancestral or derived for each
    
    new_df_RNA = pd.DataFrame(index=no_dups_loci)  # initialize two dfs, with each row being a 'locus', no columns yet. base for MPRAnalyze dataframe
    new_df_DNA = pd.DataFrame(index=no_dups_loci)
    
    meta_info = data_df.loc[:, ['locus', 'allele_source', 'bc_seq', 'oligo']]
    
    # loop over each cell type + rep combo, make new df and then concat together at the end
    
    final_dict = {}
    
    for each_cell_type in ['Hob', 'ES', 'NPC']:
        cell_data = data_df.copy().filter(regex=each_cell_type)   # outputs dataframe only with the columns containing the regex string, in this case all for a given cell type (RNA and DNA, all reps)
        cell_dfs = {}
        sanity_check_1 = {}
        for each_rep in ['rep1', 'rep2', 'rep3']:
            rep_data = cell_data.copy().filter(regex=each_rep)   # outputs dataframe for each rep/celltype, includes RNA and DNA
            merged_data = rep_data.merge(meta_info, on='oligo_bc')
            
            for each_allele in ['Ancestral', 'Derived']:
                allele_df = merged_data[merged_data['allele_source'] == each_allele]   # only the rows where the allele_source col matches the looped allele
                num_cols = max_bc_info_dict[each_allele]
                identifier = each_cell_type + '_' + each_rep + '_' + each_allele
                col_names = [identifier+'_'+str(num) for num in list(range(1,num_cols+1))]   # add columns with number from 1-max col number, + identifier that includes cell type, rep and allele
                RNA_df = new_df_RNA.copy().reindex(columns=col_names)  # make new dfs for RNA and DNA, add the number of columns found earlier of max number of ancestral or derived barcodes
                DNA_df = new_df_DNA.copy().reindex(columns=col_names)
                all_loci = allele_df['locus'].tolist()  # get list of all the loci to loop over that we have data for (there will be repeats since AH and MH oligos)
                all_loci_unique = list(dict.fromkeys(all_loci))  # unique list of all the loci we have data for
                sanity_check = {}
                for each_locus in all_loci_unique:
                    locus_df = allele_df[allele_df['locus'] == each_locus]   # get the rows of the df that are from that locus
                    sanity_check[each_locus] = len(locus_df.index)
                    RNA_data = locus_df.filter(regex='RNA').iloc[:,0].tolist()  # all barcodes RNA values per locus is converted into a list    iloc[:,0] is all the rows of the first column (there's only one column from the filtering)
                    RNA_data.extend([0] * (num_cols - len(RNA_data)))  # pad the list with zeroes to the same length as the number of columns in the df
                    RNA_df.loc[each_locus,:] = RNA_data  # slot the list as a row in the new df
                    
                    DNA_data = locus_df.filter(regex='DNA').iloc[:,0].tolist()  # same for DNA
                    DNA_data.extend([0] * (num_cols - len(DNA_data)))  
                    DNA_df.loc[each_locus,:] = DNA_data
                sanity_check_1[identifier] = sanity_check
                cell_dfs[identifier] = [RNA_df, DNA_df]   # save the dfs in a dict under cell_type, rep# and allele
        
        final_dict[each_cell_type] = cell_dfs
        print(sanity_check_1)
        if each_cell_type == 'Hob':
            test1 = sanity_check_1['Hob_rep1_Derived']
            test2 = sanity_check_1['Hob_rep2_Derived']
            len1 = len(test1)
            len2 = len(test2)
            if len1 != len2:
                print('dicts not the same length')
            for num in range(1,10):
                locus = random.choice(list(test1))
                num1 = test1[locus]
                num2 = test2[locus]
                if num1 != num2:
                    print('something went wrong')
                if num1 == num2:
                    print('everything fine...for now')

    with open(pickle_path, 'wb') as handle:
        pickle.dump(final_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    return final_dict
    
    
def cell_types_separate(final_dict):    
    
    def return_rep(f):
        split = f.split("_") 
        rep = split[1]
        return rep
    def return_allele(f):
        split = f.split("_")
        allele = split[2]
        return allele
    def return_barcode(f):
        split = f.split("_")
        barcode = split[3]
        return barcode
    
    def make_annot_table(df, cell_type):
        col_names = list(df.columns.values)   # each_cell_type + '_' + each_rep + '_' + each_allele + '_' num
        annot_table = pd.DataFrame(index = col_names, columns = ['replicate', 'allele', 'barcode', 'barcode_allele'])
        annot_table['replicate'] = annot_table.index.map(return_rep)
        annot_table['allele'] = annot_table.index.map(return_allele)
        annot_table['barcode'] = annot_table.index.map(return_barcode)
        annot_table['barcode_allele'] = annot_table.allele.str.cat(annot_table.barcode, sep='.')
        annot_table.to_csv(f'{project_dir}/{cell_type}_col_annot_1_2.txt', sep='\t', header=True)
    
    for each_cell_type, dfs in final_dict.items():  # iterate over each cell type and fuse all the dfs together
        to_fuse_RNA = []
        to_fuse_DNA = []
        for each_id, list_of_dfs in dfs.items():
            to_fuse_RNA.append(list_of_dfs[0])
            to_fuse_DNA.append(list_of_dfs[1])
        concat_df_RNA = pd.concat(to_fuse_RNA, axis=1)
        concat_df_DNA = pd.concat(to_fuse_DNA, axis=1)
        concat_df_RNA.fillna(0, inplace=True)
        concat_df_DNA.fillna(0, inplace=True)
        
        final_df_RNA = concat_df_RNA
        final_df_DNA = concat_df_DNA
        
        shape_RNA = final_df_RNA.shape
        shape_DNA = final_df_DNA.shape
        
        if not shape_RNA == shape_DNA:   # the RNA and DNA matrices need to be the exact same size (same number of columns and rows), check here
            print('something went wrong and you need to fix that bug')
        
        make_annot_table(final_df_RNA, each_cell_type)   # the same for both, only needs to be done once, just feed RNA in (could have been the DNA table too)
        
        final_df_RNA.to_csv(f'{project_dir}/{each_cell_type}_RNA_MPRAnalyze_1_2.txt', sep='\t', header=True, index_label='locus')
        final_df_DNA.to_csv(f'{project_dir}/{each_cell_type}_DNA_MPRAnalyze_1_2.txt', sep='\t', header=True, index_label='locus')

def reps_divided(final_dict):
    
    def return_rep(f):
        split = f.split("_") 
        rep = split[1]
        return rep
    def return_allele(f):
        split = f.split("_")
        allele = split[2]
        return allele
    def return_barcode(f):
        split = f.split("_")
        barcode = split[3]
        return barcode
    
    def make_annot_table(df, cell_type, rep_num):
        col_names = list(df.columns.values)   # each_cell_type + '_' + each_rep + '_' + each_allele + '_' num
        annot_table = pd.DataFrame(index = col_names, columns = ['replicate', 'allele', 'barcode', 'barcode_allele'])
        annot_table['replicate'] = annot_table.index.map(return_rep)
        annot_table['allele'] = annot_table.index.map(return_allele)
        annot_table['barcode'] = annot_table.index.map(return_barcode)
        annot_table['barcode_allele'] = annot_table.allele.str.cat(annot_table.barcode, sep='.')
        annot_table.to_csv(f'{project_dir}/{cell_type}_{rep_num}_col_annot.txt', sep='\t', header=True)
    
    for each_cell_type, dfs in final_dict.items():  # iterate over each cell type and fuse all the dfs together
        
        rep_dict = {'rep1': [[],[]], 'rep2': [[],[]], 'rep3': [[],[]]}
        
        for each_id, list_of_dfs in dfs.items():
            split = each_id.split("_")
            rep = split[1]
            rep_dict[rep][0].append(list_of_dfs[0])
            rep_dict[rep][1].append(list_of_dfs[1])
        
        rep_finals = {}
        
        for each_rep, list_of_lists_of_dfs in rep_dict.items():
            RNA_dfs = list_of_lists_of_dfs[0]
            #print(len(RNA_dfs), 'should be 2')

            DNA_dfs = list_of_lists_of_dfs[1]
            concat_df_RNA = pd.concat(RNA_dfs, axis=1)
            concat_df_DNA = pd.concat(DNA_dfs, axis=1)
            concat_df_RNA.fillna(0, inplace=True)
            concat_df_DNA.fillna(0, inplace=True)
            
            final_df_RNA = concat_df_RNA
            final_df_DNA = concat_df_DNA
            
            shape_RNA = final_df_RNA.shape
            shape_DNA = final_df_DNA.shape
            
            if not shape_RNA == shape_DNA:   # the RNA and DNA matrices need to be the exact same size (same number of columns and rows), check here
                print('something went wrong and you need to fix that bug')
            
            make_annot_table(final_df_RNA, each_cell_type, each_rep)   # the same for both, only needs to be done once, just feed RNA in (could have been the DNA table too)
            
            final_df_RNA.to_csv(f'{project_dir}/{each_cell_type}_{each_rep}_RNA_MPRAnalyze_1_2.txt', sep='\t', header=True, index_label='locus')
            final_df_DNA.to_csv(f'{project_dir}/{each_cell_type}_{each_rep}_DNA_MPRAnalyze.txt_1_2', sep='\t', header=True, index_label='locus')



if os.path.exists(pickle_path):
    final_dict = pickle.load(open(pickle_path, 'rb'))
    print('found the pickle')
else:
    final_dict =  main()
print('onto the second step')
# run one of these  
cell_types_separate(final_dict)
#reps_divided(final_dict)




    
    

