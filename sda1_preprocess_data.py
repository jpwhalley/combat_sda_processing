def compare():
    """Gets the data in a format for Victoria Hoare's tensor and matrix decomposition, as well as the command on the rescomp to decompose the tensors. 
        INPUT: The counts for bulk-RNA-seq, the scRNA-seq counts classified into major lineage pseudo-bulk-RNA-seq is needed. As well as the cell counts from the CITEseq, CyTOF all cells, and CYTOF depleted cells and the normalised data from Luminex and timsTOF for proteins. Finally the clinical data to match up the different IDs used for each of these datasets.
        FUNCTION: compare()
        OUTPUT: The 3D tensor of cell type, donors and gene expression, and the 2D matrices of the other modalities decomposed into a series of components. All normalised (if not already normalised before - the protein data). Including missing data. All formatted for input into SDA algorithm (space delimited, NA for missing data) and a batch file to submit this to the BMRC cluster here in Oxford.
        CONTACT: Justin Whalley; currently jpw@well.ox.ac.uk ; 2021/06/28
        Time taken: When submitted to the cluster, each of the 10 iterations took around 1 day 7 hours 40 minutes, running on 24 nodes (256GB of RAM) of Intel Haswell​ CPUs (2.5 GHz). 
    """
    
    import numpy as np
    import scipy as sp
    import pandas as pd
    from natsort import natsorted, index_natsorted, order_by_index
    from os import path, makedirs
    
    # Get the clinical data in and lets restrict to COVID-19
    hosp_data = pd.read_csv('data/Final_COMBAT_basic_clinical_data_freeze_170820.txt',index_col=0,sep='\t')
    print('Start: ', len(hosp_data), '\n\n')
    
    # Expression data
    bulkRNA = pd.read_csv('data/Counts_143_23063.txt',index_col=0,sep='\t')
    bulkRNA.rename(columns=lambda x: x.replace('.','-'), inplace=True)
    scRNA = pd.read_csv('data/pseudobulk_counts_lin.txt',index_col=0,sep='\t')
    # scRNA = pd.read_csv('data/pseudobulk_counts_narrow.txt',index_col=0,sep='\t')
    
    # Cell count data
    citeseq = pd.read_csv('data/comp_total_name.unique_counts.csv', index_col=0)
    cytof_non = pd.read_csv('data/COVID_CyTOFnondepleted_CellCounts.csv', index_col=0).T
    cytof_dep = pd.read_csv('data/COVID_CyTOFdepleted_Cellcounts_unique_total.csv', index_col=0).T
    
    # Proteomic Data
    masspec = pd.read_csv('data/ints.tsv',index_col=0,sep='\t')
    luminex = pd.read_csv('data/luminex_data_processed_filtered.txt',index_col=0, sep='\t')
    
    pseudobulk_ids = []
    for col in scRNA.columns.values:
        if '_B' in col:
            temp = (col.replace('_B', '').replace('.', '-'))
            if temp not in ['U00601-Ua005E-PBUa', 'U00613-Ua005E-PBUa', 'N00049-Ja005E-PBGa', 'G05092-Ja005E-PBCa', 'S00030-Ja003E-PBCa', 'U00504-Ua005E-PBUa', 'U00515-Ua005E-PBUa', 'U00519-Ua005E-PBUa']: # Excluded due to mismatch in clinical data (x2), too few cells in citeSEQ (x3), patients assayed in too few modalities (x2) 
                pseudobulk_ids.append(temp)
    
    common_samps = {}
    for ind in hosp_data.index.values:
        if hosp_data.RNASeq_sample_ID.loc[ind] in bulkRNA.columns.values:
            if ind not in common_samps:
                common_samps[ind] = ['RNA']
            else:
                common_samps[ind].append('RNA')
        
        if hosp_data.scRNASeq_sample_ID.loc[ind] in pseudobulk_ids:
            if ind not in common_samps:
                common_samps[ind] = ['RNA']
            else:
                common_samps[ind].append('RNA')
                
        if hosp_data.scRNASeq_sample_ID.loc[ind] in citeseq.columns.values:
            if ind not in common_samps:
                common_samps[ind] = ['Cells']
            else:
                common_samps[ind].append('Cells')
                
        if hosp_data.COMBAT_participant_timepoint_ID.loc[ind] in cytof_dep.columns.values:
            if ind not in common_samps:
                common_samps[ind] = ['Cells']
            else:
                common_samps[ind].append('Cells')
                
        if hosp_data.COMBAT_proteomics_sample_ID.loc[ind] in masspec.columns.values:
            if ind not in common_samps:
                common_samps[ind] = ['Proteins']
            else:
                common_samps[ind].append('Proteins')
                
        if hosp_data.COMBAT_Luminex_sample_ID.loc[ind] in luminex.index.values:
            if ind not in common_samps:
                common_samps[ind] = ['Proteins']
            else:
                common_samps[ind].append('Proteins')
    print(len(common_samps))
    
    # Delete samples which we do not have RNA-seq data for (either bulk or pseudobulk)
    to_del = []
    samps_in_common = []
    for ind in common_samps:
        if 'RNA' not in common_samps[ind]:
            to_del.append(ind)
        else:
            samps_in_common.append(ind)
    
    for key in to_del:
        del(common_samps[key])
        
    print(len(common_samps))

    hosp_data = hosp_data.loc[samps_in_common]
    
    # For scRNA in pseudobulk format
    scRNA_samps = []
    for samp in list(scRNA.columns.values):
        temp = samp.split('_')
        if temp[0].replace('.', '-') not in scRNA_samps:
            scRNA_samps.append(temp[0].replace('.', '-'))
    
    genes_in_common = set(bulkRNA.index.values).intersection(scRNA.index.values)
    
    scRNA = scRNA.loc[genes_in_common]
    print('Any NaNs in scRNA: ', scRNA.isnull().sum().sum())
    pseudo_samps = []
    pseudo_clust = []
    for col in scRNA:
        # temp = col.split('-')
        semp = col.split('_')
        if semp[1] not in pseudo_clust:
            pseudo_clust.append(semp[1])
        if semp[0].replace('.', '-') not in pseudo_samps:
            pseudo_samps.append(semp[0].replace('.', '-'))
    
    for ident in hosp_data.index.values:
        if hosp_data.scRNASeq_sample_ID.loc[ident] == hosp_data.scRNASeq_sample_ID.loc[ident]:
            for clust in pseudo_clust:
                if hosp_data.scRNASeq_sample_ID.loc[ident].replace('-', '.') + '_' + clust not in scRNA:
                    scRNA[hosp_data.scRNASeq_sample_ID.loc[ident].replace('-', '.') + '_' + clust] = np.nan
                elif hosp_data.scRNASeq_sample_ID.loc[ident] not in citeseq.columns.values:
                    for clust in pseudo_clust:
                        scRNA[hosp_data.scRNASeq_sample_ID.loc[ident].replace('-', '.') + '_' + clust] = np.nan
        else:
            hosp_data.scRNASeq_sample_ID.loc[ident] = hosp_data.COMBAT_participant_timepoint_ID.loc[ident]
            for clust in pseudo_clust:
                scRNA[hosp_data.scRNASeq_sample_ID.loc[ident].replace('-', '.') + '_' + clust] = np.nan
            
    pseudobulk = {}
    for clust in pseudo_clust:
        pseudo_samps = []
        for samp in hosp_data.scRNASeq_sample_ID.values:
            pseudo_samps.append(samp.replace('-', '.') + '_' + clust)
        pseudobulk[clust] = scRNA[pseudo_samps]
        pseudobulk[clust][pseudobulk[clust].columns.values] = (pseudobulk[clust][pseudobulk[clust].columns.values] / (pseudobulk[clust][pseudobulk[clust].columns.values].sum()+1)) * 1000000
        pseudobulk[clust] = np.log2(pseudobulk[clust]+1)
    
        print('We have scRNA '+ clust + ': ', pseudobulk[clust].iloc[0].notnull().sum())
        print('We have not scRNA'+ clust + ': ', pseudobulk[clust].iloc[0].isnull().sum())
    print('\n')
    
    # citeseq
    we_have = 0
    we_have_not = 0
    cell_counts = pd.DataFrame(np.nan, index=citeseq.index.values, columns=hosp_data.scRNASeq_sample_ID.values)
    for samp in hosp_data.scRNASeq_sample_ID.values: 
        if samp in citeseq.columns.values:
            we_have += 1
            if citeseq[samp].sum() > 500:
                cell_counts[samp] = citeseq[samp]
            else:
                print('Citeseq ', samp, ': has ', citeseq[samp].sum(), ' cells\n')
        else:
            we_have_not += 1
    cell_counts[cell_counts.columns.values] = (cell_counts[cell_counts.columns.values] / (cell_counts[cell_counts.columns.values].sum()+1)) * 1000000
    cell_counts = np.log2(cell_counts+1)
    
    print('citeseq; we have ', we_have, 'samples')
    print('citeseq; we have not', we_have_not, 'samples\n\n')
    
    # bulkRNA - Want to keep Flu samples for this
    bulkRNA = bulkRNA.loc[genes_in_common]
    print('Any NaNs in whole blood: ', bulkRNA.isnull().sum().sum())
    samps_to_use = []
    we_have = 0
    we_have_not = 0
    for ident in hosp_data.index.values:
        if hosp_data.RNASeq_sample_ID.loc[ident] == hosp_data.RNASeq_sample_ID.loc[ident]:
            samps_to_use.append(hosp_data.RNASeq_sample_ID.loc[ident])
            we_have += 1
        else:
            samps_to_use.append(hosp_data.COMBAT_participant_timepoint_ID.loc[ident])
            we_have_not += 1
    print('Whole blood; we have ', we_have, 'samples')
    print('Whole blood; we have not', we_have_not, 'samples\n\n')
    
    whole_blood = pd.DataFrame(np.nan, index=genes_in_common, columns=samps_to_use)
    for col in whole_blood:
        if col in bulkRNA:
            if bulkRNA[col].isnull().sum() > 0:
                print('Whole blood ', col, ': has ', bulkRNA[col].isnull().sum(), ' NaNs')
            elif bulkRNA[col].sum() == 0:
                print('Whole blood ', col, ': is zero')
            
            whole_blood[col] = bulkRNA[col]
    
    whole_blood[whole_blood.columns.values] = (whole_blood[whole_blood.columns.values] / (whole_blood[whole_blood.columns.values].sum()+1)) * 1000000
    whole_blood = np.log2(whole_blood+1)
    
    # luminex - had to change the ids from Ja to Ha or vice versa in the raw data - already normalised and log scaled
    luminex = luminex.T
    print('Any NaNs in Luminex: ', luminex.isnull().sum().sum())
    samps_to_use = []
    we_have = 0
    we_have_not = 0
    for ident in hosp_data.index.values:
        if hosp_data.COMBAT_Luminex_sample_ID.loc[ident] == hosp_data.COMBAT_Luminex_sample_ID.loc[ident]:
            samps_to_use.append(hosp_data.COMBAT_Luminex_sample_ID.loc[ident])
            we_have += 1
        else:
            samps_to_use.append(hosp_data.COMBAT_participant_timepoint_ID.loc[ident])
            we_have_not += 1
    print('Luminex; we have ', we_have, 'samples')
    print('Luminex; we have not', we_have_not, 'samples\n\n')
    
    luminex_df = pd.DataFrame(np.nan, index=luminex.index.values, columns=samps_to_use)
    for col in luminex_df:
        if col in luminex:
            if luminex[col].isnull().sum() > 0:
                print('Luminex ', col, ': has ', luminex[col].isnull().sum(), ' NaNs')
            elif luminex[col].sum() == 0:
                print('Luminex ', col, ': is zero')
            
            luminex_df[col] = luminex[col]
    
    # masspec - already normalised and log scaled
    print('Any NaNs in masspec: ', masspec.isnull().sum().sum())
    samps_to_use = []
    we_have = 0
    we_have_not = 0
    for ident in hosp_data.index.values:
        if hosp_data.COMBAT_proteomics_sample_ID.loc[ident] == hosp_data.COMBAT_proteomics_sample_ID.loc[ident]:
            samps_to_use.append(hosp_data.COMBAT_proteomics_sample_ID.loc[ident])
            we_have += 1
        else:
            samps_to_use.append(hosp_data.COMBAT_participant_timepoint_ID.loc[ident])
            we_have_not += 1
    print('masspec; we have ', we_have, 'samples')
    print('masspec; we have not', we_have_not, 'samples\n\n')
    
    masspec_df = pd.DataFrame(np.nan, index=masspec.index.values, columns=samps_to_use)
    for col in masspec_df:
        if col in masspec:
            if masspec[col].isnull().sum() > 0:
                print('masspec ', col, ': has ', masspec[col].isnull().sum(), ' NaNs')
            elif masspec[col].sum() == 0:
                print('masspec ', col, ': is zero')
            
            masspec_df[col] = masspec[col]
    
    # cytof non
    print('Any NaNs in cytof_non: ', cytof_non.isnull().sum().sum())
    samps_to_use = []
    we_have = 0
    we_have_not = 0
    for ident in hosp_data.index.values:
        if hosp_data.COMBAT_participant_timepoint_ID.loc[ident] in cytof_non.columns.values:
            samps_to_use.append(hosp_data.COMBAT_participant_timepoint_ID.loc[ident])
            we_have += 1
        else:
            samps_to_use.append(hosp_data.COMBAT_participant_timepoint_ID.loc[ident])
            we_have_not += 1
    print('cytof_non; we have ', we_have, 'samples')
    print('cytof_non; we have not', we_have_not, 'samples\n\n')
    
    cytof_non_df = pd.DataFrame(np.nan, index=cytof_non.index.values, columns=samps_to_use)
    for col in cytof_non_df:
        if col in cytof_non:
            if cytof_non[col].sum() > 500:
                cytof_non_df[col] = cytof_non[col]
            else:
                print('cytof_non ', samp, ': has ', cytof_non[col].sum(), ' cells\n')
            
            
            
    cytof_non_df[cytof_non_df.columns.values] = (cytof_non_df[cytof_non_df.columns.values] / (cytof_non_df[cytof_non_df.columns.values].sum()+1)) * 1000000
    cytof_non_df = np.log2(cytof_non_df+1)
    
    # cytof depleted
    print('Any NaNs in cytof_dep: ', cytof_dep.isnull().sum().sum())
    samps_to_use = []
    we_have = 0
    we_have_not = 0
    for ident in hosp_data.index.values:
        if hosp_data.COMBAT_participant_timepoint_ID.loc[ident] in cytof_dep.columns.values:
            samps_to_use.append(hosp_data.COMBAT_participant_timepoint_ID.loc[ident])
            we_have += 1
        else:
            samps_to_use.append(hosp_data.COMBAT_participant_timepoint_ID.loc[ident])
            we_have_not += 1
    print('cytof_dep; we have ', we_have, 'samples')
    print('cytof_dep; we have not', we_have_not, 'samples\n\n')
    
    cytof_dep_df = pd.DataFrame(np.nan, index=cytof_dep.index.values, columns=samps_to_use)
    for col in cytof_dep_df:
        if col in cytof_dep:
            if cytof_dep[col].sum() > 500:
                cytof_dep_df[col] = cytof_dep[col]
            else:
                print('cytof_dep ', samp, ': has ', cytof_dep[col].sum(), ' cells\n')
            
    cytof_dep_df[cytof_dep_df.columns.values] = (cytof_dep_df[cytof_dep_df.columns.values] / (cytof_dep_df[cytof_dep_df.columns.values].sum()+1)) * 1000000
    cytof_dep_df = np.log2(cytof_dep_df+1)
    
    # Final steps, combine RNA, set up the clinical data and go
    to_concat = [whole_blood]
    for clust in pseudo_clust:
        to_concat.append(pseudobulk[clust])
    rna = pd.concat(to_concat, axis=1)
    
    hosp_data = hosp_data.set_index('COMBAT_participant_timepoint_ID')
    
    
    # Wrtie it out and see what happens
    directory = 'c_new_pseudo/'

    if not path.exists(directory):
        makedirs(directory)
    
    rna.T.to_csv('c_new_pseudo/rna.txt', header=False, index=False, sep=' ', na_rep='NA')
    cell_counts.T.to_csv('c_new_pseudo/citeseq.txt', header=False, index=False, sep=' ', na_rep='NA')
    cytof_non_df.T.to_csv('c_new_pseudo/cytof_non.txt', header=False, index=False, sep=' ', na_rep='NA')
    cytof_dep_df.T.to_csv('c_new_pseudo/cytof_dep.txt', header=False, index=False, sep=' ', na_rep='NA')
    luminex_df.T.to_csv('c_new_pseudo/luminex.txt', header=False, index=False, sep=' ', na_rep='NA')
    masspec_df.T.to_csv('c_new_pseudo/masspec.txt', header=False, index=False, sep=' ', na_rep='NA')
    
    f = open('c_new_pseudo/md_samples.txt', 'w')
    for samp in hosp_data.index.values:
        f.write(samp + '\n')
    f.close()
    
    f = open('c_new_pseudo/md_tissues.txt', 'w')
    f.write('whole_blood' + '\n')
    for clust in pseudo_clust:
        f.write(clust + '\n')
    f.close()
    
    f = open('c_new_pseudo/md_genes.txt', 'w')
    for samp in rna.index.values:
        f.write(samp + '\n')
    f.close()
    
    f = open('c_new_pseudo/md_citeseq_cells.txt', 'w')
    for samp in cell_counts.index.values:
        f.write(samp + '\n')
    f.close()
    
    f = open('c_new_pseudo/md_cytof_non_cells.txt', 'w')
    for samp in cytof_non_df.index.values:
        f.write(samp + '\n')
    f.close()
    
    f = open('c_new_pseudo/md_cytof_dep_cells.txt', 'w')
    for samp in cytof_dep_df.index.values:
        f.write(samp + '\n')
    f.close()
    
    f = open('c_new_pseudo/md_luminex_proteins.txt', 'w')
    for samp in luminex_df.index.values:
        f.write(samp + '\n')
    f.close()
    
    f = open('c_new_pseudo/md_masspec_proteins.txt', 'w')
    for samp in masspec_df.index.values:
        f.write(samp + '\n')
    f.close()
    
    hosp_data.to_csv('c_new_pseudo/clinical_data_for_miss.csv')
    
    
    # Sort the cluster commands out
    f = open('c_new_pseudo/matd.sh', "w")

    f.write('#!/bin/bash\n\n')

    f.write('#$ -cwd\n')
    f.write('#$ -N cov -j y\n')
    f.write('#$ -q covid.q -P combat.prjc\n')
    f.write('#$ -t 1-10:1\n')
    f.write('#$ -pe shmem 24\n')
    f.write('#$ -r y\n\n')

    f.write('echo "#######################################################################"\n')
    f.write('echo "SGE job id: "$JOB_ID\n')
    f.write('echo "SGE task id: "$SGE_TASK_ID\n')
    f.write('echo "Run on host: "`hostname`\n')
    f.write('echo "Operating system: "`uname -s`\n')
    f.write('echo "Username: "`whoami`\n')
    f.write('echo "Started at: "`date`\n')
    f.write('echo "#######################################################################"\n\n')

    f.write('mkdir /well/combat/users/bru597/td/sda/c_new_pseudo/results${SGE_TASK_ID}\n')
    f.write('export OPENMP_NUM_THREADS=${NSLOTS:-1}\n')
    
    to_write = '/apps/well/sda/1.1/./sda --data '
    
    to_write = to_write + '/well/combat/users/bru597/td/sda/c_new_pseudo/rna.txt '
    to_write = to_write + '/well/combat/users/bru597/td/sda/c_new_pseudo/citeseq.txt '
    to_write = to_write + '/well/combat/users/bru597/td/sda/c_new_pseudo/cytof_non.txt '
    to_write = to_write + '/well/combat/users/bru597/td/sda/c_new_pseudo/cytof_dep.txt '
    to_write = to_write + '/well/combat/users/bru597/td/sda/c_new_pseudo/luminex.txt '
    to_write = to_write + '/well/combat/users/bru597/td/sda/c_new_pseudo/masspec.txt '
    
    to_write = to_write + '--N '+str(len(hosp_data.index.values))+' --out /well/combat/users/bru597/td/sda/c_new_pseudo/results${SGE_TASK_ID} --num_openmp_threads ${NSLOTS:-1} --num_blocks ${NSLOTS:-1} --eigen_parallel true --num_comps 1000 --save_freq 500 --remove_zero_comps true --max_iter 3000 --ignore_missing true \n\n'
    
    f.write(to_write)

    f.write('echo "#########################################################################"\n')
    f.write('echo "Finished at: "`date`\n')
    f.write('echo "#########################################################################"\n')
    f.write('exit 0')

    f.close()
    