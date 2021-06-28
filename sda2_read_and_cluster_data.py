def cluster():
    """Takes the decomposed tensors from sda and clusters them together
		INPUT: Names for treatment, genes, samples as well as the decomposed tensor
		FUNCTION: cluster()
		OUTPUT: The final loading_scores and posterior inclusion probabilities.
        CONTACT: Justin Whalley; currently jpw@well.ox.ac.uk ; 2021/06/28
        Time taken: 25min 9s"""
    
    import numpy as np
    import scipy as sp
    import pandas as pd
    import scipy.cluster.hierarchy as hac
    from scipy import stats
    from collections import Counter
    from natsort import natsorted, index_natsorted, order_by_index
    from os import path, makedirs
    
    def my_metric(x, y):
        r = abs(stats.pearsonr(x, y)[0])
        return 1 - r # absolute correlation to distance: range 0 to 1
        
    # Get the clinical data in and lets restrict to COVID-19
    samps_in_common = open('md_samples.txt').read().splitlines()
    
    tissues = open('md_tissues.txt').read().splitlines()
    
    rna_genes = open('md_genes.txt').read().splitlines()
    
    cell_counts = open('md_citeseq_cells.txt').read().splitlines()
    
    cytof_non_df = open('md_cytof_non_cells.txt').read().splitlines()
    
    cytof_dep_df = open('md_cytof_dep_cells.txt').read().splitlines()
    
    luminex_proteins = open('md_luminex_proteins.txt').read().splitlines()
    
    masspec_proteins = open('md_masspec_proteins.txt').read().splitlines()
    
    
    col_names = {1:rna_genes, 2:cell_counts, 3:cytof_non_df, 4:cytof_dep_df, 5:luminex_proteins, 6:masspec_proteins}
            
    
    # Now do all
    indices = []
    start = True
    directory = 'c_new_pseudo/results'
    no_of_iterations = range(1,11)
    
    for i in no_of_iterations:
        if start:
            # Get the names for samples
            samples = pd.read_csv(directory+str(i)+'/it3000/A', sep=' ', header=None)
            samples = samples.T
            samples.columns = samps_in_common
            sind_to_del = samples[(samples.T == 0).all()].index.values
            sind_to_keep = samples[(samples.T != 0).any()].index.values
            samples = samples[(samples.T != 0).any()]
            index = list(range(len(samples))) # Set up the correct index for all
            samples.index = index
            indices.append(index)
            
            # Matrices
            tind_to_del = {}
            genes = {}
            gind_to_del = {}
            pip = {}
            for j in range(1,7):
                # genes
                genes[j] = pd.read_csv(directory+str(i)+'/it3000/X' + str(j), sep=' ', header=None)
                genes[j].columns = col_names[j]
                # print(j, len(genes[j]))
                # genes[j] = genes[j].reindex(corr_info['ProbeID'], axis=1)
                gind_to_del[j] = genes[j][(genes[j].T == 0).all()].index.values
                # ind_to_keep = (genes[j].T != 0).any()
                genes[j] = genes[j].loc[sind_to_keep]
                genes[j].index = index
                # pip
                pip[j] = pd.read_csv(directory+str(i)+'/it3000/S' + str(j), sep=' ', header=None)
                pip[j].columns = col_names[j]
                # expr_info = expr_info.reindex(index=order_by_index(expr_info.index, index_natsorted(zip(expr_info.chrom, expr_info.chromStart))))
                # pip[j] = pip[j].reindex(expr_info['ProbeID'], axis=1)
                pip[j] = pip[j].loc[sind_to_keep]
                pip[j].index = index
                
                if j == 1:
                    treatment = pd.read_csv(directory+str(i)+'/it3000/B' + str(j), sep=' ', header=None)
                    treatment = treatment.T
                    treatment.columns = tissues
                    tind_to_del = treatment[(treatment.T == 0).all()].index.values
                    treatment = treatment[(treatment.T != 0).any()]
                    treatment.index = index
            
            print(i)
            for j in range(1,7):
                if j == 1:
                    print('Results ' + str(j) + ': ' +str([len(sind_to_del),len(gind_to_del[j]), len(tind_to_del)]) + ' zero rows dropped')
                else:
                    print('Results ' + str(j) + ': ' +str([len(sind_to_del),len(gind_to_del[j])]) + ' zero rows dropped')
            print('-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-')
            # To keep on adding to these dataframes
            start = False
        else:
            # samples
            samp_temp = pd.read_csv(directory+str(i)+'/it3000/A', sep=' ', header=None)
            samp_temp = samp_temp.T
            samp_temp.columns = samps_in_common
            sind_to_del = samp_temp[(samp_temp.T == 0).all()].index.values
            sind_to_keep = samp_temp[(samp_temp.T != 0).any()].index.values
            samp_temp = samp_temp[(samp_temp.T != 0).any()]
            index = list(range(max(samples.index),max(samples.index)+len(samp_temp)))
            samp_temp.index = index
            samples = samples.append(samp_temp)
            indices.append(index)
            
            # Matrices
            for j in range(1,7):
                genes_temp = pd.read_csv(directory+str(i)+'/it3000/X' + str(j), sep=' ', header=None)
                genes_temp.columns = col_names[j]
                # genes_temp = genes_temp.reindex(corr_info['ProbeID'], axis=1)
                gind_to_del[j] = genes_temp[(genes_temp.T == 0).all()].index.values
                # ind_to_keep = (genes_temp.T != 0).any()
                genes_temp = genes_temp.loc[sind_to_keep]
                genes_temp.index = index
                genes[j] = genes[j].append(genes_temp)
                # pip
                pip_temp = pd.read_csv(directory+str(i)+'/it3000/S' + str(j), sep=' ', header=None)
                pip_temp.columns = col_names[j]
                # expr_info = expr_info.reindex(index=order_by_index(expr_info.index, index_natsorted(zip(expr_info.chrom, expr_info.chromStart))))
                # pip_temp = pip_temp.reindex(expr_info['ProbeID'], axis=1)
                pip_temp = pip_temp.loc[sind_to_keep]
                pip_temp.index = index
                pip[j] = pip[j].append(pip_temp)
                
                if j == 1:
                    treat_temp = pd.read_csv(directory+str(i)+'/it3000/B' + str(j), sep=' ', header=None)
                    treat_temp = treat_temp.T
                    treat_temp.columns = tissues
                    tind_to_del = treat_temp[(treat_temp.T == 0).all()].index.values
                    treat_temp = treat_temp[(treat_temp.T != 0).any()]
                    treat_temp.index = index
                    treatment = treatment.append(treat_temp)
        
            print(i)
            for j in range(1,7):
                if j == 1:
                    print('Results ' + str(j) + ': ' +str([len(sind_to_del),len(gind_to_del[j]), len(tind_to_del)]) + ' zero rows dropped')
                else:
                    print('Results ' + str(j) + ': ' +str([len(sind_to_del),len(gind_to_del[j])]) + ' zero rows dropped')
            print('-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-')
    
    
    Z = hac.linkage(samples, method='complete', metric=my_metric) # 12min 12s
    clusters = hac.fcluster(Z,0.4,criterion='distance')
    clust_size = Counter(clusters)

    samples_final = pd.DataFrame(columns=samples.columns.values)
    treatment_final = pd.DataFrame(columns=treatment.columns.values)
    genes_final = {}
    pip_final = {}
    for j in range(1,7):
        genes_final[j] = pd.DataFrame(columns=genes[j].columns.values)
        pip_final[j] = pd.DataFrame(columns=genes[j].columns.values)
            
    index = 1
    for clust in clust_size:
        if clust_size[clust] >= len(no_of_iterations)//2:
            counted = 0
            for ind in indices:
                if len(set(np.nonzero(clusters == clust)[0]).intersection(ind)) > 0:
                    counted += 1
            if counted >= len(no_of_iterations)//2:
                temp_samples = samples[pd.Series(clusters==clust, index=samples.index.values)]
                for ind in temp_samples.index.values:
                    # print(stats.pearsonr(temp_samples.iloc[0], temp_samples.loc[ind])[0])
                    if stats.pearsonr(temp_samples.iloc[0], temp_samples.loc[ind])[0] < 0:
                        temp_samples.loc[ind] = temp_samples.loc[ind] * -1
                samples_final.loc[index] = temp_samples.mean(axis=0)
                
                temp_treatment = treatment[pd.Series(clusters==clust, index=treatment.index.values)]
                for ind in temp_treatment.index.values:
                    # print(stats.pearsonr(temp_treatment.iloc[0], temp_treatment.loc[ind])[0])
                    if stats.pearsonr(temp_treatment.iloc[0], temp_treatment.loc[ind])[0] < 0:
                        temp_treatment.loc[ind] = temp_treatment.loc[ind] * -1
                treatment_final.loc[index] = temp_treatment.mean(axis=0)
                
                temp_genes = {}
                for j in range(1,7):
                    temp_genes[j] = genes[j][pd.Series(clusters==clust, index=samples.index.values)]
                    for ind in temp_genes[j].index.values:
                        # print(stats.pearsonr(temp_genes[j].iloc[0], temp_genes[j].loc[ind])[0])
                        if stats.pearsonr(temp_genes[j].iloc[0], temp_genes[j].loc[ind])[0] < 0:
                            temp_genes[j].loc[ind] = temp_genes[j].loc[ind] * -1
                    genes_final[j].loc[index] = temp_genes[j].mean(axis=0)
                    pip_final[j].loc[index] = pip[j][pd.Series(clusters==clust, index=samples.index.values)].median(axis=0)
                index += 1
            else:
                print("Not enough clusters for: " + str(clust))
        
    # Save the data for further use
    samples_final.to_csv('samples_final.csv')
    treatment_final.to_csv('treatment_final.csv')
    pos2treat = {1:'RNA', 2:'citeseq_cell_counts', 3:'cytof_non_cell_counts', 4:'cytof_dep_cell_counts', 5:'luminex_proteins', 6:'masspec_proteins'}
    for j in pos2treat:
        genes_final[j].to_csv('genes_final_'+pos2treat[j]+'.csv')
        pip_final[j].to_csv('pip_final_'+pos2treat[j]+'.csv')
    