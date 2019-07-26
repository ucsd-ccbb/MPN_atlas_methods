import numpy as np
#import matplotlib
#import matplotlib.pyplot as plt
#import seaborn as sns
import networkx as nx
import pandas as pd
import random
import string
import scipy.stats
import network_prop
import sys
#import visJS2jupyter.visJS_module
#import visJS2jupyter.visualizations

def main(num_reps=10,int_name='GIANT_p3',rand_method = 'degree_ks_test',single_or_double='single'):
    '''
    
    Calculate z-scores for heat propagation
    
    example:
    python AWS_zscore_jamieson.py 1000 STRING degree_binning single
    
    '''
    
    
    print('number of randomizations = '+str(num_reps))
    print('background interactome = ' + int_name)
    print('randomization method = ' + rand_method)
    print('single or double = ' + single_or_double)
    
    num_reps = int(num_reps)
    
    # load interactomes and select focal interactome
    Gint = nx.Graph()
    if int_name=='GIANT_p3':
        G_giant = nx.read_gpickle('interactomes/G_giant_.3.gpickle')
        Gint = G_giant
    elif int_name=='GIANT_p25':
        G_giant = nx.read_gpickle('interactomes/G_giant_.25.gpickle')
        Gint = G_giant
    elif int_name=='GIANT_p2':
        G_giant = nx.read_gpickle('../interactomes/G_giant_.2.gpickle')
        Gint = G_giant
    elif int_name=='GIANT_p15':
        G_giant = nx.read_gpickle('interactomes/G_giant_.15.gpickle')
        Gint = G_giant
    elif int_name=='GIANT_full_p2':
        G_giant = nx.read_gpickle('interactomes/G_giant_full.2.gpickle')
        Gint = G_giant
    elif int_name=='GIANT_kidney_p2':
        G_giant = nx.read_gpickle('interactomes/kidney_0.2.gpickle')
        Gint = G_giant
    elif int_name=='GIANT_heart_p2':
        G_giant = nx.read_gpickle('interactomes/G_giant_heart.2.gpickle')
        Gint = G_giant
    elif int_name=='GIANT_neuron_p3':
        G_giant = nx.read_gpickle('interactomes/G_giant_neuron.3.gpickle')
        Gint = G_giant
    elif int_name=='STRING':
        G_str = nx.read_gpickle('../interactomes/G_str_181022.gpickle')
        Gint = G_str
    elif int_name=='Menche':
        G_menche = nx.read_gpickle('interactomes/G_menche.gpickle')
        Gint = G_menche
    elif int_name=='inBIO':
        G_inbio = nx.read_gpickle('interactomes/G_inbio.gpickle')
        Gint = G_inbio
    elif int_name=='PCnet':
        G_pcnet = nx.read_gpickle('../interactomes/G_PCnet.gpickle')
        Gint = G_pcnet
    elif int_name=='PCnet_pruned':
        G_pcnet = nx.read_gpickle('../interactomes/PCnet_pruned_5MB.gpickle')
        Gint = G_pcnet
    elif int_name=='G_coexp_preA_p8':
        G_coexp = nx.read_gpickle('../interactomes/G_coexp_preA.8.gpickle')
        Gint = G_coexp
    elif int_name=='G_coexp_preU_p8':
        G_coexp = nx.read_gpickle('../interactomes/G_coexp_preU.8.gpickle')
        Gint = G_coexp
    if 'None' in Gint.nodes():
        Gint.remove_node('None')


    # load HC genes
#    AMLvsAN_genes = pd.read_csv('AMLvsAN_seed_genes_190417.tsv',sep='\t',index_col='Unnamed: 0')
#    MFvsAN_genes = pd.read_csv('MFvsAN_seed_genes_190417.tsv',sep='\t',index_col='Unnamed: 0')
#    AMLvsMF_genes = pd.read_csv('AMLvsMF_seed_genes_190417.tsv',sep='\t',index_col='Unnamed: 0')
    
    # updated 5/9/19
    AMLvsAN_genes = pd.read_csv('AMLvsAN_seed_genes_190509.tsv',sep='\t',index_col='Unnamed: 0')
    MFvsAN_genes = pd.read_csv('MFvsAN_seed_genes_190509.tsv',sep='\t',index_col='Unnamed: 0')
    AMLvsMF_genes = pd.read_csv('AMLvsMF_seed_genes_190509.tsv',sep='\t',index_col='Unnamed: 0')
    
    seed_AMLvsAN = AMLvsAN_genes['seed_genes'].tolist()[0].translate(None,string.punctuation).split(' ')
    print(len(seed_AMLvsAN))
    seed_AMLvsAN = list(np.intersect1d(seed_AMLvsAN,Gint.nodes())) # only keep seed genes in the interactome
    print(len(seed_AMLvsAN))
    
    seed_MFvsAN = MFvsAN_genes['seed_genes'].tolist()[0].translate(None,string.punctuation).split(' ')
    print(len(seed_MFvsAN))
    seed_MFvsAN = list(np.intersect1d(seed_MFvsAN,Gint.nodes())) # only keep seed genes in the interactome
    print(len(seed_MFvsAN))
    
    seed_AMLvsMF = AMLvsMF_genes['seed_genes'].tolist()[0].translate(None,string.punctuation).split(' ')
    print(len(seed_AMLvsMF))
    seed_AMLvsMF = list(np.intersect1d(seed_AMLvsMF,Gint.nodes())) # only keep seed genes in the interactome
    print(len(seed_AMLvsMF))
    
    # calculate the z-score
    # calc Wprime from Gint
    Wprime = network_prop.normalized_adj_matrix(Gint,conserve_heat=True)
    
    if single_or_double=='single': # calculate z-scores from a single set of seed genes

        print('calculating z-scores AMLvsAN')
        z_seed_AMLvsAN,Fnew_rand_AMLvsAN = calc_zscore_heat(Gint,Wprime,seed_AMLvsAN,num_reps=num_reps,rand_method=rand_method)
        z_seed_AMLvsAN.to_csv('z_seed_AMLvsAN'+str(num_reps)+'_reps'+int_name+'_'+rand_method+'.tsv',sep='\t')
        
        print('calculating z-scores MFvsAN')
        z_seed_MFvsAN,Fnew_rand_MFvsAN = calc_zscore_heat(Gint,Wprime,seed_MFvsAN,num_reps=num_reps,rand_method=rand_method)
        z_seed_MFvsAN.to_csv('z_seed_MFvsAN'+str(num_reps)+'_reps'+int_name+'_'+rand_method+'.tsv',sep='\t')
        
        print('calculating z-scores AMLvsMF')
        z_seed_AMLvsMF,Fnew_rand_AMLvsMF = calc_zscore_heat(Gint,Wprime,seed_AMLvsMF,num_reps=num_reps,rand_method=rand_method)
        z_seed_AMLvsMF.to_csv('z_seed_AMLvsMF'+str(num_reps)+'_reps'+int_name+'_'+rand_method+'.tsv',sep='\t')


# ----- this part is obsolete -----
#    elif single_or_double=='double': # calculate z-scores from two sets of seed genes:
        
        # --- CNV-DIFF_SPLICE and CNV-METH temorarily commented out ----
        
#        print('calculating CNV-DIFF_SPLICE z-scores')
#        z_CNV_DIFF_SPLICE,Fnew_rand_CNV_DIFF_SPLICE = calc_zscore_heat_double(Gint,Wprime,CNV_HC,DIFF_SPLICE_HC,num_reps=num_reps,rand_method = rand_method)
#        z_CNV_DIFF_SPLICE.to_csv('z_CNV_DIFF_SPLICE'+str(num_reps)+'_reps'+int_name+'_'+rand_method+'.tsv',sep='\t')
#        
#        print('calculating CNV-METH z-scores')
#        z_CNV_METH,Fnew_rand_CNV_METH = calc_zscore_heat_double(Gint,Wprime,CNV_HC,METH_HC,num_reps=num_reps,rand_method = rand_method)
#        z_CNV_METH.to_csv('z_CNV_METH'+str(num_reps)+'_reps'+int_name+'_'+rand_method+'.tsv',sep='\t')
        
#        print('calculating METH-DIFF_SPLICE z-scores')
#        z_METH_DIFF_SPLICE,Fnew_rand_METH_DIFF_SPLICE = calc_zscore_heat_double(Gint,Wprime,METH_HC,DIFF_SPLICE_HC,num_reps=num_reps,rand_method = rand_method)
#        z_METH_DIFF_SPLICE.to_csv('z_METH_DIFF_SPLICE'+str(num_reps)+'_reps'+int_name+'_'+rand_method+'.tsv',sep='\t')
    
    
def calc_zscore_heat(Gint,Wprime,genes_D1,num_reps=10,ks_sig = 0.3,rand_method = 'degree_ks_test'):
    '''
    Helper function to calculate the z-score of heat values from one input seet of genes
    
    rand_method = 'degree_ks_test', or 'degree_binning'.  select the type of randomization
    '''
    seed_D1 = list(np.intersect1d(list(genes_D1),Gint.nodes()))
    Fnew_D1 = network_prop.network_propagation(Gint,Wprime,seed_D1,alpha=.5,num_its=20)
    
    
    Fnew_rand_D1 = np.zeros([num_reps,len(Fnew_D1)])
    if rand_method == 'degree_ks_test':
        for r in range(num_reps):
            if (r%50)==0:
                print(r)
            # UPDATE 8/23/17 -- replace with randomly selecting seed nodes, checking for degree distribution equivalence

            p=0
            # resample until degree distributions are not significantly different
            while p<ks_sig:
                seed_D1_random = Gint.nodes()
                np.random.shuffle(seed_D1_random)
                seed_D1_random = seed_D1_random[0:len(seed_D1)]
                ks_stat,p=scipy.stats.ks_2samp(pd.Series(Gint.degree(seed_D1)),pd.Series(Gint.degree(seed_D1_random)))


            Fnew_rand_tmp = network_prop.network_propagation(Gint,Wprime,seed_D1_random,alpha=.5,num_its=20)
            Fnew_rand_tmp.loc[seed_D1_random]=np.nan # set seeds to nan so they don't bias results
            Fnew_rand_D1[r] = Fnew_rand_tmp.loc[Fnew_D1.index.tolist()]

    elif rand_method == 'degree_binning':
        bins = get_degree_binning(Gint,10)
        min_degree, max_degree, genes_binned = zip(*bins)
        bin_df = pd.DataFrame({'min_degree':min_degree,'max_degree':max_degree,'genes_binned':genes_binned})
        for r in range(num_reps):
            if (r%50)==0:
                print(r)
            # UPDATE 1/30/18 -- sample from degree bins
            seed_D1_random = []
            for g in seed_D1:
                degree_temp = nx.degree(Gint,g)
                # find genes with similar degrees to focal gene degree
                genes_temp = bin_df[(bin_df['min_degree']<=degree_temp)&(bin_df['max_degree']>=degree_temp)]['genes_binned'].tolist()[0]
                np.random.shuffle(genes_temp) # shuffle them
                seed_D1_random.append(genes_temp[0]) # build the seed_D1_random list
                

            Fnew_rand_tmp = network_prop.network_propagation(Gint,Wprime,seed_D1_random,alpha=.5,num_its=20)
            Fnew_rand_tmp.loc[seed_D1_random]=np.nan # set seeds to nan so they don't bias results
            Fnew_rand_D1[r] = Fnew_rand_tmp.loc[Fnew_D1.index.tolist()]


    z_score_D1 = (np.log(Fnew_D1)-np.nanmean(np.log(Fnew_rand_D1),axis=0))/np.nanstd(np.log(Fnew_rand_D1),axis=0)
    
    return z_score_D1, Fnew_rand_D1

# this function obsolete
#def calc_zscore_heat_double(Gint,Wprime,genes_D1,genes_D2,num_reps=10,ks_sig = 0.3,rand_method = 'degree_binning'):
#    '''
#    Helper function to calculate the z-score of heat values from two input sets of genes
#    
#    rand_method = 'degree_binning'.  (this is the only option for now 'degree_ks_test' is deprecated) 
#    '''
#    seed_D1 = list(np.intersect1d(list(genes_D1),Gint.nodes()))
#    Fnew_D1 = network_prop.network_propagation(Gint,Wprime,seed_D1,alpha=.5,num_its=20)
#    
#    seed_D2 = list(np.intersect1d(list(genes_D2),Gint.nodes()))
#    Fnew_D2 = network_prop.network_propagation(Gint,Wprime,seed_D2,alpha=.5,num_its=20)
#    
#    Fnew_both = Fnew_D1*Fnew_D2
#    
#    
#    Fnew_rand_both = np.zeros([num_reps,len(Fnew_both)])
#    if rand_method == 'degree_binning':
#        bins = get_degree_binning(Gint,10)
#        min_degree, max_degree, genes_binned = zip(*bins)
#        bin_df = pd.DataFrame({'min_degree':min_degree,'max_degree':max_degree,'genes_binned':genes_binned})
#        for r in range(num_reps):
#            if (r%50)==0:
#                print(r)
#            # UPDATE 1/30/18 -- sample from degree bins
#            seed_D1_random = []
#            for g in seed_D1:
#                degree_temp = nx.degree(Gint,g)
#                # find genes with similar degrees to focal gene degree
#                genes_temp = bin_df[(bin_df['min_degree']<=degree_temp)&(bin_df['max_degree']>=degree_temp)]['genes_binned'].tolist()[0]
#                np.random.shuffle(genes_temp) # shuffle them
#                seed_D1_random.append(genes_temp[0]) # build the seed_D1_random list
#
#            seed_D2_random = []
#            for g in seed_D2:
#                degree_temp = nx.degree(Gint,g)
#                # find genes with similar degrees to focal gene degree
#                genes_temp = bin_df[(bin_df['min_degree']<=degree_temp)&(bin_df['max_degree']>=degree_temp)]['genes_binned'].tolist()[0]
#                np.random.shuffle(genes_temp) # shuffle them
#                seed_D2_random.append(genes_temp[0]) # build the seed_D1_random list
#                
#
#            Fnew_rand_1 = network_prop.network_propagation(Gint,Wprime,seed_D1_random,alpha=.5,num_its=20)
#            Fnew_rand_1.loc[seed_D1_random]=np.nan # set seeds to nan so they don't bias results
#            
#            Fnew_rand_2 = network_prop.network_propagation(Gint,Wprime,seed_D2_random,alpha=.5,num_its=20)
#            Fnew_rand_2.loc[seed_D2_random]=np.nan # set seeds to nan so they don't bias results
#            
#            Fnew_rand_both[r] = (Fnew_rand_1*Fnew_rand_2).loc[Fnew_D1.index.tolist()]
#
#
#    z_score_both = (np.log(Fnew_both)-np.nanmean(np.log(Fnew_rand_both),axis=0))/np.nanstd(np.log(Fnew_rand_both),axis=0)
#    
#    return z_score_both, Fnew_rand_both

def get_degree_binning(g, bin_size, lengths=None):
    '''
    This function comes from network_utilities.py of emregtoobox.  
    '''
    degree_to_nodes = {}
    for node, degree in g.degree().iteritems():
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    values = degree_to_nodes.keys()
    values.sort()
    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1 
        #print low, high, len(val)
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins



if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])