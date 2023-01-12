# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd
import copy
import re


import matplotlib.colors as c
import matplotlib.pyplot as plt
import seaborn as sns

def split_chromosomes(df):
    chr_uniq = df.iloc[:,0].value_counts()
    chrs = chr_uniq.keys().to_list()
    count = [chr_uniq[key] for key in chrs]
    print("\n=======================================================================\n")
    print("Number of chromosomes detected : ",len(chrs))
    print("\n=======================================================================\n")
    print("Number of locations in each chromosomes ")
    print(pd.DataFrame(chr_uniq))
    return(chrs,count)

'''
Read genotype data
'''
def ReadGenotype(GenotypeFile):
    df = pd.read_csv(GenotypeFile)
    print("\n=======================================================================\n")
    print("Number of individuals detected : ", df.shape[1]-5)
    print("\n=======================================================================\n")
    print("Individuals: ")
    print(df.columns[5:].to_list())

    chrs,count = split_chromosomes(df)
    return (df,chrs,count)
'''
Count homozygosity and heterozygosity one sample
'''
def homo_count_one_sample(val):
    homo_ref_count = sum([1 for x in val if re.search(r"[A]",x) ])
    homo_alt_count = sum([1 for x in val if re.search(r"[B]",x) ])
    miss_count = sum([1 for x in val if re.search(r"[-U]",x) ])
    hetero_count = sum([1 for x in val if re.search(r"[H]",x) ])
    singleton_count = sum([1 for x in val if re.search(r"[S]",x) ])
    return(homo_ref_count, homo_alt_count, hetero_count, miss_count,singleton_count)

'''
Count homozygosity and heterozygosity all samples
'''
def homo_count_all_samples(snp):

    homo_stats_df = pd.DataFrame(columns=snp.columns[5:],index=['A_count', 'B_count', 'H_count', 'U_count','S_count'])
    for col in snp.columns[5:]:
        homo_stats_df[col] = homo_count_one_sample(snp[col].tolist())
        
    homo_stats_df_T = homo_stats_df.T
    homo_stats_df_T['A_count_perc'] = homo_stats_df_T['A_count']/homo_stats_df_T.sum(axis=1)
    homo_stats_df_T['B_count_perc'] = homo_stats_df_T['B_count']/homo_stats_df_T.sum(axis=1)
    homo_stats_df_T['H_count_perc'] = homo_stats_df_T['H_count']/homo_stats_df_T.sum(axis=1)
    homo_stats_df_T['U_count_perc'] = homo_stats_df_T['U_count']/homo_stats_df_T.sum(axis=1)
    homo_stats_df_T['S_count_perc'] = homo_stats_df_T['S_count']/homo_stats_df_T.sum(axis=1)
    return(homo_stats_df_T)

'''
Display % of different genotypes detected
'''
def plot_homo_count_stats(homo_stats_df_T):
    ax = homo_stats_df_T.loc[:,'A_count_perc':'S_count_perc'].plot.bar(stacked=True,figsize=(20,10),color=['r','b','y','k','pink'])
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    ax.set_xlabel("Samples")
    ax.set_ylabel("Fraction")

'''
Convert genotype to numerical coded matrix
homo parent 1 -> 3
homo parent 2 -> 2
hetero -> 1
missing -> 0
'''
def smooth_genotype_to_numpy_new(snp):
    #samples start at column 5
    snp_tmp = snp.iloc[:,5:].copy(deep=True)
    snp_tmp = snp_tmp.replace(['A'],3)
    snp_tmp = snp_tmp.replace(['B'],2)
    snp_tmp = snp_tmp.replace(['H'],1)
    snp_tmp = snp_tmp.replace(['U'],0)  
    snp_tmp = snp_tmp.replace(['S'],4) # singletons
    return(snp_tmp.to_numpy())


def plotHeatMap(df):
    #'r','b','y','k'
    colors = {"black": 0,"yellow":1, "blue":2, "red":3}
    l_colors = sorted(colors, key=colors.get)
    cMap = c.ListedColormap(l_colors)
    plt.figure(figsize = (20,15))
    ax = sns.heatmap(df,cmap=l_colors, vmin=0, vmax=len(colors))
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([0.5, 1.5, 2.5, 3.5])
    colorbar.set_ticklabels(['Miss', 'HET', 'HOMO1', 'HOMO2'])

def plotHeatMap_with_singletons(df):
    #'r','b','y','k'
    colors = {"black": 0,"yellow":1, "blue":2, "red":3,"pink": 4}
    l_colors = sorted(colors, key=colors.get)
    cMap = c.ListedColormap(l_colors)
    plt.figure(figsize = (20,15))
    ax = sns.heatmap(df,cmap=l_colors, vmin=0, vmax=len(colors))
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([0.5,1.5,2.5,3.5,4.5])
    colorbar.set_ticklabels(['Miss', 'HET', 'HOMO1', 'HOMO2','Singleton'])

