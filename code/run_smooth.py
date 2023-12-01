# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 14:40:05 2022

@author: soibamb
"""
import utilities
import smooth
import ImputeMissingGenotype
import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt
from optparse import OptionParser


'''
Parsing input section
'''
parser = OptionParser()
parser.add_option("-i", "--input", action="store", type="string", dest="genotype_file")
parser.add_option("-o", "--output", action="store", type="string", dest="output_prefix",default="test")
parser.add_option("-c", "--chr", action="store", type="string", dest="chr_list",default="all")
parser.add_option("-l","--lower", type="float", dest="lower_th",default = 0.70)
parser.add_option("-u","--upper",type="float", dest="upper_th",default = 0.98)
parser.add_option("-g","--gap", type="float", dest="gap",default = 0.02)
parser.add_option("-k",type="int", dest="n_neighbors",default = 30)

(options, args) = parser.parse_args()
print(options.genotype_file)
print(options.output_prefix)
print(options.chr_list)
print(options.lower_th)
print(options.upper_th)
print(options.gap)

genotype_file = options.genotype_file
output_prefix = options.output_prefix
chr_list = options.chr_list
lower_th = options.lower_th
upper_th = options.upper_th
gap = options.gap
n_neighbors = options.n_neighbors



#Read Raw data
df,chrs,count = utilities.ReadGenotype(genotype_file) # SMOOTH_df_85_flies_ref_corrected_A_B_H_U.csv
column_names = df.columns.to_list()
#print(chrs)
#print(count)
if chr_list == "all":
    chr_list = chrs
else:
    chr_list = chr_list.split(",")
    
print("\n=======================================================================\n")
print("user chose to analyze the chromosomes: " , chr_list)
print("\n=======================================================================\n")
'''
For each chromosome of interest, report statistics and plots and heatmaps of the raw data
'''
print("Generating genotype statistics, plots, and heatmaps ")

for chromosome in chr_list:
    print(" doing " + chromosome)
    df1 = copy.deepcopy(df.loc[df[column_names[0]] == chromosome,:])
    #print(df1)
    df_stats = utilities.homo_count_all_samples(df1)
    utilities.plot_homo_count_stats(df_stats)
    output_file = output_prefix + "_" + chromosome
    plt.savefig(output_file+".stats.png") 
    df_stats.to_csv(output_file+".stats.csv",index=True) 
    X1 = utilities.smooth_genotype_to_numpy_new(df1)
    
    X1 = pd.DataFrame(X1,columns=df1.columns[5:],index=df1[column_names[1]])
    utilities.plotHeatMap(X1)
    plt.savefig(output_file+".heatmap.png")


    print("\n=======================================================================\n")

    '''
    Remove singletons, report statistics and plots and heatmaps
    '''
    print("Identifying Singletons\n")
    print("Missing values and Singletons will be imputed\n")


    
#for chromosome in chr_list:
    
#    print("doing " + chromosome)
#    print("\n")
    input_file = genotype_file
    output_file = output_prefix + "_" + chromosome
    
    
    SINGLETONS_df,SMOOTHED_df = smooth.run_smooth_new(input_file,chromosome,upper_th,lower_th,gap)
    
    
    
    #SINGLETONS_df.to_csv(output_file +"_singletons_stats.csv", index=False)
    SINGLETONS_df_ = SINGLETONS_df.drop(columns=['threshold'])
    SINGLETONS_df_ = SINGLETONS_df_.groupby(['Singletons']).sum().T
    
    SMOOTHED_df.to_csv(output_file+"_singletons.csv",index=False)
    
    SINGLETONS_df_['S_count'] = SINGLETONS_df_.sum(axis=1)
    SINGLETONS_df_['A_count'] = df_stats['A_count'] - SINGLETONS_df_['A -> Singletons']
    SINGLETONS_df_['B_count'] = df_stats['B_count'] - SINGLETONS_df_['B -> Singletons']
    SINGLETONS_df_['H_count'] = df_stats['H_count'] - SINGLETONS_df_['H -> Singletons']
    SINGLETONS_df_['U_count'] = df_stats['U_count']
    
    test_df = SINGLETONS_df_.loc[:,'S_count':'U_count']

    SINGLETONS_df_['A_count_perc'] = SINGLETONS_df_['A_count']/test_df.sum(axis=1)
    SINGLETONS_df_['B_count_perc'] = SINGLETONS_df_['B_count']/test_df.sum(axis=1)
    SINGLETONS_df_['H_count_perc'] = SINGLETONS_df_['H_count']/test_df.sum(axis=1)   
    SINGLETONS_df_['U_count_perc'] = SINGLETONS_df_['U_count']/test_df.sum(axis=1)
    SINGLETONS_df_['S_count_perc'] = SINGLETONS_df_['S_count']/test_df.sum(axis=1)

    SINGLETONS_df_.to_csv(output_file +"_singletons_stats.csv")
    utilities.plot_homo_count_stats(SINGLETONS_df_)
    plt.savefig(output_file +"_singletons_stats.png")
    '''
    Generate plots, statistics after smoothing
    '''
    X1 = utilities.smooth_genotype_to_numpy_new(SMOOTHED_df)
    X1=pd.DataFrame(X1,columns=SMOOTHED_df.columns[5:],index=SMOOTHED_df.iloc[:,1])
    utilities.plotHeatMap_with_singletons(X1)
    plt.savefig(output_file+"_singletons_heatmap.png")
    ##############################################################################
    
    
    X,XX = ImputeMissingGenotype.FillMissingGenotype(SMOOTHED_df,num_neighbors = n_neighbors)
    X.to_csv(output_file +"_imputed.csv",index=False)
    
    X1 = utilities.smooth_genotype_to_numpy_new(X)
    X1=pd.DataFrame(X1,columns=X.columns[5:],index=X.iloc[:,1])
    
    utilities.plotHeatMap(X1)
    plt.savefig(output_file +"_imputed_heatmap.png")
    
    df_stats = utilities.homo_count_all_samples(X)
    df_stats.to_csv(output_file +"_imputed_stats.csv",index=False)
    
    utilities.plot_homo_count_stats(df_stats)
    plt.savefig(output_file +"_imputed_stats.png")

    