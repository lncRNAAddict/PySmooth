# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 12:34:23 2022

@author: soibamb
"""

from sklearn.neighbors import KNeighborsClassifier
import pandas as pd
import numpy as np
import copy

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()
        
def smooth_numpy_to_genotype(snp):
    snp_tmp_1 = snp.iloc[:,0:5].copy(deep=True)
    snp_tmp_2 = snp.iloc[:,5:].copy(deep=True)
    snp_tmp_2 = snp_tmp_2.replace([3],'A')
    snp_tmp_2 = snp_tmp_2.replace([2],'B')
    snp_tmp_2 = snp_tmp_2.replace([1],'H')
    snp_tmp_2 = snp_tmp_2.replace([0],'U')  # missing items
    snp_tmp_2 = snp_tmp_2.replace([-1],'U') # singletons
    snp_tmp = pd.concat([snp_tmp_1,snp_tmp_2],axis=1)
    return(snp_tmp)



'''
smooth_df: chrom,pos,locus,ref base, alt base, samples .....
'''
def FillMissingGenotype(SMOOTH_df,num_neighbors = 30):
    
    
    #SMOOTH_df = smooth_numpy_to_genotype(SMOOTH_df)
    X=copy.deepcopy(SMOOTH_df)
    XX = copy.deepcopy(X)
    REF=[ref + ref for ref in XX['REF'].tolist()]
    ALT=[alt + alt for alt in XX['ALT'].tolist()]
    REF_ALT = [ ref + alt for ref,alt in zip(XX['REF'].tolist(),XX['ALT'].tolist() ) ]
    #print(REF)
    #print(ALT)
    #print(REF_ALT)
    #printProgressBar(0, l, prefix = 'Number of Samples:', suffix = 'Complete', length = 50)

    for col in SMOOTH_df.iloc[:,5:].columns:
        #print 
        # non missing items
        X_train = copy.deepcopy(SMOOTH_df.loc[ (X[col]!='U') & (X[col]!='S'),:].iloc[:,1].to_numpy())
        #print(X_train.shape)
        y_train = copy.deepcopy(SMOOTH_df.loc[ (X[col]!='U') & (X[col]!='S'),col])
        #print(y_train)
        #missing items
        X_test = copy.deepcopy(SMOOTH_df.loc[ (X[col]=='U') | (X[col]=='S'),:].iloc[:,1].to_numpy())
        #print(X_test.shape)
        if X_test.shape[0] > 0:
          #print(X_train.shape)
          #print(X_test.shape)
          #knn model based on non missing items
          neigh = KNeighborsClassifier(n_neighbors=num_neighbors)
          neigh.fit(np.reshape(X_train,(X_train.shape[0],1)), y_train)
          #predict the genotype of missing items
          y_test = neigh.predict(np.reshape(X_test,(X_test.shape[0],1)))
        
          #combine non missing and predicted missing items
          pos = np.concatenate((X_train,X_test))
          genotype = np.concatenate((y_train,y_test))
          df_tmp = pd.DataFrame(columns=['POS','genotype'])
          df_tmp['POS']=pos
          df_tmp['genotype'] = genotype
          #order based on position in the chromosome
          df_tmp = df_tmp.sort_values(by='POS')
          X[col] = df_tmp['genotype'].to_list()
          val = df_tmp['genotype'].to_list()
          ## get double notation AA, AT, etc
          '''
          for i in range(0,len(val)):
            if val[i] == 'A':
                val[i] = REF[i]
            elif val[i] == 'B':
                val[i] = ALT[i]
            else:
                val[i] = REF_ALT[i]
          XX[col] = val        
'''
    return (X,XX)


    
