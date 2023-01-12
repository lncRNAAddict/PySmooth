# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 14:32:09 2022

@author: soibamb
"""

import pandas as pd
import numpy as np
import copy

w1= [0.059,0.082,0.112,0.151,0.202,0.265,0.342,0.433,0.537,0.647,0.758,0.857,0.934,0.981,0.998] 
w1=np.array(w1+list(reversed(w1)))

# In[ ]:

def split_chromosomes(df):
    chr_uniq = df.iloc[:,0].value_counts()
    chrs = chr_uniq.keys().to_list()
    count = [chr_uniq[key] for key in chrs]
    print("Number of chromosomes detected : ",len(chrs))
    #display(pd.DataFrame(chr_uniq,columns=['chromosome','count']))
    return(chrs,count)

# In[ ]:

def smooth_numpy_to_genotype(snp):
    snp_tmp_1 = snp.iloc[:,0:5].copy(deep=True)
    snp_tmp_2 = snp.iloc[:,5:].copy(deep=True)
    snp_tmp_2 = snp_tmp_2.replace([3],'A')
    snp_tmp_2 = snp_tmp_2.replace([2],'B')
    snp_tmp_2 = snp_tmp_2.replace([1],'H')
    snp_tmp_2 = snp_tmp_2.replace([0],'U')  # missing items
    snp_tmp_2 = snp_tmp_2.replace([-1],'S') # singletons
    snp_tmp = pd.concat([snp_tmp_1,snp_tmp_2],axis=1)
    return(snp_tmp)

# In[ ]:

def smooth_genotype_to_numpy_new(snp):
    #samples start at column 5
    snp_tmp = snp.iloc[:,5:].copy(deep=True)
    snp_tmp = snp_tmp.replace(['A'],3)
    snp_tmp = snp_tmp.replace(['B'],2)
    snp_tmp = snp_tmp.replace(['H'],1)
    snp_tmp = snp_tmp.replace(['U'],0)  
    return(snp_tmp.to_numpy())


# In[ ]:

def smooth_genotype_one_sample_new(val,threshold,w):
    '''
    ref homo = 3
    alt homo = 2
    hetero = 1
    U = 0
    '''
    val1 = copy.deepcopy(val)
    #val2 = copy.deepcopy(val) # to store
    #REMOVED = 0
    REMOVED={}
    genotypes =[3,2,1,0]
    for key in genotypes:
        REMOVED[key] = 0
        
    for i in range(15,val1.shape[0]-15):
        ori_val = val1[i]
        #print("ori_val",ori_val)
        

        if ori_val != 0:
            t = copy.deepcopy(val1[i-15:i+16])
            t = np.delete(t,15)
            #print("t",t)
            w2=np.where(t!=0,w,0)
            #print("w",w2)
            y = np.sum(w2)

            w2=np.where(t == ori_val,w2,0)
            #print("w",w2)
            t = np.where(t == ori_val ,1,0)
            #print("t",t)
            x = np.sum(np.multiply(t,w2))
            res = x
            if y == 0:
                res = 0
            else:
                res = x/y
            
            if abs(1 - res) > threshold:
                '''
                print(val1[i-15:i+16])
                print("t",t)
                print("w",w2)
                print("res",res)
                print("1-res",1-res)
                '''
                val1[i] = 0 
                #val2[i] = -1 # singleton is assigned -1
                #print("removed")
                
                REMOVED[ori_val] = REMOVED[ori_val] + 1
        else:
            REMOVED[0] = REMOVED[0] + 1
    #print("removed entries ", REMOVED)
    Singletons=[]
    for key in [3,2,1]:
        Singletons.append(REMOVED[key])
    return (val1,Singletons)


import copy
import time

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

def smooth_genotype_all_samples(VAL,threshold,w):
    VAL1 = copy.deepcopy(VAL)
    REMOVED=np.zeros((3,VAL1.shape[1]))
    #genotypes =['A','B','H']
    #for key in genotypes:
    #    REMOVED[key] = 0
        
    l=VAL1.shape[1]
    #printProgressBar(0, l, prefix = 'Number of Samples:', suffix = 'Complete', length = 100)
    for col in range(0,VAL1.shape[1]):
        V = copy.deepcopy(VAL1[:,col])
        (T,removed) = smooth_genotype_one_sample_new(V,threshold,w)
        VAL1[:,col] = copy.deepcopy(T)
       # VAL1[:,col] = copy.deepcopy(T1)
        REMOVED[:,col] = np.array(removed) 
    
        #progress bar
        #printProgressBar(col + 1, l, prefix = 'Number of Samples:', suffix = 'Complete', length = 100)
    return (VAL1,REMOVED)

def generate_thresholds(begin_threshold,ending_threshold,gap):
    thresholds = [begin_threshold]
    th = begin_threshold 
    while th > ending_threshold:
        th = th - gap
        thresholds.append(th)
    print("Thresholds for SMOOTHING:\n",thresholds)
    return(thresholds)

# In[1]:
'''

SMOOTH_df = pd.read_csv('SMOOTH_df_85_flies_A_B_H_U.csv')
columns=SMOOTH_df.columns
display(SMOOTH_df)

chrs,count=split_chromosomes(SMOOTH_df)

thresholds = generate_thresholds(0.98,0.7,0.02)

SMOOTH_df_2L = SMOOTH_df.loc[SMOOTH_df.iloc[:,0] == 'CM002910.1',:]
SMOOTH_df_2L_num = smooth_genotype_to_numpy_new(SMOOTH_df_2L)
SMOOTH_df_2L_num_iter=copy.deepcopy(SMOOTH_df_2L_num)

th=thresholds[0]
print("Smoothing with threshold ",th)
SMOOTH_df_2L_num_iter,SINGLETONS = smooth_genotype_all_samples(SMOOTH_df_2L_num_iter,th,w1)
display(pd.DataFrame(SINGLETONS,index=['3','2','1','0'],columns=SMOOTH_df_2L.columns[5:]))

for th in thresholds[1:]:
    print("Smoothing with threshold ",th)
    #removed_df = pd.DataFrame(index=SMOOTH_df_2L.columns[5:],columns=['A','B','H','U'])
    SMOOTH_df_2L_num_iter, REMOVED = smooth_genotype_all_samples(SMOOTH_df_2L_num_iter,th,w1)
    SINGLETONS = np.concatenate((SINGLETONS,REMOVED)) 
'''
# In[3]:    
def run_smooth(genotypeFile,chromosome_list,begin_threshold,end_threshold,threshold_gap):
    SMOOTH_df = pd.read_csv(genotypeFile)
    #display(SMOOTH_df)
    print("\n=======================================================================\n")

    chrs,count=split_chromosomes(SMOOTH_df)
    
    print("\n=======================================================================\n")

    thresholds = generate_thresholds(begin_threshold,end_threshold,threshold_gap)
    print("\n=======================================================================\n")

    print("selected chromosomes are ",'\n'.join(chromosome_list) )
    for CHR in chromosome_list:
        print("\n=======================================================================\n")
        print("Doing chromsome ",CHR)

        SMOOTH_df_CHR = copy.deepcopy(SMOOTH_df.loc[SMOOTH_df.iloc[:,0] == CHR,:])
        SMOOTH_df_CHR_num = smooth_genotype_to_numpy_new(SMOOTH_df_CHR)
        SMOOTH_df_CHR_num_iter=copy.deepcopy(SMOOTH_df_CHR_num)
        
        th=thresholds[0]
        print("Smoothing with threshold ",th)
        SMOOTH_df_CHR_num_iter,SINGLETONS = smooth_genotype_all_samples(SMOOTH_df_CHR_num_iter,th,w1)

        #display(pd.DataFrame(SINGLETONS,index=['3','2','1','0'],columns=SMOOTH_df_CHR.columns[5:]))

        for th in thresholds[1:]:
            print("Smoothing with threshold ",th)
            #removed_df = pd.DataFrame(index=SMOOTH_df_2L.columns[5:],columns=['A','B','H','U'])
            SMOOTH_df_CHR_num_iter, REMOVED = smooth_genotype_all_samples(SMOOTH_df_CHR_num_iter,th,w1)
            SINGLETONS = np.concatenate((SINGLETONS,REMOVED))
            #display(pd.DataFrame(SINGLETONS,columns=SMOOTH_df_CHR.columns[5:]))
        np.savetxt(genotypeFile+"_singletons_"+CHR + ".csv", SINGLETONS, delimiter=",")
        SMOOTHED_df = pd.DataFrame(SMOOTH_df_CHR_num_iter,columns=SMOOTH_df_CHR.columns[5:])
        for col in reversed(SMOOTH_df_CHR.columns[0:5]):
            SMOOTHED_df.insert(0,col,SMOOTH_df_CHR.loc[:,col].to_list())
        SMOOTHED_df.to_csv(genotypeFile+"_SMOOTH_"+CHR + ".csv",index=False)
    return (SINGLETONS,SMOOTHED_df)
    


def run_smooth_new(genotypeFile,CHR,begin_threshold,end_threshold,threshold_gap):
    

    
    
    SMOOTH_df = pd.read_csv(genotypeFile)
    thresholds = generate_thresholds(begin_threshold,end_threshold,threshold_gap)
    print("\n=======================================================================\n")
    #print("Doing chromsome ",CHR)

    SMOOTH_df_CHR = copy.deepcopy(SMOOTH_df.loc[SMOOTH_df.iloc[:,0] == CHR,:])
    SMOOTH_df_CHR_num = smooth_genotype_to_numpy_new(SMOOTH_df_CHR) # change A,B,H,U to 3,2,1,0 code
    
    ''' 
    padding with zeros at the ends of the chromosome
    '''
    SMOOTH_df_CHR_num = np.pad(SMOOTH_df_CHR_num,((15,15),(0,0)))
    SMOOTH_df_CHR_num_iter=copy.deepcopy(SMOOTH_df_CHR_num)
        
    th=thresholds[0]
    print("Smoothing with threshold ",th)
    SMOOTH_df_CHR_num_iter,SINGLETONS = smooth_genotype_all_samples(SMOOTH_df_CHR_num_iter,th,w1)
    test = ['A -> Singletons','B -> Singletons','H -> Singletons']
    column_1 = test
    column_threshold=[th,th,th]
    for th in thresholds[1:]:
            print("Smoothing with threshold ",th)
            #removed_df = pd.DataFrame(index=SMOOTH_df_2L.columns[5:],columns=['A','B','H','U'])
            SMOOTH_df_CHR_num_iter, REMOVED =smooth_genotype_all_samples(SMOOTH_df_CHR_num_iter,th,w1)
            SINGLETONS = np.concatenate((SINGLETONS,REMOVED))
            column_1 = column_1 + test
            column_threshold = column_threshold + [th,th,th]
            #display(pd.DataFrame(SINGLETONS,columns=SMOOTH_df_CHR.columns[5:]))
    
    SINGLETONS_df = pd.DataFrame(SINGLETONS,columns=SMOOTH_df.columns[5:])
    SINGLETONS_df.insert(loc = 0, column='threshold',value=column_threshold)
    SINGLETONS_df.insert(loc = 0, column='Singletons',value=column_1)
    #SINGLETONS_df.to_csv(genotypeFile+"_singletons_"+ CHR + ".csv", index=False)
    #np.savetxt(genotypeFile+"_singletons_"+ CHR + ".csv", SINGLETONS, delimiter=",")
     
    #convert 3,2,1,0 to A,B,H,U, S(singleton)
    SMOOTH_df_CHR_num_singletons=copy.deepcopy(SMOOTH_df_CHR_num)
    x=np.not_equal(SMOOTH_df_CHR_num,SMOOTH_df_CHR_num_iter)
    SMOOTH_df_CHR_num_singletons[np.where(x)]=-1
    #SMOOTH_df_CHR_num_singletons=smooth_numpy_to_genotype(SMOOTH_df_CHR_num_singletons)    
    
    '''
    remove padding
    '''
    SMOOTH_df_CHR_num_singletons = SMOOTH_df_CHR_num_singletons[15:SMOOTH_df_CHR_num_singletons.shape[0]-15,:]
    
    SMOOTHED_df = pd.DataFrame(SMOOTH_df_CHR_num_singletons,columns=SMOOTH_df_CHR.columns[5:])
    
    for col in reversed(SMOOTH_df_CHR.columns[0:5]):
            SMOOTHED_df.insert(0,col,SMOOTH_df_CHR.loc[:,col].to_list())

    SMOOTHED_df=smooth_numpy_to_genotype(SMOOTHED_df)
    #SMOOTHED_df.to_csv(genotypeFile+"_SMOOTH_"+CHR + ".csv",index=False)
    return (SINGLETONS_df,SMOOTHED_df)
    
    
#SINGLETONS,SMOOTH_df_CHR_num_iter = run_smooth('SMOOTH_df_85_flies_A_B_H_U.csv',['CM002916.1'],0.98,0.70,0.02)










