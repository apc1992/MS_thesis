#!usr/bin/env/python3

# this script is intended to shortly normalize an input count matrix by applying Z-score to each value of the matrix. Input matrix must be in the working directory #

import numpy as np
import scipy.stats as stats

file=open("tpm_significant.csv","r")
new=open("Z-score.txt","w")

def zscore():
    data=np.genfromtxt(file,delimiter='\t', usecols = np.arange(1,24))
    return stats.zscore(data, axis=1)
listZscore = zscore().tolist()
for itemList in listZscore:
    new.write('\n')
    for item in itemList:
        new.write(str(item)+'\t')


