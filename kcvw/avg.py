#!/usr/bin/python

import sys
import numpy as np

#print("There are %d files." % len(sys.argv[1:]))

alldata = [np.loadtxt(f,delimiter=',') for f in sys.argv[1:]]

cols = alldata[0].shape[1]
maxentries =  max( [len(d) for d in alldata])
outdata = np.zeros([maxentries,cols])
#print(outdata.shape)

for i in range(0,maxentries):
    print(",".join(map(str,np.average([d[i] for d in filter( lambda d : len(d) > i,alldata)],axis=0))))
    
