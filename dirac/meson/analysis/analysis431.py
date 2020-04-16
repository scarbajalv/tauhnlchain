#!/usr/bin/env python
# coding: utf-8


from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import ROOT


ifile = ROOT.TFile.Open("/home/sane/mainhnlchain/dirac/120/meson2hnl/data/120-431.root","READ")
tree = ifile.Get("nu")
treesize = tree.GetEntries()
print (tree.GetEntries()/1000000)

data0_12e = np.array([])
data0_14e = np.array([])
data0_16e = np.array([])
for i in range (0,tree.GetEntries()):
    tree.GetEntry(i)
    if getattr(tree,"id") == 12 and getattr(tree,"det_id") == 0 :
        data0_12e = np.append(data0_12e,getattr(tree,"e"))
    if i % 10000 == 0:
    	print (str(1.*i/treesize*100),"%")


datatoy = np.array([1,1,2,1,2,34,1,2,1,2,3,1,1,2,3,1,1])

nbins = 4

plt.hist(datatoy, bins=nbins, histtype="step", density=True, label='lel')
plt.savefig("01.pdf")










