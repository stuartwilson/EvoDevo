import numpy as np
import pylab as pl
import h5py
import sys

if (len(sys.argv)<2):
    fname = 'logs/out.h5'
else:
    fname = sys.argv[1]
    
F = h5py.File(fname)
meanF = F['meanF'][:]
maxF = F['maxF'][:]
minF = F['minF'][:]
varF = F['varF'][:]
AmeanF = F['AmeanF'][:]
AmaxF = F['AmaxF'][:]
AminF = F['AminF'][:]
AvarF = F['AvarF'][:]
BmeanF = F['BmeanF'][:]
BmaxF = F['BmaxF'][:]
BminF = F['BminF'][:]
BvarF = F['BvarF'][:]
CmeanF = F['CmeanF'][:]
CmaxF = F['CmaxF'][:]
CminF = F['CminF'][:]
CvarF = F['CvarF'][:]
DmeanF = F['DmeanF'][:]
DmaxF = F['DmaxF'][:]
DminF = F['DminF'][:]
DvarF = F['DvarF'][:]
G = F['G'][:]
S = F['species'][:]
SC = F['SC'][:]
F.close()

gens = len(meanF)
p=1000
G = G.reshape(p,int(len(G)/p))
SC = SC.reshape(gens,int(len(SC)/gens))


U = np.zeros([1,G.shape[1]])
U[0,:] = G[0]
counts = [1]
for i in range(1,p):
    different = True
    for j in range(U.shape[0]):
        identical = True
        for k in range(G.shape[1]):
            if(G[i,k] != U[j][k]):
                identical=False
                break
        if(identical):
            different = False
            counts[j] +=1
            break
    if(different):
        U = np.vstack([U,G[i]])
        counts = np.hstack([counts,1])

print(U)
print(counts)
print(np.sum(counts))




all = False



F = pl.figure()
f = F.add_subplot(111)
f.pcolor(U)
f.set_title('Unique genomes')

F = pl.figure(figsize=(7,9))
f = F.add_subplot(511)
f.plot(meanF)
if(all):
    f.plot(maxF)
    f.plot(minF)
    f.plot(varF)
    f.legend(['mean','max','min','var'],loc='right')
f.set_ylabel('fitness')


f = F.add_subplot(512)
f.plot(AmeanF)
if(all):
    f.plot(AmaxF)
    f.plot(AminF)
    f.plot(AvarF)
f.set_ylabel('genotype bias')

f = F.add_subplot(513)
f.plot(DmeanF)
if(all):
    f.plot(DmaxF)
    f.plot(DminF)
    f.plot(DvarF)
    f.legend(['mean','max','min','var'],loc='right')
f.set_ylabel('genotype complexity')


f = F.add_subplot(514)
f.plot(BmeanF)
if(all):
    f.plot(BmaxF)
    f.plot(BminF)
    f.plot(BvarF)
f.set_ylabel('phenotype bias')


f = F.add_subplot(515)
f.plot(CmeanF)
if(all):
    f.plot(CmaxF)
    f.plot(CminF)
    f.plot(CvarF)
f.set_ylabel('cycle length')

F = pl.figure()
f = F.add_subplot(211)
f.plot(S)
f.set_ylabel('species ('+str(S[-1])+')')

f = F.add_subplot(212)
f.plot(SC)
f.set_ylabel('species size')


pl.show()

