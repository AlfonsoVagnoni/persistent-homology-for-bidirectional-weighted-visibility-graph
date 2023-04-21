# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 13:45:48 2022

@author: Ryl
"""
import scipy.io as scio
import pyflagser as pf
import numpy as np
from scipy import sparse
import math
import os

def get_loc(i,li):
    if np.isinf(i):
        return len(li)
    for index, value in enumerate(li):
        if i < value:
            return index
    else:
        return index + 1

file_dir = 'E:\\mat'
path = []
name = []

for root, dirs, files in os.walk(file_dir):
    for file in files:
        if os.path.splitext(file)[1] == '.mat':
            path.append(os.path.join(root, file))
            name.append(os.path.splitext(file)[0])  

for iii in range(len(path)):

    dataMat = scio.loadmat(path[iii])
    VG=dataMat['WVG']
    numb=len(VG[0])
    daolian=len(VG)
    threshold=np.arange(math.pi/36,math.pi/2+0.01,math.pi/36)
    num_sets = len(threshold)

    homology0 = np.zeros(shape=(numb,daolian,1, num_sets + 1))
    homology1 = np.zeros(shape=(numb,daolian,num_sets, num_sets + 1))
    homology1p = np.zeros(shape=(numb,daolian,num_sets, num_sets + 1))
    En = np.zeros(shape=(numb,daolian,3))

    for e in range(numb):
        for w in range(daolian):
            M = VG[w,e]
            A = sparse.coo_matrix((M[:, 2], (M[:, 0], M[:, 1])),
                                shape=(int(max(M[:, 0]))+1, int(max(M[:, 0]))+1))
            homology = pf.flagser_weighted(A, directed=False, max_dimension=1)
        
            sets0 = homology['dgms'][0]
            plist0 = [] 
            for i in range(len(sets0)):   
                if sets0[i,0]<math.pi/2 and sets0[i,1]<math.pi/2:
                    plist0.append(i)

            sets0=sets0[plist0,:]
            
            for i in range(len(sets0)):
                index0=get_loc(sets0[i][1],threshold)
                homology0[e,w,0,index0] = homology0[e,w,0,index0] + 1
                
            pres0= sets0[:,1]-sets0[:,0]  
            pres0[np.isinf(pres0)]=0
            pres0 = pres0[pres0!=0]
            L0=np.sum(pres0);
            p0=pres0/L0;
            En[e,w,0] = -np.sum(p0*np.log(p0))
            
            sets1 = homology['dgms'][1]
            
            plist = [] 
            for i in range(len(sets1)):   
                if sets1[i,0]<math.pi/2 and sets1[i,1]<math.pi/2:
                    plist.append(i)

            sets1=sets1[plist,:]
            
            pres = sets1[:,1]-sets1[:,0]
            
            for i in range(len(sets1)):
                index1=get_loc(sets1[i][0],threshold)
                index2=get_loc(sets1[i][1],threshold)
                index3=get_loc(pres[i],threshold)
                homology1[e,w,index1,index2] = homology1[e,w,index1,index2] + 1
                homology1p[e,w,index1,index3] = homology1p[e,w,index1,index3] + 1
                
            pres[np.isinf(pres)]=0
            pres = pres[pres!=0]
            L=np.sum(pres);
            p1=pres/L;
            En[e,w,1] = -np.sum(p1*np.log(p1))
            
            juzhen=homology1[e,w,0:18,0:18]
            guiyi=juzhen/np.sum(np.sum(juzhen))
            guiyi=guiyi.flatten()
            guiyinozero=guiyi[guiyi!=0]
            logguiyi=np.log(guiyinozero);
            En[e,w,2] = -np.sum(guiyinozero*logguiyi);
            
    matfile = 'E://mat//ph//{}.mat'.format(name[iii])

    scio.savemat(matfile, {'homology0':homology0,'homology1':homology1,'homology1p':homology1p,'E':En})
