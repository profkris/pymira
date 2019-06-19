# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 07:41:20 2019

@author: simon
"""
import numpy as np
import cv2

path = r'C:\Users\simon\Dropbox\Rift\Royal Society Scripted\Panorama_grabs'
ofile = r'C:\Users\simon\Dropbox\output.mp4'

nloop = 10 # number of loops
rf = 2 #None # size sclaing factor

import os
imFiles = []
testSize = True
for file in os.listdir(path):
    if file.endswith(".bmp"):
        imFiles.append(os.path.join(path,file))
        if testSize:
            frame0 = cv2.imread(os.path.join(imFiles[0]))
            testSize = False

orSize = np.asarray(frame0.shape)
newSize = orSize.copy()
if rf is not None:
    newSize[0:2] = (newSize[0:2] / float(rf)).astype('int')

fourcc = cv2.VideoWriter_fourcc(*'MP4V')
#ofile = os.path.join(path,'output.mp4')
out = cv2.VideoWriter(ofile, fourcc, 20.0, (newSize[0],newSize[1]))

nfile = len(imFiles)
for j in range(nloop):
    for i,imFile in enumerate(imFiles):
        print('Loop {} of {}, reading/writing image {} of {} (size {}x{})'.format(j+1,nloop,i+1,nfile,newSize[0],newSize[1]))
        frame = cv2.imread(os.path.join(imFile))
        if rf is not None:
            frame = cv2.resize(frame, (newSize[0],newSize[1]))#, interpolation = inter)
        out.write(frame)

print('Written to {}'.format(ofile))
out.release()
cv2.destroyAllWindows()