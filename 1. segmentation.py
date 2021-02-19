#import libraries
import os
import skimage as sk
from skimage import io
import numpy as np
from sklearn import mixture
from skimage.morphology import disk
from skimage.filters import median

#*this process can be iterated*

#image import
os.chdir("/path") #Directory path where stacks of u-CT reconstructions with VOI defined are located
sample_name="sample.tif" #name of the sample to read (tif stack format)
I_stack = sk.io.imread(sample_name) # read image using imread function

#segmentation
I_stack_seg=[]
for i in range(I_stack.shape[0]):

    I=I_stack[i]

    # reshape training image from matrix to 1d array:
    I_resh = I.reshape(I.shape[0]*I.shape[1],1).astype(np.uint8) #we are including the 3 dimensions MxNxL in a 1D array

    #delete all the 0 values in the 1D array
    Iresh1=I_resh[I_resh>0] #valors on iresh >0

    #Fitting a gaussian mixture model (GMM) with 4 groups
    gmm4=mixture.GaussianMixture(n_components=4)
    y_pred=gmm4.fit_predict(I_resh.reshape(-1, 1))
    I_segmented=y_pred.reshape(I.shape[0],I.shape[1])

    # Relabeling the pixels according to the sorted means of the GMM model so that 0 corresponds to the lowest mean and so on
    index = np.argsort(gmm4.means_,axis=0)
    I_segmented1 = np.zeros((I.shape[0],I.shape[1]))
    I_segmented1[I_segmented==0]=np.where(index==0)[0][0]
    I_segmented1[I_segmented==1]=np.where(index==1)[0][0]
    I_segmented1[I_segmented==2]=np.where(index==2)[0][0]
    I_segmented1[I_segmented==3]=np.where(index==3)[0][0]

    #append the segmented slice (I_segmented1) to the global segmented matrix(I_stack_seg)
    I_stack_seg.append(I_segmented1)


#image filtering
I_stack_seg_filt=[]

for i in range(I_stack.shape[0]):
    f=I_stack_seg[i].astype(np.uint8)
    f=median(f,disk(5)) # 5pixel (=50um) disc median filter to eliminate the segmentation artefacts
    I_stack_seg_filt.append(f)
