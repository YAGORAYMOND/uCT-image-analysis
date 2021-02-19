#import libraries
import numpy as np
import matplotlib.pyplot as plt
import glob, os


#Import segmented files
os.chdir("/Users/yago/Desktop/python segmentation CONILLS/04. SEGMENTED/NUMPY/filtered")
#os.chdir("/Path") #Directory path containing all the segmented stacks recorded in numpy format for one single condition
samples=[]
for sample_name in glob.glob("*.npy"):
        #load every sample and store it in the samples array
        #sample = np.load(f"/Path/{sample_name}")
        sample = np.load(f"/Users/yago/Desktop/python segmentation CONILLS/04. SEGMENTED/NUMPY/filtered/{sample_name}")
        samples.append(sample)


#find the dimensions of the biggest sample
max_x=0
max_y=0
for sample in samples:
    if sample.shape[1]>max_x:
        max_x=sample.shape[1]
    if sample.shape[2]>max_y:
        max_y=sample.shape[2]


#create a new np array
void_px_sum= np.zeros([max_x,max_y])
bone_px_sum= np.zeros([max_x,max_y])


#add in the array all the centered samples
for sample in samples:

    void=np.where(sample==1,1,0)
    bone=np.where(sample==2,1,0)

    for i in range (sample.shape[1]):
        for j in range (sample.shape[2]):
            i0=int((max_x-sample.shape[1])/2)
            j0=int((max_y-sample.shape[2])/2)
            void_px_sum [i0+i,j0+j] += void[:,i,j].sum()
            bone_px_sum [i0+i,j0+j] += bone[:,i,j].sum()



#calculate for every "pixel" in the 2D Z projection arry calculate the bone fraction
colormap= np.zeros([max_x,max_y])

for i in range (colormap.shape[0]):
    for j in range(colormap.shape[1]):

        px_void=void_px_sum[i,j]
        px_bone=bone_px_sum[i,j]

        if px_bone!=0 and px_void!=0:
            colormap[i,j]=(px_bone/(px_void+px_bone))*100


#resulting colormap
plt.imshow(colormap,cmap=plt.cm.jet, vmin=0, vmax=100)
