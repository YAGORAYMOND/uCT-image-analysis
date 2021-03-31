#Import libraries
import pandas as pd
from skimage.morphology import erosion, dilation,disk #per fer els discs
from operator import truediv #to operate with array divisions
import math

#Create a pandas that collects each analysed sample
df=pd.DataFrame(columns=[
    'Sample family',
    'Sample name', 
    'peel layer', 
    'peel distance (um)', 
    'peel px','bone px',  
    'bone faction (%)'
])

#*this process can be iterated*

#Use the segmented files from section 1 (variable name: I_stack_seg_filt)
sample=I_stack_seg_filt
sample=np.int8(sample) #convert the sample from int64 (64 bit integrer) to int8 (8 bit integrer)

#Create binary arrays with the tissue of interest converted to 1 and the rest to 0
bone=sample==2 
bone=binary_dilation(binary_erosion(bone,ball(2)), ball(2))
biomat=sample==3
voi=sample!=0

#Define peel settings
n_peels=20 #Define here the number of layers. With a inter strand gap of 250um and a pixel size of 10um we want to analise 125um which corresponds to 13 peels
px_peel=1 #Define the number of pixels of each peel

#Create two empty arrays where adding the pixels counted in every z. Each index of the list corresponds to a different peel. 
peel_count=[0]*n_peels
bone_peel_count=[0]*n_peels

for z in range(sample.shape[0]): #Pass through all Zs of the stack

    #Select the images
    bone_z=bone[z]
    biomat_z=biomat[z]
    voi_z=voi[z]

    biomat_dilate=biomat_z #biomat_dilate variable is defined as biomat_z only for the first iteration

    for p in range(n_peels):

        biomat_dilate_new=dilation(biomat_dilate, disk(px_peel))#dilates biomat_dilate per crear un nou disc

        peel_and_borders= np.logical_xor(biomat_dilate_new,biomat_dilate) # boolean subtraction
        peel=np.logical_and(voi_z,peel_and_borders)#boolean and (i.e., everything inside the voi)

        bone_peel=np.logical_and(bone_z, peel) #determines which part of the peel corresponds to bone

        #Add to the count arrays the new pixels calculations for each peel (Array index) and at each Z.
        peel_count[p]+=peel.sum()
        bone_peel_count[p]+=bone_peel.sum()

        #The outer border of the peel is now defined as the inner border of the followig peel
        biomat_dilate=biomat_dilate_new        

#Calculate the bone fraction for each peel
bone_fraction =list(map(truediv,bone_peel_count,peel_count))
bone_fraction = [0 if math.isnan(x) else x for x in bone_fraction] #and just in case of having NaN values result of 0/0, we replace them by 0    

#Register the results in the pandas
for p in range(n_peels):
    df=df.append({
        'Sample family':"sample family",
        'Sample name': "sample name",
        'peel layer': p ,
        'peel distance (um)':(p+0.5)*px_peel*10, #defines the peel midpoint in um (knowing that the CT resolution in 10um/px)
        'peel px': peel_count[p],
        'bone px':bone_peel_count[p], #each disc has 10 px (at a resolution of 0.01mm/px corresponds to 0.1mm) 
        'bone faction (%)':(bone_peel_count[p]/peel_count[p])
        },ignore_index=True)
    
df
