#Import libraries
import numpy as np
from skimage.morphology import binary_erosion, binary_dilation, ball
import pandas as pd

#Create a pandas for storing all the calculations
df = pd.DataFrame(columns=['sample name','bone surface #px', 'available biomat surface #px', 'bone surface in available surface'])

#*this process can be iterated*

#Use the segmented files from section 1 (variable name: I_stack_seg_filt)
sample=I_stack_seg_filt
sample=np.int8(sample) #convert the sample from int64 (64 bit integrer) to int8 (8 bit integrer)

#Create binary arrays with the tissue of interest converted to 1's and the rest to 0's
bone=np.where(sample==2,1,0) 
bone_filtered=binary_dilation(binary_erosion(bone,ball(2)), ball(2)) 
biomat=np.where(sample==3,1,0) 

#Calculate the surface (number of px of the interface) bone/scaffold
contact_region=(binary_dilation(bone_filtered, ball(1)))*biomat #first we dilate the bonefiltered 1px then we multipy times the biomat (only px in common remain)

#Calculate the total surface available in the biomaterial
#  1.Dilate bone or biomat and substract the original biomat to the result 
surface_biomat_available=(binary_dilation(biomat))-biomat
#  2.Quantify the percentages of bone surface/surface available
pixels_in_contact_region=contact_region.sum()#count the number of ones in a binary np array
pixels_available=surface_biomat_available.sum()
percentage_of_bone_surface_in_contact_with_biomaterial=pixels_in_contact_region/pixels_available

#Print the results
print(f"Percentage of biomaterial surface covered by bone: {round(percentage_of_bone_surface_in_contact_with_biomaterial*100,2)} %")

#Append the results to the pandas
df = df.append({
    'sample name': "sample name", 
    'bone surface #px': pixels_in_contact_region, 
    'available biomat surface #px': pixels_available, 
    'bone surface in available surface': percentage_of_bone_surface_in_contact_with_biomaterial
}, ignore_index=True)

df