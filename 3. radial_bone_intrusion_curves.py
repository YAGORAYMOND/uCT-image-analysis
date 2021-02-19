#import libraries
import pandas as pd
from skimage.morphology import erosion,disk
from operator import truediv
import math

#create a dataframe to record the outcome
df=pd.DataFrame(columns=[
        'Sample family',
        'Sample name',
        'peel layer',
        'Radium (mm)',
        'available volume (px)',
        'bone volume (px)',
        'bone faction (%)'
    ])

sample=np.array(I_stack_seg_filt)
VOL = sample!=0 #define an initial region coincident with the VOI to start the onion peeling quantification
z_total= sample.shape[0]
n_peels=21 #Define total the number of peels.

#create empty arrays to record the data
available_vol_count=[0]*n_peels
bone_vol_count=[0]*n_peels


for z in range(sample.shape[0]): #iterate in each Z

        print(f'{sample_name[:-21]}: {round((z/z_total)*100,1)}%') #plot the Z progress (%) in the terminal

        for p in range(n_peels): #successively erode-peel layers in each Z
            if p == 0:
                IE0 = VOL[z]
            else:
                IE0 = IE
            selem = disk(10) #peel size:10px (100um)
            IE = erosion(IE0,selem) #exterior circle of the disc
            IE1 = erosion(IE,selem) #interior circle of the disc
            IEband = np.logical_xor(IE,IE1) #determined the analysis band/peel as the difference between exterior and interior disc

            # quantify the bone in each band/peel:
            IR_band = sample[z]*IEband
            IR_available = (IR_band==1)|(IR_band==2) #volume availabe in each peel
            IR_bone = (IR_band==2) #bone volume in each peel

            #append the number of pixels of each group to the global radial peel count
            available_vol_count[p] += IR_available.sum()
            bone_vol_count[p] += IR_bone.sum()

#calculate the bone fraction for each peel
bone_fraction =list(map(truediv,bone_vol_count,available_vol_count)) #difives element by element in the arrays of equal length
bone_fraction = [0 if math.isnan(x) else x for x in bone_fraction] #in case of having NaN values result of 0/0, these are replaced by 0

#add the results to the data frame
for p in range(n_peels):
    df=df.append({
        'Sample family':sample_name[:1],
        'Sample name': sample_name[:-21],
        'peel layer': p,
        'Radium (mm)':(round(p*0.1,1)), #10px/peel*0.01mm/px=0.1mm/peel
        'available volume (px)': available_vol_count[p],
        'bone volume (px)': bone_vol_count[p],
        'bone faction (%)':bone_fraction[p]
        },ignore_index=True)
