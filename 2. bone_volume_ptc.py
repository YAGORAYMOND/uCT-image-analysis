#import libraries
import numpy as np

#quantification
seg= np.array(I_stack_seg_filt)

bkg = seg[seg==0].size       #Background (i.e. outside the VOI)
void = seg[seg==1].size      #Void (i.e. non mineralised tissue)
bone = seg[seg==2].size      #Bone
scaffold = seg[seg==3].size  #Scaffold

bv_av=(bone/(void+bone))*100 #Bone volume over avilable volume

print("Bone Volume / Available Volume = %d %%\n" %bv_av)
