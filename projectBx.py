# Projected components of B-field as arrays

import yt 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import matplotlib.image as mpimg

yt.enable_parallelism()

#Project B-field in a plane - derive these

def _Bxy(field, data): #define a new field in x-y plane
	return np.sqrt(data["X-magnfield"]**2+data["Y-magnfield"]**2)

def _Byz(field, data): #define a new field in y-z plane
	return np.sqrt(data["Y-magnfield"]**2+data["Z-magnfield"]**2)
	
def _Bxz(field, data): #define a new field in x-z plane
	return np.sqrt(data["X-magnfield"]**2+data["Z-magnfield"]**2)
    
fn="data.0510.3d.hdf5" 
ds=yt.load(fn) #load dataset

ds.add_field(("gas", "Bxy"), function=_Bxy, units="G", force_override=True) #add the new field
ds.add_field(("gas", "Byz"), function=_Byz, units="G", force_override=True)
ds.add_field(("gas", "Bxz"), function=_Bxz, units="G", force_override=True)

#Make projection plots of the components of the field. Turn these into 2D arrays

mppx=yt.ProjectionPlot(ds, "z", "X-magnfield", weight_field="density") #Project X-component of B-field from z-direction
#mppx.save("projected-Bx.png") #save
imgx=mpimg.imread("projected-Bx.png") #imread turns an image into a 2D array
#print imgx

mppy=yt.ProjectionPlot(ds, "z", "Y-magnfield", weight_field="density") #Project Y-component of B-field from z-direction
#mppy.save("projected-By.png")
imgy=mpimg.imread("projected-By.png")
#print imgy

#mppz=yt.ProjectionPlot(ds, "x", "Z-magnfield", weight_field="density")
#mppz.save("projected-Bz.png")
#imgz=mpimg.imread("projected-Bz.png")
#print imgz

#thetalist=[]

all_data_lvl0 = ds.smoothed_covering_grid(level=0, left_edge=[0,0.0,0.0], dims=ds.domain_dimensions)
X_mag_lvl0=all_data_lvl0["X-magnfield"]
Y_mag_lvl0=all_data_lvl0["Y-magnfield"]
Z_mag_lvl0=all_data_lvl0["Z-magnfield"]
flatx0=np.sum(X_mag_lvl0,axis=1)

all_data_lvl5 = ds.smoothed_covering_grid(level=5, left_edge=[0,0.0,0.0], dims=[2048,2048,2048])
X_mag_lvl5=all_data_lvl5["X-magnfield"]
Y_mag_lvl5=all_data_lvl5["Y-magnfield"]
Z_mag_lvl5=all_data_lvl5["Z-magnfield"]
flatx5=np.sum(X_mag_lvl5,axis=1)
print flatx5