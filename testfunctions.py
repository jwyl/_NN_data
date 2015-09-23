# Projected components of B-field as arrays

import yt 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import matplotlib.image as mpimg
from yt.visualization import fixed_resolution
from mpl_toolkits.axes_grid1 import make_axes_locatable
#yt.enable_parallelism()

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

all_data_lvl3 = ds.smoothed_covering_grid(level=3, left_edge=[0.0,0.0,0.0], dims=[512,512,512])
X_mag_lvl3=all_data_lvl3["X-magnfield"]
#Y_mag_lvl3=all_data_lvl3["Y-magnfield"]
#Z_mag_lvl3=all_data_lvl3["Z-magnfield"]
#flatx31=np.sum(X_mag_lvl3,axis=1)
#flatx32=np.sum(X_mag_lvl3,axis=2)

def Crop3DBfield(lvl): #make a 3D covering grid of all 3 of the B-field components, a function of the AMRgrid level (resolution)
    dims=64*2**lvl
    crop_data_lvl= ds.smoothed_covering_grid(level=lvl, left_edge=[-4.0e17,-4.0e17,-4.0e17], dims=[dims,dims,dims])
    Xmag3d=crop_data_lvl["X-magnfield"]
    Ymag3d=crop_data_lvl["Y-magnfield"]
    Zmag3d=crop_data_lvl["Z-magnfield"]
    return Xmag3d,Ymag3d,Zmag3d

def Crop3Ddensity(lvl): #make a 3D covering grid of the density
    dims=64*2**lvl
    crop_data_lvl=ds.smoothed_covering_grid(level=lvl, left_edge=[-4.0e17,-4.0e17,-4.e17], dims=[dims,dims,dims])
    density3d=crop_data_lvl["density"]
    return density3d

def Density2D(lvl,axis):
    density3d=Crop3Ddensity(lvl)
    density2d=np.sum(density3d,axis)
    return density2d
        
def Flattenx(lvl,axis):
    (Xmag3d,Ymag3d,Zmag3d)=Crop3DBfield(lvl)
    Bx=np.sum(Xmag3d,axis)
    return Bx
    
def Flatteny(lvl,axis):
    (Xmag3d,Ymag3d,Zmag3d)=Crop3DBfield(lvl)
    By=np.sum(Ymag3d,axis)
    return By
    
def Flattenz(lvl,axis):
    (Xmag3d,Ymag3d,Zmag3d)=CropThreeDarray(lvl)
    Bz=np.sum(Zmag3d,axis)
    return Bz

fig=plt.figure()
ppx=yt.ProjectionPlot(ds, "z", "Bxy", weight_field="density") #Project X-component of B-field from z-direction
B=ppx._frb["Bxy"]
ax=fig.add_subplot(111)
tick_locs=np.linspace(0,800,9)
tick_lbls=np.array(tick_locs*64e-4)
plt.xticks(tick_locs,tick_lbls)
plt.yticks(tick_locs,tick_lbls)
mag=ax.pcolormesh(np.log10(B),cmap="bone")
cbar_m=plt.colorbar(mag)
cbar_m.set_label("Bxy projected field strength")
res=800

densxy=Density2D(0,2) #integrated density along given axis
x2=Flattenx(0,2) #X-magnetic field integrated along given axis
y2=Flatteny(0,2) #Y-magnetic field
U=np.asarray(zip(*x2)[::-1]) #rotate the matrix 90 degrees to correct orientation to match projected plots
V=np.asarray(zip(*y2)[::-1])
norm=np.sqrt(U**2+V**2) #magnitude of the vector
Unorm=U/norm #normalise vectors 
Vnorm=V/norm
mask_Unorm=np.ma.masked_where(densxy<np.mean(densxy),Unorm) #create a masked array of Unorm values only in high density regions
mask_Vnorm=np.ma.masked_where(densxy<np.mean(densxy),Vnorm)
X,Y=np.meshgrid(np.linspace(0,res,64, endpoint=True),np.linspace(0,res,64,endpoint=True))
quivers=ax.quiver(X,Y,Unorm,Vnorm,norm*1e6,scale=50)
cbar=plt.colorbar(quivers, orientation="horizontal")
cbar.set_label('Bxy vectors (uG)')
plt.title("Bxy")
plt.xlabel("x (1e4 AU)")
plt.ylabel("y (1e4 AU)")
#plt.savefig("quiver_Bxy.png")


fig=plt.figure()
ppx=yt.ProjectionPlot(ds, "z", "Byz", weight_field="density") #Project X-component of B-field from z-direction
B=ppx._frb["Byz"]
ax=fig.add_subplot(111)
tick_locs=np.linspace(0,800,9)
tick_lbls=np.array(tick_locs*64e-4)
plt.xticks(tick_locs,tick_lbls)
plt.yticks(tick_locs,tick_lbls)
mag=ax.pcolormesh(np.log10(B))
#divider = make_axes_locatable(ax)
#ax_1 = divider.append_axes("right", size="5%", pad=0.05)
cbar_m=plt.colorbar(mag, fraction=0.046, pad=0.04)
cbar_m.set_label("Byz")
res=800

densxy=Density2D(0,0) #integrated density along given axis
x2=Flattenx(0,0) #X-magnetic field integrated along given axis
y2=Flatteny(0,0) #Y-magnetic field
U=np.asarray(zip(*x2)[::-1]) #rotate the matrix 90 degrees to correct orientation to match projected plots
V=np.asarray(zip(*y2)[::-1])
norm=np.sqrt(U**2+V**2) #magnitude of the vector
Unorm=U/norm #normalise vectors 
Vnorm=V/norm
mask_Unorm=np.ma.masked_where(densxy<np.mean(densxy),Unorm) #create a masked array of Unorm values only in high density regions
mask_Vnorm=np.ma.masked_where(densxy<np.mean(densxy),Vnorm)
X,Y=np.meshgrid(np.linspace(0,res,64, endpoint=True),np.linspace(0,res,64,endpoint=True))
quivers=ax.quiver(X,Y,Unorm,Vnorm,norm*1e6,scale=50,cmap="bone")
#cax_2 = divider.append_axes("top", size="5%", pad=0.05)
cbar=plt.colorbar(quivers, orientation="horizontal")
cbar.set_label('Byz vectors (uG)')
plt.title("Byz vectors on weighted Byz Projection")
plt.xlabel("x (1e4 AU)")
plt.ylabel("y (1e4 AU)")
#plt.savefig("quiversByz_colorsonbone.png")