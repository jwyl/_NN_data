# Projected components of B-field as arrays

import yt 
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import matplotlib.image as mpimg
from yt.visualization import fixed_resolution
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#yt.enable_parallelism()
import pickle

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

def Density2D(lvl,axis): #flatten the 3D covering grid of density to make a 2D grid. Specify level of resolution and axis as arguements
    density3d=Crop3Ddensity(lvl)
    density2d=np.sum(density3d,axis)
    return density2d
    
def pickledensity(lvl):
    xarray=Density2D(lvl, 0)
    xarray=np.asarray(xarray)
    fileObject1=open("yzdensitypickle%s" % lvl, 'wb') 
    pickle.dump(xarray,fileObject1)
    fileObject1.close()
    
    yarray=Density2D(lvl, 1)
    yarray=np.asarray(yarray)
    fileObject2=open("xzdensitypickle%s" % lvl, 'wb')
    pickle.dump(yarray,fileObject2)
    fileObject2.close()
    
    zarray=Density2D(lvl, 2)
    zarray=np.asarray(zarray)
    fileObject3=open("xydensitypickle%s" %lvl, 'wb')
    pickle.dump(zarray,fileObject3)
    fileObject3.close()
        
def Flattenx(lvl,axis):
    (Xmag3d,Ymag3d,Zmag3d)=Crop3DBfield(lvl)
    Bx=np.sum(Xmag3d,axis) #axis, 0 corresponds to x, 1 corresponds to y, 2 corresponds to z??
    return Bx
    
def Flatteny(lvl,axis):
    (Xmag3d,Ymag3d,Zmag3d)=Crop3DBfield(lvl)
    By=np.sum(Ymag3d,axis)
    return By
    
def Flattenz(lvl,axis):
    (Xmag3d,Ymag3d,Zmag3d)=Crop3DBfield(lvl)
    Bz=np.sum(Zmag3d,axis)
    return Bz

bg_color="YlGn"
quiv_color="autumn"
tick_locs=np.linspace(0,800,9)
tick_lbls=np.array(tick_locs*64e-4)

#from yt import YTArray

def normBxy(): #Create a normalised array of B vectors
    #densxy=Density2D(0,2) #integrated density along given axis
    U=Flattenx(0,2) #X-magnetic field integrated along given axis
    V=Flatteny(0,2) #Y-magnetic field
    U=np.asarray(U)#(zip(*x2)[::-1]) #rotate the matrix 90 degrees counter clockwise to correct orientation to match projected plots
    V=np.asarray(V)#(zip(*y2)[::-1])
    norm=np.sqrt(U**2+V**2) #magnitude of the vector
    Unorm=U/norm #normalise vectors 
    Vnorm=V/norm
    #mask_Unorm=np.ma.masked_where(densxy<np.mean(densxy),Unorm) #create a masked array of Unorm values only in high density regions
    #mask_Vnorm=np.ma.masked_where(densxy<np.mean(densxy),Vnorm)
    return Unorm, Vnorm, norm
    
def normBxz(): #Create a normalised array of B vectors
    #densxy=Density2D(0,2) #integrated density along given axis
    U=Flattenx(0,1) #X-magnetic field integrated along given axis
    V=Flattenz(0,1) #Y-magnetic field
    U=np.asarray(U)#(zip(*x2)[::-1]) #rotate the matrix 90 degrees counter clockwise to correct orientation to match projected plots
    V=np.asarray(V)#(zip(*y2)[::-1])
    norm=np.sqrt(U**2+V**2) #magnitude of the vector
    Unorm=U/norm #normalise vectors 
    Vnorm=V/norm
    #mask_Unorm=np.ma.masked_where(densxy<np.mean(densxy),Unorm) #create a masked array of Unorm values only in high density regions
    #mask_Vnorm=np.ma.masked_where(densxy<np.mean(densxy),Vnorm)
    return Unorm, Vnorm, norm
    
def normByz(): #Create a normalised array of B vectors
    #densxy=Density2D(0,2) #integrated density along given axis
    U=Flatteny(0,0) #X-magnetic field integrated along given axis
    V=Flattenz(0,0) #Y-magnetic field
    U=np.asarray(U)#(zip(*x2)[::-1]) #rotate the matrix 90 degrees counter clockwise to correct orientation to match projected plots
    V=np.asarray(V)#(zip(*y2)[::-1])
    norm=np.sqrt(U**2+V**2) #magnitude of the vector
    Unorm=U/norm #normalise vectors 
    Vnorm=V/norm
    #mask_Unorm=np.ma.masked_where(densxy<np.mean(densxy),Unorm) #create a masked array of Unorm values only in high density regions
    #mask_Vnorm=np.ma.masked_where(densxy<np.mean(densxy),Vnorm)
    return Unorm, Vnorm, norm
    
def pickleBxy():
    Unorm, Vnorm, norm=normBxy() #choose the z axis as direction 
    
    fileObject1=open("UnormBxypickle", 'wb') 
    pickle.dump(Unorm,fileObject1)
    fileObject1.close()
    
    fileObject2=open("VnormBxypickle", 'wb')
    pickle.dump(Vnorm,fileObject2)
    fileObject2.close()
    
    fileObject3=open("normBxypickle", 'wb')
    pickle.dump(norm,fileObject3)
    fileObject3.close()
    
def pickleBxz():
    Unorm, Vnorm, norm=normBxz() #choose the z axis as direction 
    
    fileObject1=open("UnormBxzpickle", 'wb') 
    pickle.dump(Unorm,fileObject1)
    fileObject1.close()
    
    fileObject2=open("VnormBxzpickle", 'wb')
    pickle.dump(Vnorm,fileObject2)
    fileObject2.close()
    
    fileObject3=open("normBxzpickle", 'wb')
    pickle.dump(norm,fileObject3)
    fileObject3.close()
    
def pickleByz():
    Unorm, Vnorm, norm=normByz() #choose the z axis as direction 
    
    fileObject1=open("UnormByzpickle", 'wb') 
    pickle.dump(Unorm,fileObject1)
    fileObject1.close()
    
    fileObject2=open("VnormByzpickle", 'wb')
    pickle.dump(Vnorm,fileObject2)
    fileObject2.close()
    
    fileObject3=open("normByzpickle", 'wb')
    pickle.dump(norm,fileObject3)
    fileObject3.close()
    
def pickleytdensity():
    ppz=yt.ProjectionPlot(ds, "z", "Bxy", weight_field="density") #Project X-component of B-field from z-direction
    Bz=ppz._frb["density"]
    Bz=np.asarray(Bz)
    file_name="testpickle"
    fileObject=open(file_name, 'wb')
    pickle.dump(Bz,fileObject)
    fileObject.close()
    
def quivBxy_dens():
    fig=plt.figure()
    ppz=yt.ProjectionPlot(ds, "z", "Bxy", weight_field="density") #Project X-component of B-field from z-direction
    Bz=ppz._frb["density"]
    Bz=np.asarray(Bz)
    ax=fig.add_subplot(111)
    plt.xticks(tick_locs,tick_lbls)
    plt.yticks(tick_locs,tick_lbls)
    Bzmag=ax.pcolormesh(np.log10(Bz),cmap=bg_color)
    cbar_m=plt.colorbar(Bzmag)
    cbar_m.set_label("density")
    res=800
    Unorm,Vnorm,norm=normBxy()
    X,Y=np.meshgrid(np.linspace(0,res,64, endpoint=True),np.linspace(0,res,64,endpoint=True))
    quivers=ax.quiver(X,Y,Unorm,Vnorm,norm*1e6,scale=50,cmap=quiv_color)
    cbar=plt.colorbar(quivers, orientation="horizontal")
    cbar.set_label('Bxy vectors (uG)')
    plt.title("Bxy on density projection")
    plt.xlabel("x (1e4 AU)")
    plt.ylabel("y (1e4 AU)")
    #plt.savefig("quiver_Bxy.png")

def quivByz_dens():
    fig=plt.figure()
    ppx=yt.ProjectionPlot(ds, "x", "Byz", weight_field="density") #Project X-component of B-field from z-direction
    Bx=ppx._frb["density"]
    ax=fig.add_subplot(111)
    plt.xticks(tick_locs,tick_lbls)
    plt.yticks(tick_locs,tick_lbls)
    Bxmag=ax.pcolormesh(np.log10(Bx), cmap=bg_color)
    cbar_m=plt.colorbar(Bxmag)
    cbar_m.set_label("density")
    res=800

    Unorm, Vnorm, norm=normBxy()
    X,Y=np.meshgrid(np.linspace(0,res,64,endpoint=True),np.linspace(0,res,64,endpoint=True))
    quivers=ax.quiver(X,Y,Unorm,Vnorm,norm*1e6,scale=50,cmap=quiv_color)
    cbar=plt.colorbar(quivers, orientation="horizontal")
    cbar.set_label('Byz vectors (uG)')
    plt.title("Byz vectors on weighted density Projection")
    plt.xlabel("(1e4 AU)")
    plt.ylabel("(1e4 AU)")
    #plt.savefig("quiversByz_colorsonbone.png")

def quivBxz_dens():
    fig=plt.figure()
    ppy=yt.ProjectionPlot(ds, "y", "Bxz", weight_field="density") #Project X-component of B-field from z-direction
    By=ppy._frb["density"]
    ax=fig.add_subplot(111)
    plt.xticks(tick_locs,tick_lbls)
    plt.yticks(tick_locs,tick_lbls)
    Bymag=ax.pcolormesh(np.log10(By), cmap=bg_color)
    cbar_m=plt.colorbar(Bymag)
    cbar_m.set_label("density")
    res=800

    Unorm, Vnorm, norm=normBxy()
    X,Y=np.meshgrid(np.linspace(0,res,64, endpoint=True),np.linspace(0,res,64,endpoint=True))
    quivers=ax.quiver(X,Y,Unorm,Vnorm,norm*1e6,scale=50,cmap=quiv_color)
    cbar=plt.colorbar(quivers, orientation="horizontal")
    cbar.set_label('Bxz vectors (uG)')
    plt.title("Bxz vectors on density projection")
    plt.xlabel("(1e4 AU)")
    plt.ylabel("(1e4 AU)")
    #plt.savefig("quiversByz_colorsonbone.png")
    
def plt3quiv(): #Plot quiver plots of all three planes
    quivBxy_dens()
    quivByz_dens()
    quivBxz_dens()

streamdensity=(2,2)

def streamlineBxy_dens():
    fig=plt.figure()
    ppz=yt.ProjectionPlot(ds, "z", "Bxy", weight_field="density") #Project X-component of B-field from z-direction
    Bz=ppz._frb["density"]
    ax=fig.add_subplot(111)
    tick_locs=np.linspace(0,800,9)
    tick_lbls=np.array(tick_locs*64e-4)
    plt.xticks(tick_locs,tick_lbls)
    plt.yticks(tick_locs,tick_lbls)
    Bzmag=ax.pcolormesh(np.log10(Bz),cmap="YlGn")
    cbar_m=plt.colorbar(Bzmag)
    cbar_m.set_label("density")
    res=800

    #densxy=Density2D(0,2) #integrated density along given axis
    U=Flattenx(0,2) #X-magnetic field integrated along given axis
    V=Flatteny(0,2) #Y-magnetic field
    #U=np.asarray(zip(*x2)[::-1]) #rotate the matrix 90 degrees to correct orientation to match projected plots
    #V=np.asarray(zip(*y2)[::-1])
    norm=np.sqrt(U**2+V**2) #magnitude of the vector
    Unorm=U/norm #normalise vectors 
    Vnorm=V/norm
    #mask_Unorm=np.ma.masked_where(densxy<np.mean(densxy),Unorm) #create a masked array of Unorm values only in high density regions
    #mask_Vnorm=np.ma.masked_where(densxy<np.mean(densxy),Vnorm)
    X,Y=np.meshgrid(np.linspace(0,res,64, endpoint=True),np.linspace(0,res,64,endpoint=True))
    streams=plt.streamplot(X,Y,Unorm,Vnorm,color=norm*1e6,density=streamdensity,cmap=plt.cm.autumn)
    cbar=plt.colorbar(orientation="horizontal")
    cbar.set_label('Bxy (uG)')
    plt.title("Bxy streamlines on density projection")
    plt.xlabel("(1e4 AU)")
    plt.ylabel("(1e4 AU)")
    
def streamlineByz_dens():
    fig=plt.figure()
    ppx=yt.ProjectionPlot(ds, "x", "Byz", weight_field="density") #Project X-component of B-field from z-direction
    Bx=ppx._frb["density"]
    ax=fig.add_subplot(111)
    plt.xticks(tick_locs,tick_lbls)
    plt.yticks(tick_locs,tick_lbls)
    Bxmag=ax.pcolormesh(np.log10(Bx), cmap="YlGn")
    cbar_m=plt.colorbar(Bxmag)
    cbar_m.set_label("density")
    res=800

    #densxy=Density2D(0,0) #integrated density along given axis
    U=Flatteny(0,0) #X-magnetic field integrated along given axis
    V=Flattenz(0,0) #Y-magnetic field
    #U=np.asarray(zip(*x2)[::-1]) #rotate the matrix 90 degrees to correct orientation to match projected plots
    #V=np.asarray(zip(*y2)[::-1])
    norm=np.sqrt(U**2+V**2) #magnitude of the vector
    Unorm=U/norm #normalise vectors 
    Vnorm=V/norm
    #mask_Unorm=np.ma.masked_where(densxy<np.mean(densxy),Unorm) #create a masked array of Unorm values only in high density regions
    #mask_Vnorm=np.ma.masked_where(densxy<np.mean(densxy),Vnorm)
    X,Y=np.meshgrid(np.linspace(0,res,64, endpoint=True),np.linspace(0,res,64,endpoint=True))
    streams=plt.streamplot(X,Y,Unorm,Vnorm,color=norm*1e6,density=(3,3),cmap=plt.cm.autumn)
    cbar=plt.colorbar(orientation="horizontal")
    cbar.set_label('Byz (uG)')
    plt.title("Byz streamlines on weighted density projection")
    plt.xlabel("(1e4 AU)")
    plt.ylabel("(1e4 AU)")
    #plt.savefig("quiversByz_colorsonbone.png")
    
def streamlineBxz_dens():
    fig=plt.figure()
    ppy=yt.ProjectionPlot(ds, "y", "Bxz", weight_field="density") #Project X-component of B-field from z-direction
    By=ppy._frb["density"]
    ax=fig.add_subplot(111)
    plt.xticks(tick_locs,tick_lbls)
    plt.yticks(tick_locs,tick_lbls)
    Bymag=ax.pcolormesh(np.log10(By), cmap="YlGn")
    cbar_m=plt.colorbar(Bymag)
    cbar_m.set_label("density")
    res=800

    #densxy=Density2D(0,1) #integrated density along given axis
    U=Flattenx(0,1) #X-magnetic field integrated along given axis
    V=Flattenz(0,1) #Z-magnetic field
    #U=np.asarray(zip(*x2)[::-1]) #rotate the matrix 90 degrees to correct orientation to match projected plots
    #V=np.asarray(zip(*y2)[::-1])
    norm=np.sqrt(U**2+V**2) #magnitude of the vector
    Unorm=U/norm #normalise vectors 
    Vnorm=V/norm
    #mask_Unorm=np.ma.masked_where(densxy<np.mean(densxy),Unorm) #create a masked array of Unorm values only in high density regions
    #mask_Vnorm=np.ma.masked_where(densxy<np.mean(densxy),Vnorm)
    X,Y=np.meshgrid(np.linspace(0,res,64, endpoint=True),np.linspace(0,res,64,endpoint=True))
    streams=plt.streamplot(X,Y,Unorm,Vnorm,color=norm*1e6,density=(3,3),cmap=plt.cm.autumn)
    cbar=plt.colorbar(orientation="horizontal")
    cbar.set_label('Bxz streamlines (uG)')
    plt.title("Bxz streamlines on weighted density projection")
    plt.xlabel("(1e4 AU)")
    plt.ylabel("(1e4 AU)")
    #plt.savefig("quiversByz_colorsonbone.png")

def plt3stream():
    streamlineBxy_dens()
    streamlineByz_dens()
    streamlineBxz_dens()
    
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
def interact3D(): #create a 3D projection of the magnetic field vectors
    X=Flattenx(0,2)
    Y=Flatteny(0,2)
    Z=Flattenz(0,2)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    logX=np.log10(X)
    logY=np.log10(Y)
    logZ=np.log10(Z)
    ax.plot_surface(X,Y,Z,rstride=1,cstride=1, alpha=0.3)
    cset = ax.contour(X,Y,Z, zdir='z', offset=0.)
    #cset = ax.contour(logX,logY,logZ, zdir='x', offset=0.)
    #cset = ax.contour(logX,logY,logZ, zdir='y', offset=0.)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()


def pltytfield():
    fig=plt.figure()
    ppy=yt.ProjectionPlot(ds, "x", "Bxy", weight_field="density") #Project X-component of B-field from z-direction
    By=ppy._frb["Bxy"]
    ax=fig.add_subplot(111)
    plt.xticks(tick_locs,tick_lbls)
    plt.yticks(tick_locs,tick_lbls)
    Bymag=ax.pcolormesh(np.log10(By))
    cbar_m=plt.colorbar(Bymag)
    cbar_m.set_label("Bxy")
    plt.title("Bxy in yz plane")

    fig=plt.figure()
    ppy=yt.ProjectionPlot(ds, "y", "Bxy", weight_field="density") #Project X-component of B-field from z-direction
    By=ppy._frb["Bxy"]
    ax=fig.add_subplot(111)
    plt.xticks(tick_locs,tick_lbls)
    plt.yticks(tick_locs,tick_lbls)
    Bymag=ax.pcolormesh(np.log10(By))
    cbar_m=plt.colorbar(Bymag)
    cbar_m.set_label("Bxy")
    plt.title("Bxy in xz plane")

    fig=plt.figure()
    ppy=yt.ProjectionPlot(ds, "z", "Bxy", weight_field="density") #Project X-component of B-field from z-direction
    By=ppy._frb["Bxy"]
    ax=fig.add_subplot(111)
    plt.xticks(tick_locs,tick_lbls)
    plt.yticks(tick_locs,tick_lbls)
    Bymag=ax.pcolormesh(np.log10(By))
    cbar_m=plt.colorbar(Bymag)
    cbar_m.set_label("Bxy")
    plt.title("Bxy in xy plane")