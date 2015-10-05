#plotting from arrays
import matplotlib.pyplot as plt
import scipy as sp
from scipy import ndimage
import pickle
import numpy as np

def unpickleBxy():
    file_name="xydensitypickle3" #test pickle is the array of the density projection
    fileObject=open(file_name,'r')
    density_xy=pickle.load(fileObject)
    fileObject.close()

    file_name1="UnormBxypickle" #the 64 by 64 array of x-component normalised
    fileObject1=open(file_name1, 'r')
    Unorm=pickle.load(fileObject1)
    fileObject1.close()

    fileObject2=open("VnormBxypickle", 'r') #the 64 by 64 array of the y-componenet normalised
    Vnorm=pickle.load(fileObject2)
    fileObject2.close()

    fileObject3=open("normBxypickle", 'r') #an array of the magnitude of the vectors
    norm=pickle.load(fileObject3)
    fileObject3.close()
    return density_xy, Unorm, Vnorm, norm
    
tick_locs=np.linspace(0,800,9)
tick_lbls=np.array(tick_locs*64e-4)
bg_color="YlGn"
quiv_color="autumn"

def pltdensityarray(array): # the background and plane of projection
    fig=plt.figure()
    ax=fig.add_subplot(111)
    plt.xticks(tick_locs,tick_lbls)
    plt.yticks(tick_locs,tick_lbls)
    densitycplt=ax.pcolormesh(np.log10(array),cmap=bg_color)
    cbar_m=plt.colorbar(densitycplt)
    cbar_m.set_label("density")
    
def quivBxy_dens():
    density_xy, Unorm, Vnorm, norm = unpickleBxy()
    pltdensityarray(density_xy) #either xydensity, xzdensity or yzdensity
    res=512
    X,Y=np.meshgrid(np.linspace(0,res,64, endpoint=True),np.linspace(0,res,64,endpoint=True))
    quivers=plt.quiver(X,Y,Unorm,Vnorm,norm*1e6,scale=50,cmap=quiv_color)
    cbar=plt.colorbar(quivers, orientation="horizontal")
    cbar.set_label('Bxy vectors (uG)')
    plt.title("Bxy on density projection")
    plt.xlabel("x (1e4 AU)")
    plt.ylabel("y (1e4 AU)")
    #plt.savefig("quiver_Bxy.png")

def quivByz_dens():
    pltdens3array(bg_orientaion)
    
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

    Unorm, Vnorm, norm=normBxy(1)
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
    