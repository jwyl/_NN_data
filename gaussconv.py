#Unpickle arrays from projections and reproject them
import matplotlib.pyplot as plt
import scipy as sp
from scipy import ndimage
import pickle
import numpy as np

file_name="testpickle" #test pickle is the array of the density projection
fileObject=open(file_name,'r')
densityarr=pickle.load(fileObject)

file_name1="UnormBxypickle" #the 64 by 64 array of x-component normalised
fileObject1=open(file_name1, 'r')
Unorm=pickle.load(fileObject1)

fileObject2=open("VnormBxypickle", 'r') #the 64 by 64 array of the y-componenet normalised
Vnorm=pickle.load(fileObject2)

fileObject3=open("normBxypickle", 'r') #an array of the magnitude of the vectors
norm=pickle.load(fileObject3)

tick_locs=np.linspace(0,800,9)
tick_lbls=np.array(tick_locs*64e-4)

def gaussconv(sigma): #a function to convolve the arrays of the x and y components
    Ugauss=sp.ndimage.gaussian_filter(Unorm,sigma,mode="wrap")
    Vgauss=sp.ndimage.gaussian_filter(Vnorm,sigma,mode="wrap")
    return Ugauss, Vgauss
    
def plotgaussquiv(sigma):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    plt.xticks(tick_locs,tick_lbls)
    plt.yticks(tick_locs,tick_lbls)
    densitymap=ax.pcolormesh(np.log10(densityarr),cmap='YlGn')
    cbar_m=plt.colorbar(densityarr)
    cbar_m.set_label("density")
    res=800
    Ugauss,Vgauss=gaussconv(sigma)
    X,Y=np.meshgrid(np.linspace(0,res,64, endpoint=True),np.linspace(0,res,64,endpoint=True))
    quivers=ax.quiver(X,Y,Ugauss,Vgauss,norm*1e6,scale=50,cmap='autumn')
    cbar=plt.colorbar(quivers, orientation="horizontal")
    cbar.set_label('Bxy (uG)')
    plt.title("Gaussian Convolved Bxy on density projection")
    plt.xlabel("x (1e4 AU)")
    plt.ylabel("y (1e4 AU)")
    
def plotgaussstream(sigma):
    fig=plt.figure()
    #ax=fig.add_subplot(111)
    plt.xticks(tick_locs,tick_lbls)
    plt.yticks(tick_locs,tick_lbls)
    Bzmag=plt.pcolormesh(np.log10(Bz),cmap='YlGn')
    cbar_m=fig.colorbar(Bzmag)
    cbar_m.set_label("density")
    res=800
    Ugauss,Vgauss=gaussconv(sigma)
    X,Y=np.meshgrid(np.linspace(0,res,64, endpoint=True),np.linspace(0,res,64,endpoint=True))
    streams=plt.streamplot(X,Y,Ugauss,Vgauss,color=norm*1e6,density=(2,2),cmap='autumn')
    cbar=plt.colorbar(orientation='horizontal')
    cbar.set_label('Bxy (uG)')
    plt.title("Gaussian Convolved Bxy on density projection")
    plt.xlabel("x (1e4 AU)")
    plt.ylabel("y (1e4 AU)")