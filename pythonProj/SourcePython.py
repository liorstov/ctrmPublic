import ctrm,time,os,numpy as np, pandas as pd
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from plotnine import *
from dplython import * 

def plot(df,number): 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(df.x,df.y,df.radius,c = df.semb)
    colorb = plt.colorbar(scatter)
    colorb.ax.set_title("Semblance [-]")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Radius')
    ax.invert_zaxis()
    ax.set_title('CTRM\nsample number: %d' % number)
    
    plt.show()


a = ctrm.box(100,20,1,35,1,5)
a.readCoord("")
a.readVelo("")
a.readEnergy("")
a.createImageSpace()
a.CalcSurfaceDist()
a.corrolationOnGeo()
#a.writeIP();
#a.writeSemblence();
b= a.getSample()

samples = pd.DataFrame({"x": b[:,0],"y":b[:,1], "t" :b[:,2], "radius":b[:,3],"semb":b[:,4]})
#samples = samples[samples>0].dropna()
samples
sampleNumber = 0

grouped = samples.groupby(["x" ,"y","radius"]).mean()

plot(samples[samples.t==sampleNumber], sampleNumber)

ggplot(samples[samples.t==sampleNumber],aes(x="x", y="y"))+geom_point()