import ctrm,time,os,numpy as np, pandas as pd
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from plotnine import *
from dplython import * 


a = ctrm.box(80,25,1,35,1,1,1000,0)
a.readCoord("C:\\Users\\liors\\geoData\\for_lior\\memidion1\\aquisition38chan.txt")
a.readVelo("C:\\Users\\liors\\geoData\\for_lior\\memidion1\\vrs.txt")
a.readEnergy("C:\\Users\\liors\\geoData\\for_lior\\memidion1\\LB_40")
a.getEnergy()
a.createImageSpace()
a.CalcSurfaceDist()
a.corrolationOnGeo()
pd.DataFrame(a.getEnergy()).groupby(0).count()
a.writeIP();
#a.writeSemblence();
b= a.getSample()

samples = pd.DataFrame({"index": b[:,0],"x": b[:,1],"y":b[:,2], "t" :b[:,3], "radius":b[:,4],"semb":b[:,5]})
#samples = samples[samples>0].dropna()
samples
sampleNumber = 0

grouped = samples.groupby(["x" ,"y","radius"]).mean()

plot(samples[samples.t==sampleNumber], sampleNumber)

ggplot(samples[samples.t==sampleNumber],aes(x="x", y="y"))+geom_point()