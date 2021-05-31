import ctrm,time,os,numpy as np, pandas as pd
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from plotnine import *
from dplython import * 

#energy = np.fromfile("C:\\Users\\liors\\geoData\\for_lior\\memidion1\\LB_40",dtype = np.single).reshape(38,1000).transpose()
energyOrig = pd.DataFrame(np.fromfile("C:\\Users\\liors\\source\\repos\\grail exp\\red_grail_2020_10_22-23\\100106800.1.0.2020.10.23.05.47.00.000.Z.samp",dtype = np.single).reshape(44,60001).transpose())
target = (energyOrig.loc[58500:59500])
RGlocations = pd.read_csv("C:\\Users\\liors\\source\\repos\ctrm\\for_lior\\red_grail_locations.csv").sort_values(["X"],ignore_index=True)
RGlocations =RGlocations.loc[RGlocations.name.str.contains("GA")].reset_index()
RGlocations["Xfix"] = RGlocations.X-RGlocations.X.min()
RGlocations["Yfix"] = RGlocations.Y-RGlocations.Y.min()
RGlocations["Zfix"] = RGlocations.alt-RGlocations.alt.min()

a = ctrm.box(300,200,1,30,20,20,1000,0,150,30)
a.readVelo("C:\\Users\\liors\\source\\repos\\ctrm\\ctrm\\vrs_grail.txt")
#a.readEnergy("C:\\Users\\liors\\geoData\\for_lior\\memidion1\\LB_40")
a.setCoord(RGlocations[['Xfix','Yfix','Zfix']])
a.setEnergy(target)
a.getEnergy()

a.createImageSpace()
a.CalcSurfaceDist()
a.corrolationOnGeo()
pd.DataFrame(a.getEnergy()).groupby(0).count()
a.writeIP();
#a.writeSemblence();
b= a.getSample()

samples = pd.DataFrame({"IPindex": b[:,0],"x": b[:,1],"y":b[:,2], "t" :b[:,3], "radius":b[:,4],"semb":b[:,5]})
samples.to_pickle("grail.pkl")

#samples = samples[samples>0].dropna()
samples
sampleNumber = 0

grouped = samples.groupby(["x" ,"y","radius"]).mean()

plot(samples[samples.t==sampleNumber], sampleNumber)

ggplot(samples[samples.t==sampleNumber],aes(x="x", y="y"))+geom_point()