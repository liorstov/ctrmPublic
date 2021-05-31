
RGlocations = pd.read_csv("C:\\Users\\liors\\source\\repos\ctrm\\for_lior\\red_grail_locations.csv").sort_values(["X"],ignore_index=True)
distmat = np.ones([RGlocations.shape[0],RGlocations.shape[0]])*100000
namelist = ['']*RGlocations.shape[0]
namelist[0] = RGlocations.iloc[0][3]

RGlocations["Xfix"] = RGlocations.X-RGlocations.X.min()
RGlocations["Yfix"] = RGlocations.Y-RGlocations.Y.min()
RGlocations["Zfix"] = RGlocations.alt-RGlocations.alt.min()

cols = np.linspace(0,RGlocations["Xfix"].max(), num=250)
rows =  np.linspace(0,RGlocations["Yfix"].max(), num=250)
RGlocations['col'] = np.searchsorted(cols, RGlocations['Xfix'])
RGlocations['row'] = np.searchsorted(rows, RGlocations['Yfix'])
(RGlocations.groupby(['col','row']).count()['name'] >1).any()

energyOrig = pd.DataFrame(np.fromfile("C:\\Users\\liors\\source\\repos\\grail exp\\red_grail_2020_10_22-23\\100106800.1.0.2020.10.23.05.47.00.000.Z.samp",dtype = np.single).reshape(44,60001).transpose())
energy = energyOrig.groupby(np.arange(len(energyOrig))//1000)
target = (energyOrig.loc[58500:59500])
px.imshow(energy.max(),aspect='auto').show()

RGlocations["data"] = energy.loc[25]
px.scatter(RGlocations, x = "col", y = "row", color = 'data',hover_name = 'name', size = np.ones(49)*5).show()

target = (energy.apply(lambda x: x.abs().max())/energy.mean().abs())

test  = RGlocations.join(energyOrig.transpose())
test2 =pd.melt(test,RGlocations.columns)
px.scatter(test2, x = "col", y = "row", color = 'value',hover_name = 'name', size = np.ones(2940049)*5,animation_frame="variable").show()
