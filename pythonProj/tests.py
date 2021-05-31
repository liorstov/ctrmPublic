px.line(target.groupby(['x']).std(), y='semb', labels = dict(x = "time", semb  = "semb STD")).show()
px.scatter(target,'t', 'semb',color = 'x').show()
target = samples.loc[(samples.radius == 3) & (samples.y == 0) & (samples.x == 10) & (samples.t < 800)]
target = samples.loc[(samples.radius == 28) & (samples.y == 60)  & (samples.t < 800)]
melted = pd.melt(grouped, id_vars = ['location'], value_vars = ['mean', 'std' , 'max', 'min','median'] )
melted = pd.melt(grouped, ignore_index = False)
px.scatter(melted,x='location', y='value',color = 'variable').show()
px.scatter(grouped,x='mean', y='std').show()
go.Figure(go.Scatter(y=[1,2] ,customdata = [1,2], hovertemplate="%{customdata}<extra></extra>")).show()
px.imshow(lol.pivot(index = 'x', columns = 'y')).show()


np.array([range()])
for indexi,i in RGlocations.iterrows():
    for indexj,j in RGlocations.iterrows():
        xdist = i.X-j.X 
        ydist = i.Y-j.Y
        distmat[indexi,indexj] = np.sqrt( xdist**2 + ydist**2)
distmat[distmat==0]=np.Inf
distmat[:,0] = np.Inf
row = 0
for index in range(distmat.shape[0]-1):
    closest = np.argmin(distmat[row,:])
    print(closest)
    distmat[:,row] = np.Inf
    namelist[index+1]=RGlocations.iloc[closest][3]
    row = closest