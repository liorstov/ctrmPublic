px.line(target.groupby(['x']).std(), y='semb', labels = dict(x = "time", semb  = "semb STD")).show()
px.scatter(target,'t', 'semb',color = 'x').show()
target = samples.loc[(samples.radius == 25) & (samples.y == 60) & (samples.x == 10) & (samples.t < 800)]
target = samples.loc[(samples.radius == 28) & (samples.y == 60)  & (samples.t < 800)]
melted = pd.melt(grouped, id_vars = ['location'], value_vars = ['mean', 'std' , 'max', 'min','median'] )
melted = pd.melt(grouped, ignore_index = False)
px.scatter(melted,x='location', y='value',color = 'variable').show()
px.scatter(grouped,x='mean', y='std').show()
go.Figure(go.Box(y=target.semb, x = target.get_level_values(2), boxmean='sd')).show()
