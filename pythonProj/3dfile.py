import plotly.express as px
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output,State
import plotly.graph_objects as go , plotly.io as pio# or plotly.express as px
import pandas as pd
import numpy as np
from mpl_toolkits.mplot3d import axes3d

pio.templates.default = 'plotly_white' 
#samples = pd.read_pickle("c:\\users\\liors\\source\\repos\\ctrm\\ctrm\\sample.pkl")
samples = pd.read_pickle("c:\\users\\liors\\source\\repos\\ctrm\\ctrm\\grail.pkl")
RGlocations = pd.read_csv("C:\\Users\\liors\\source\\repos\ctrm\\for_lior\\red_grail_locations.csv").sort_values(["X"],ignore_index=True)
RGlocations =RGlocations.loc[RGlocations.name.str.contains("GA")].reset_index()
RGlocations["Xfix"] = RGlocations.X-RGlocations.X.min()
RGlocations["Yfix"] = RGlocations.Y-RGlocations.Y.min()
RGlocations["Zfix"] = RGlocations.alt-RGlocations.alt.min()


grouped = samples.groupby(["IPindex","x" ,"y","radius"]).agg({"semb": ['mean','min','max','median', 'std' ]}).droplevel(axis=1, level=0).reset_index()

#data = [go.Scatter3d(x = grouped.x,y=grouped.y,z=grouped.radius, marker=dict(size=12, color=grouped.semb, colorscale='Reds', opacity=0.5))]
app = dash.Dash(__name__, external_stylesheets=["https://codepen.io/chriddyp/pen/bWLwgP.css"])


layout = go.Layout(
    #scene=go.layout.Scene(
    #    xaxis=go.layout.scene.XAxis(
    #        showspikes=True,
    #        spikecolor='#1fe5bd',
    #        spikethickness=10,
    #    ),
    #    yaxis=go.layout.scene.YAxis(
    #        showspikes=False,
    #        spikecolor='#1fe5bd',
    #        spikethickness=6,range=[0,200],
    #    ),
    #    zaxis=go.layout.scene.ZAxis(
    #        showspikes=False,
    #        spikecolor='black',
    #        spikethickness=10,
    #    ),
    #),
    height = 700,width = 700
)
#grouped.loc[(grouped.x==96 )&( grouped.y==140),'semb'] = 1
#fig = go.Figure(data =data, layout=layout)
#fig = px.scatter_3d(grouped, x = 'x', y='y',z='radius', color = 'semb',color_continuous_scale = 'Reds', opacity = 0.3)
#fig2d = px.scatter()
#fig.update_scenes(zaxis_autorange="reversed",yaxis_autorange="reversed",xaxis_autorange="reversed")
colorsc = 0;
theme =  {
    'dark': True,
    'detail': '#007439',
    'primary': '#00EA64',
    'secondary': '#6E6E6E',
}

app = dash.Dash()
app.layout = html.Div([
    html.Div([
        html.Div(dcc.Graph(id = "2dplot"),style={'display': 'inline-block'}),
        html.Div(dcc.Graph(id = "3dplot",figure = px.scatter_3d()),style={'display': 'inline-block'}),    
        html.Div(dcc.Slider(id = "bar", value = 0,max = samples.semb.max(), step = samples.semb.max()/10,min = 0,updatemode='drag',vertical=True),style={'display': 'inline-block'}),
        html.Div(dcc.Dropdown(id = "aggregation-dropdown",options=[{'label': 'mean', 'value': 'mean'},  {'label': 'maximum', 'value': 'max'},{'label': 'std', 'value': 'std'},{'label': 'median', 'value': 'median'}], value='mean', clearable=False),style={'display': 'inline-block', 'width': "5%"}),
        html.Button('Submit', id='submit-val', n_clicks=0)
        ]),
    
    dcc.Input(id="avrWindow", type = "number", placeholder="window size",min = 1,max = 10,value = 1),
    
    dcc.Graph(id = "ImagePointView"),
    dcc.Store(id = "radius")])
    
    #html.Div([dcc.Slider(id='radius',min=1,max=35,value=1,marks={str(rad): str(rad) for rad in range(1,35)},step=None,vertical=True)]
    #         ,style={'display': 'inline-block'})],style={'display': 'block','height':'100%'})

@app.callback(
    Output("ImagePointView", "figure"),
    Input('3dplot', 'clickData'),
    Input('2dplot', 'clickData'))
def updateIP(selectedPoint,selectedPoint2):
    callbackVal = dash.callback_context.triggered[0]
    print(callbackVal)
    x=y=radius =radius= 0;
    if callbackVal['value'] is not None:   
        x,y = callbackVal['value']['points'][0]['x'],callbackVal['value']['points'][0]['y']
        if ((callbackVal["prop_id"]== '3dplot.clickData') ):
            radius = callbackVal['value']['points'][0]['customdata'][0]
        else:
            radius = callbackVal['value']['points'][0]['customdata']

        print(x,y,radius)
    target = samples[(samples.x == x)&(samples.y == y)&(samples.radius==radius)]
    title = "semblance for x: %d, y: %d, depth: %d" %(x,y,radius)
    fig = px.line(target,'t', 'semb',title = title)
    return(fig)

@app.callback(
    Output('2dplot', 'figure'),
    Output('radius', 'data'),
    Input('3dplot', 'clickData'),
    Input('3dplot',"figure"),
    State("aggregation-dropdown","value"))    
def update_figure(selected_rad,dfig,aggregation):
    colorsc = dfig['layout']['coloraxis']['colorscale']
    radius = 0;
    if selected_rad is not None:     
        radius = selected_rad['points'][0]['z']
    filtered = grouped.loc[grouped.radius==radius];
    fig2d = go.Figure(data=go.Heatmap(
        z=filtered[aggregation],
        x=filtered["x"],
        y=filtered["y"],
        zmin = min(grouped[aggregation]),
        zmax = max(grouped[aggregation]),
        colorscale=colorsc,zsmooth = 'best', customdata = np.ones(filtered.size)*radius))
    fig2d.add_trace(go.Scatter(x=RGlocations["Xfix"], y=RGlocations["Yfix"],mode = 'markers', marker = dict(color = "Blue")))
    #fig2d = px.imshow(filtered.pivot_table(index = 'y', columns = 'x',values = aggregation),color_continuous_scale =colorsc,origin = "lower",aspect = "auto", labels = {'customdata':np.ones(filtered.size)*radius})
    fig2d.update_layout(title = ("radius: %d" %radius), xaxis = {"title":"x"},yaxis = {"title":"y"})
    title="Plot Title"
    #fig.add_mesh3d(x=[0,0,int(max(grouped.x)),int(max(grouped.x))],y=[0,int(max(grouped.y)),int(max(grouped.y)),0],z=[selected_rad,selected_rad,selected_rad,selected_rad])
    return(fig2d,{'radius':radius})
@app.callback(
    Output("3dplot", "figure"),
    Output("bar", "max"),
    Output("bar", "step"),
    Input("submit-val","n_clicks"),
    Input("bar","value"),
    State("aggregation-dropdown","value"),
    State("avrWindow","value"))
def update_3dplot(clicks,bar, selectedAgg,window):   

    title = "Results are aggregated with: %s\n semblance > %s " %(selectedAgg, bar)

    fig = px.scatter_3d(grouped, x = 'x', y='y',z='radius', color = selectedAgg,color_continuous_scale = 'Reds', opacity = 0.3, title = title ,template="plotly_white",size = (grouped[selectedAgg].ravel()>bar)*1, hover_data = {'customdata':grouped.radius})
    #fig = go.Figure(go.Scatter3d(x = grouped.x, y=grouped.y,z=grouped.radius ,marker=dict(size = (grouped[selectedAgg].ravel()>bar)*5,color = grouped[selectedAgg], colorscale= "Reds"),mode = "markers"))
    fig.update_scenes(zaxis_autorange="reversed",yaxis_autorange="reversed",xaxis_autorange="reversed")
    fig.update_traces(marker = dict(line=dict(width=0)))
    return (fig, grouped[selectedAgg].max(),grouped[selectedAgg].max()/10 )


app.run_server(debug=True, use_reloader=True)  # Turn off reloader if inside Jupyter
fig = px.scatter_3d(grouped, x = 'x', y='y',z='radius', color = selectedAgg,color_continuous_scale = 'Reds', opacity = 0.3, title = title,template="plotly_white",size = (grouped[selectedAgg].ravel()>0.5)*1)
fig.update_scenes(zaxis_autorange="reversed",yaxis_autorange="reversed",xaxis_autorange="reversed")
fig.update_xaxis(range = [0,2])
fig.update_layout(xaxis=dict(rangeslider=dict(visible=True)))
fig.show()