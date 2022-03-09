import plotly.express as px
import plotly.graph_objs as go

import pandas as pd
import numpy as np

from collections import Counter
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import plot, iplot, init_notebook_mode
import plotly.graph_objs as go
init_notebook_mode(connected=True)



somali = pd.read_csv('somalier-ancestry.tsv', sep='\t')
somali = somali[somali['#sample_id'].str.startswith("P00")]
somali.loc[somali['predicted_ancestry'] == 'AFR', 'predicted_ancestry'] = 'African'
somali.loc[somali['predicted_ancestry'] == 'AMR', 'predicted_ancestry'] = 'Ad Mixed American'
somali.loc[somali['predicted_ancestry'] == 'EAS', 'predicted_ancestry'] = 'East Asian'
somali.loc[somali['predicted_ancestry'] == 'EUR', 'predicted_ancestry'] = 'European'
somali.loc[somali['predicted_ancestry'] == 'SAS', 'predicted_ancestry'] = 'South Asian'

PCs = somali[[ 'PC1', 'PC2', 'PC3','PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']]

variances = np.var(PCs, axis = 0)

PC3= PC3.rename(columns={'PC1': 'PC2(38.9%)', 'PC2': 'PC1(45.8%)','PC4':'PC3(9.1%)'})

fig = px.scatter_3d(
    PC3, x='PC1(45.8%)', y='PC2(38.9%)', z='PC3(9.1%)', color=somali['predicted_ancestry'],
    color_discrete_sequence=['#4C78A8' ,'#54A24B' ,'#F58518',  '#E45756', '#EECA3B'])
fig.update_traces(marker_size=3)
fig.update_layout(
    autosize=True,
    width=1200,
    height=1000,
    font=dict(
        family='Times New Roman'),
   plot_bgcolor='rgba(0,0,0,0)',
    scene=dict(
        xaxis=dict(showticklabels=False,showline=True,linecolor = '#636363',
      linewidth = 2),
        yaxis=dict(showticklabels=False,showline=True,linecolor = '#636363',
      linewidth = 2),
        zaxis=dict(showticklabels=False,showline=True,linecolor = '#636363',
      linewidth = 2),
    )
)


fig.write_image("fig3d.pdf",engine="kaleido")


fig = px.bar(counts, x = 'count', y = 'z ' ,orientation='h', color='predicted_ancestry',color_discrete_sequence= ['#EECA3B','#E45756','#F58518','#4C78A8','#54A24B' ])
fig.update_layout(
    height=300,
    font=dict(
        family='Times New Roman'),
    xaxis =dict(showgrid = False,showticklabels=False),
    yaxis =dict( showgrid = False,showticklabels=False),
    plot_bgcolor='rgba(0,0,0,0)'
)

fig.show()
fig.write_image("fig2bar.pdf",engine="kaleido")