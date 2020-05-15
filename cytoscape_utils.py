
#!/usr/bin/env python3

from py2cytoscape.data.cynetwork import CyNetwork
from py2cytoscape.data.cyrest_client import CyRestClient
from py2cytoscape.data.style import StyleUtil
import py2cytoscape.util.cytoscapejs as cyjs
import py2cytoscape.cytoscapejs as renderer
from py2cytoscape.cyrest.base import api
from py2cytoscape import cyrest
from IPython.display import Image
from IPython.display import SVG
from time import sleep


import networkx as nx
import pandas as pd
import json
import numpy as np


basic_settings = {

    # You can set default values as key-value pairs.


    'NODE_FILL_COLOR': '#f0f0f0',
    'NODE_SIZE': 55,
    'NODE_BORDER_WIDTH': 0,
    'NODE_LABEL_COLOR': '#555555',

    'EDGE_WIDTH': 2,
    #'EDGE_TRANSPARENCY': 100,
    #'EDGE_STROKE_UNSELECTED_PAINT': '#333333',


    'NETWORK_BACKGROUND_PAINT': 'white'

}


cytoscape=cyrest.cyclient()
cytoscape.version()



def get_visualizatio(sif, node_color,node_names_file, scores_file):

    ## create network
        cy = CyRestClient()
        cy.session.delete()
        # Step 2: Load network from somewhere
        ## entire path is needed
        net = cy.network.create_from(sif)

        dz = net.get_node_table()





        names = pd.read_csv(node_names_file)
        names = names[['GeneID', 'Symbol']]
        names['name'] = names['GeneID'].astype(str)

        dz = pd.merge(names, dz, on = 'name')


        # Now update existing node table with the data frame above.
        net.update_node_table(dz, network_key_col='name', data_key_col='name')



         # Step 6: Create Visual Style as code (or by hand if you prefer)
        my_style = cy.style.create('standard')
        my_style.update_defaults(basic_settings)






        ## apply this layout to get node degrees
        cy.layout.apply(name='degree-circle', network=net)
        degrees = net.get_node_column('degree.layout')
        cy.layout.apply(name='force-directed', network=net)
        # Step 5: Apply layout

        my_style.create_passthrough_mapping(column='Symbol', col_type='String',
                                            vp='NODE_LABEL')
        cy.style.apply(my_style, net)
