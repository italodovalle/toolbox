
#!/usr/bin/python

from py2cytoscape.data.cynetwork import CyNetwork
from py2cytoscape.data.cyrest_client import CyRestClient
from py2cytoscape.data.style import StyleUtil
import py2cytoscape.util.cytoscapejs as cyjs
import py2cytoscape.cytoscapejs as renderer
from py2cytoscape.cyrest.base import api
from py2cytoscape import cyrest
#from IPython.display import Image

import networkx as nx
import pandas as pd
import json
import numpy as np

import argparse







def create_network(infile):
    cy = CyRestClient()


    # Reset
    cy.session.delete()

    # Step 2: Load network from somewhere
    ## entire path is needed
    net = cy.network.create_from(infile)

    return(net)




if __name__ == '__main__':

    ### Interactome params

    description = 'Run Proximity'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-e', required=True, dest='edge_list',
                        action='store',help='edge list file')
    parser.add_argument('-t', required=True, dest='tmp_file',action='store',
                        help='Name tmp sif file', default = 'net.sif')
    parser.add_argument('-i', required=True, dest='index_col',action='store',
                        help='index col [True, False]', default = True)
    parser.add_argument('-n', dest='node_attributes', action='store',
                        help='file with node labels')
    parser.add_argument('-k', dest='key_col', action='store',
                        help='node key')
    parser.add_argument('-l', dest='label_col', action='store',
                        help='label col')
    parser.add_argument('-s', dest='node_scores', action = 'store',
                        help = 'node scores [True, False]', default=False)
    ## parse node scores file
    parser.add_argument('-o', required=True, dest='outfile', action='store',
                        help='outfile', default = 'output.pdf')

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()

    edge_list = args.edge_list
    tmp_file = args.tmp_file
    index_col = args.index_col
    key_col = args.key_col
    label_col = args.label_col
    node_attributes = args.node_attributes
    node_scores = args.node_scores
    out_figure = args.outfile


    basic_settings = {

        # You can set default values as key-value pairs.


        'NODE_FILL_COLOR': '#6AACB8',
        'NODE_SIZE': 55,
        'NODE_BORDER_WIDTH': 2,
        'NODE_LABEL_COLOR': '#555555',

        'EDGE_WIDTH': 2,
        #'EDGE_TRANSPARENCY': 100,
        #'EDGE_STROKE_UNSELECTED_PAINT': '#333333',

        'NETWORK_BACKGROUND_PAINT': 'white'

    }


    cytoscape=cyrest.cyclient()
    cytoscape.version()


    ## create network
    net = create_network(tmp_file)
    dz = net.get_node_table()

    if node_scores:
        dz = pd.merge(scores, dz, on = 'name')

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
    cy.style.apply(my_style, net)



    ### fix later
    if node_scores:

        j = 1 ## fix_later
        #scores = net.get_node_column('name')
        #scores = scores.astype(float)

        #if scores.min() == 0:
        #    minvalue = -2
        #else:
        #    minvalue = scores.min()
        #if scores.max() == 0:
        #    maxvalue = 2
        #else:
        #    maxvalue = scores.max()
        #color_gradient = StyleUtil.create_3_color_gradient(min=minvalue, max=maxvalue, colors=('blue','white', 'red'))



    degree_to_size = StyleUtil.create_slope(min=degrees.min(), max=degrees.max(),
                                            values=(30, 80))
    my_style.create_continuous_mapping(column=col, vp='NODE_FILL_COLOR',
                                        col_type='Double',
                                        points=color_gradient)
    my_style.create_continuous_mapping(column='degree.layout', vp='NODE_SIZE',
                                        col_type='Integer',
                                        points=degree_to_size)
    cy.style.apply(my_style, net)

    fig=cytoscape.networks.getFirstImageAsPng(networkId=cytoscape.network.get()["SUID"],
                                                h=None)
    #Image(fig.content)
    #with open("../output/networks_cmap/high_%s.png"%col, "wb") as png:
    #    png.write(fig.content)

    #cytoscape.network.deselect(nodeList='all', edgeList='all')
    cytoscape.view.export(options="PDF",outputFile=out_figure)
