# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 21:12:38 2020

@author: scabini
"""

import networkx as nx
# import model as md
import numpy as np
from comunities import *
from random import choices


def createGraph(parameters):
    np.random.seed(parameters['seed'])
    random.seed(parameters['seed'])
    G = nx.Graph()
    G.add_nodes_from(range(0, parameters['n_nodes']))
    
    #### statistical distributions to model the nodes (source: IBGE 2019)
   
    age_dist = parameters['age_dist']
    
    node_ages = choices([0, 1, 2, 3, 4, 5], age_dist[:], k=parameters['n_nodes'])
    
    node_state = []
    for i in range(0,parameters['n_nodes']):
        node_state.append(0)  # state id 0 = 'susceptible'
        
    states = dict(zip(list(G.nodes),node_state))
    ages = dict(zip(list(G.nodes),node_ages))
    days_infected = dict(zip(list(G.nodes),np.zeros((G.order()))))
    
    nx.set_node_attributes(G, ages, 'age')  # add the age attribute to nodes (age range index)
    nx.set_node_attributes(G, states, 'condition')  # initial state of all nodes is susceptible
    nx.set_node_attributes(G, days_infected, 'days infected')  # days since the node is infected
    
    for layer in parameters['layers_0']:
        G = include_layer(G, layer, parameters)
    
    return G


# remove layers. "layer" is an integer between 1 and qty of layers
def remove_layer(G, layer):
    # print("Removendo camada: ", layer)
    layers = {'casas': 1, 'escolas': 2, 'igrejas': 3, 'transporte': 4, 'aleatorio': 5, 'trabalho': 6}  # <ctm> aleatorio -> random
    layer = layers[layer]    
    edges = G.edges()    
    for edge in edges:
        a=edges[edge[0], edge[1]]
        if a['layer'] == layer:
            G.remove_edge(edge[0], edge[1])
    return G


def change_layer(G, layer, new_weight):
    # print("Alterando camada: ", layer)
    layers = {'casas': 1, 'escolas': 2, 'igrejas': 3, 'transporte': 4, 'aleatorio': 5, 'trabalho': 6}
    layer = layers[layer]    
    edges = G.edges()    
    for edge in edges:
        a=edges[edge[0], edge[1]]
        if a['layer'] == layer:
            G[edge[0]][edge[1]]['weight'] = new_weight
            
    return G


def scale_layer(G, layer, scale):
    # print("Alterando camada: ", layer)
    layers = {'casas': 1, 'escolas': 2, 'igrejas': 3, 'transporte': 4, 'aleatorio': 5, 'trabalho': 6}
    layer = layers[layer]    
    edges = G.edges()    
    for edge in edges:
        a=edges[edge[0], edge[1]]
        if a['layer'] == layer:
            G[edge[0]][edge[1]]['weight'] = G[edge[0]][edge[1]]['weight']*scale
            
    return G


def include_layer(G, layer, parameters):
    layers = {'casas': 1, 'escolas': 2, 'igrejas': 3, 'transporte': 4, 'aleatorio': 5, 'trabalho': 6}
    # print("Incluindo camada: ", layer)
    if layers[layer] == 1:
        G = generate_layer1(G, parameters)  # casas : houses
    elif layers[layer] == 2:
        G = generate_layer_school(G, parameters)  # escolas : schools
    elif layers[layer] == 3:
        G = generate_layer_church(G, parameters)  # templos : churches
    elif layers[layer] == 4:
        G = generate_layer_transport(G, parameters)  # transporte : transport
    elif layers[layer] == 5:
        G = generate_layer_random(G, parameters)  # aleatorio : random
    elif layers[layer]== 6:
        G = generate_layer_work(G, parameters)  # trabalho : work
            
    return G


# isolation for mild/moderate symptoms = stay home
def isolate_node(G, node, how='hospital - total'):
    # how can be 'no', in which case it does nothing and returns the original graph
    
    edges = G.edges(node)
    edges = [i for i in edges]
    alledges = G.edges()
    # optimistic hospital case = isolates totally = removes all edges
    # this function is also used to remove recovered and dead from the network
    if 'hospital - total' in how:
        for edge in edges:
            G.remove_edge(edge[0], edge[1])
    
    # optimistic home insulation, fully insulate home= remove all but 1 layers
    elif 'home - total' in how:
        for edge in edges:
            a=alledges[edge]            
            if a['layer'] != 1:
                G.remove_edge(edge[0], edge[1])
                
    # pessimistic home isolation, keep home and random connections
    elif 'home - partial' in how:
        for edge in edges:
            a=alledges[edge]            
            if a['layer'] != 1 and a['layer'] != 5:
                G.remove_edge(edge[0], edge[1])
    
    # pessimistic hospital case, keeps random edges
    elif 'hospital - partial' in how:
        for edge in edges:
            a=alledges[edge]  
            if a['layer'] != 5:
                G.remove_edge(edge[0], edge[1])
    return G
