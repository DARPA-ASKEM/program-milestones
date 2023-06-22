# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 10:05:34 2020

@author: scabini
"""
import networkx as nx
import numpy as np
from collections import Counter
import random
from random import choices
from random import uniform
# import model 


# work-around to adjust cases where the rng is bad
def check_comunities(comunities, n):    
    chksum=0
    for comunity in comunities:
        chksum += comunity+1
    if chksum != n:
        diff = chksum - n        
        if diff < 0:
            for i in range(0, np.absolute(diff)):
                comunities.append(0)
        else:
            return False   
              
    return comunities


# layer 1 - neighborhood - groups with at least 1 adult and sizes according to
# the distribution of family sizes in Brazil
def generate_layer1(G, parameters): 
    np.random.seed(parameters['seed'])
    random.seed(parameters['seed'])
    ############### 2010 G1 DATA  [ DADOS G1 2010 ]
    fam_structure = parameters['fam_structure']
   
    prob_infection = parameters['prob_home']
    n_fam = int(np.round(parameters['n_nodes'] * 0.30))  # 30%=ratio number of families/population (estimated from ibge/g1 2010)

    comunities = False
    while comunities == False:
        comunities = choices(population=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], weights=fam_structure[:], k=n_fam)
        comunities = check_comunities(comunities, parameters['n_nodes'])
        
    if parameters['verbose']:
        print("Number of families by size: ", (Counter([i+1 for i in comunities])))  # <ctm-code> print("Quantidade de familias por tamanho: ", (Counter([i+1 for i in comunities])))
   
    nodes = list(G.nodes)  # list of nodes

    ages = nx.get_node_attributes(G, 'age')  # list of nodes with age group
    nodes_adults = []
    
    # selecting nodes with age group id >= 2 (adults)
    i = 0
    while len(nodes_adults) != len(comunities):
        if ages[i] >= 2:
            nodes_adults.append(i)
            nodes.remove(i)
        i += 1
            
    # nodes_adults = nodes_adults[0:len(comunities)]
    avg_size = 0
    localss = [0] * G.number_of_nodes()
   
    local = 0
   
    qtde_layers = (len(parameters['layers_0']))
    variable_prob = ((6-qtde_layers)*0.2)
    prob = ((1*(3*7/(24*7)))*prob_infection) 
    prob = prob+(prob*variable_prob)
    for comunity in comunities:  # for every community...
        # calculation of the probability of infection, depends on the size, time of
        # interaction and average number of people nearby (details in the worksheet)

        nodes_to_conect = []
        father = nodes_adults.pop()
        nodes_to_conect.append(father)  # select 1 adult
        for size in range(1, comunity+1):  # creates the community according to its size
            if nodes:  # while we still have...
                random_node = nodes.pop()
                nodes_to_conect.append(random_node)  # select random node
                localss[random_node] = local
            else:
                break
                
        # connect everyone with everyone in the nodes_to_connect group, edge weight=infection prob
        # this is the community/family
        avg_size += len(nodes_to_conect)
        for i in range(0, len(nodes_to_conect)):
            for j in range(i+1, len(nodes_to_conect)):
                G.add_edge(nodes_to_conect[i], nodes_to_conect[j], weight=prob, layer=1)
                
        local += 1
   
    if parameters['verbose']:
        print('Average family size: ', avg_size/len(comunities))  # <ctm-code> print('Tamanho medio de familia: ', avg_size/len(comunities))
        
    localss = dict(zip(list(G.nodes), localss))
    nx.set_node_attributes(G, localss, 'local')
    return G


# layer 2 - schools - age groups (0-17) of school sizes
# according to the average size of classrooms in the country
def generate_layer_school(G, parameters):
    np.random.seed(parameters['seed'])
    random.seed(parameters['seed'])
    # ranges to choose evenly during creation
    # prob_infection = [0.13, 0.25]  # infection prob varies with size variation
    tamanhos = [16, 30]  # Room size varies according to statistics elementary > high school  # <ctm> Portuguese: tamanhos -> English: sizes
        
    prob_infection = parameters['prob_school']
    ages = nx.get_node_attributes(G, 'age')  # list of nodes with age group
    nodes = []
    # selecting the nodes with age group id <= 1 (0-17)
    for i in range(0, len(ages)):
        if ages[i] <= 1:
            nodes.append(i)
            
    random.shuffle(nodes)
    n_salas = 0
        
    avg_size = 0    
    while nodes:
        n_salas += 1
        nodes_to_conect = []
        
        size = random.randint(tamanhos[0], tamanhos[1])
        
        # calculation of the probability of infection, depends on the size,
        # interaction time and average number of people nearby (details in the worksheet)
        prob = (5/(size))*((4*5)/(24*7))*prob_infection         
        
        for i in range(0, size):  # creates the community according to its size
            if nodes:  # while we still have...
                # as the ages are determined randomly, taking the last one on the tbm is random
                random_node = nodes.pop() 
                nodes_to_conect.append(random_node)  # select random node from nodes
            else:
                break
                
        # connect everyone with everyone in the nodes_to_connect group, edge weight=infection prob
        # this is the community
        avg_size += len(nodes_to_conect)
        for i in range(0, len(nodes_to_conect)):
            for j in range(i+1, len(nodes_to_conect)):            
                G.add_edge(nodes_to_conect[i], nodes_to_conect[j], weight=prob, layer=2)
                
    if parameters['verbose']:
        print('How many rooms:', n_salas, ' - Average size: ', avg_size / (n_salas))  # print('Qtde salas:', n_salas, ' - Tamanho medio: ', avg_size/(n_salas))
    
    return G


# layer 3 - religious activities/churches
def generate_layer_church(G, parameters):
    np.random.seed(parameters['seed'])
    random.seed(parameters['seed'])
    # ranges to choose evenly during creation
    prob_religiosos = parameters['qtde_religiao']  # half of the population of Christians in 2019, source: Veja  # <ctm> Portuguese 'qtde_religiao' -> English: 'qty_religion'
    prob_religiosos = int(np.round(prob_religiosos*G.order()))
    
    tamanhos = [(10,50), (51,80), (81,100)]  # size range of temples/churches or groups of people who attend them
    prob_sizes = [0.552786405, 0.292858506, 0.15435509]
    
    # tamanhos= [(50,100), (101,200), (201,300), (301,400), (401,500), (501,1000), (1001,5000), (5001,10000), (10001,50000), (50001,100000)]
    # prob_sizes = [0.5527864045, 0.2472135955, 0.0750908150, 0.0354664659, 0.0204130244, 0.0381586767, 0.0261059192, 0.0026340818, 0.0018020836, 0.0001818299]
    
    prob_infection = parameters['prob_religion']

    nodes = list(G.nodes)
    random.shuffle(nodes)
    nodes = nodes[0:prob_religiosos]
    n_templos = 0
        
    avg_size = 0    
    while nodes:
        n_templos += 1
        nodes_to_conect = []
        
        size = random.choices(tamanhos, prob_sizes, k=1)  # Portuguese: tamanhos -> English: sizes
        size = random.randint(size[0][0], size[0][1])
        
        # calculation of the probability of infection, depends on the size, time of
        # interaction and average number of close people
        prob = (6/(size))*(2/(24*7))*prob_infection

        for i in range(0,size):  # creates the community according to its size
            if nodes:  # while we still have...
                random_node = nodes.pop() 
                nodes_to_conect.append(random_node)  # select random node from nodes
            else:
                break
                
        # connect everyone with everyone in the nodes_to_connect group, edge weight=infection prob
        # this is the community
        avg_size += len(nodes_to_conect)
        for i in range(0, len(nodes_to_conect)):
            for j in range(i+1, len(nodes_to_conect)):              
                G.add_edge(nodes_to_conect[i], nodes_to_conect[j], weight=prob, layer=3)
    
    if parameters['verbose']:
        print('Number of temples:', n_templos, ' - Average size: ', avg_size/(n_templos))  # <ctm-code> print('Qtde de templos:', n_templos, ' - Tamanho medio: ', avg_size/(n_templos))
    
    return G


# layer 4 - public transport
def generate_layer_transport(G, parameters):
    np.random.seed(parameters['seed'])
    random.seed(parameters['seed'])

    # ranges to choose evenly during creation
    prob_tpublico = parameters['qtde_transporte']  # fraction of population using public transport
    prob_tpublico = int(np.round(prob_tpublico*G.order()))
    tamanhos = [10, 40]  # size range/qty of people traveling by bus/metro
    prob_infection = parameters['prob_transport']

    nodes = list(G.nodes)
    random.shuffle(nodes)
    nodes = nodes[0:prob_tpublico]
    n_v = 0
        
    avg_size = 0    
    while nodes:
        n_v += 1
        nodes_to_conect = []
        
        size = random.randint(tamanhos[0], tamanhos[1])

        # calculation of the probability of infection, depends on the size, time of
        # interaction and average number of close people
        prob = (8/(size))*((parameters['tempo_transporte']*7)/(24*7))*prob_infection  # <ctm> Portuguese: 'tempo_transporte' -> English: 'transport_time'
        
        for i in range(0, size):  # creates the community according to its size
            if nodes:  # while we still have...
                random_node = nodes.pop() 
                nodes_to_conect.append(random_node)  # select random node from nodes
            else:
                break
                
        # connect everyone with everyone in the nodes_to_connect group, edge weight=infection prob
        # this is the community
        avg_size += len(nodes_to_conect)
        for i in range(0, len(nodes_to_conect)):
            for j in range(i+1, len(nodes_to_conect)):
                G.add_edge(nodes_to_conect[i], nodes_to_conect[j], weight=prob, layer=4)
                
    if parameters['verbose']:                
        print('number of vehicles:', n_v, ' - Average size: ', avg_size/(n_v))  # <ctm-code> print('Qtde de veiculos:', n_v, ' - Tamanho medio: ', avg_size/(n_v))
    
    return G


# layer 5 - random social encounters, surfaces, and all that stuff
def generate_layer_random(G, parameters):
    np.random.seed(parameters['seed'])
    random.seed(parameters['seed'])
    # ranges to choose evenly during creation
    rand_interaction = [1, 10]  # range of random social interactions people may have
    qtde_interactions = int(np.round(G.order()*(sum(rand_interaction)/2)))  # average

    prob_infection = parameters['prob_random']
    # calculation of the probability of infection, depends on size, interaction time
    # and average number of close people
    prob = (1*((1/(24*7))))*prob_infection
        
    nodes = list(G.nodes)
    random.shuffle(nodes)
     
    for i in range(0, qtde_interactions):
        source = random.choice(nodes)
        target = random.choice(nodes)
        G.add_edge(source, target, weight=prob, layer=5)
        
    if parameters['verbose']:
        print('Number of random links:', qtde_interactions)  # <ctm-code> print('Qtde de links aleatorios:', qtde_interactions)
    
    return G


# layer 6 - work - working age groups of people,
# group sizes is a random [5, 50]
def generate_layer_work(G, parameters):
    np.random.seed(parameters['seed'])
    random.seed(parameters['seed'])    
    # ranges to choose evenly during creation
    # prob_infection = [0.13, 0.25]  # Infection prob varies with size variation
    tamanhos = [5, 30]  # Room size varies according to statistics elementary > high school
        
    prob_infection = parameters['prob_work']
    ages = nx.get_node_attributes(G, 'age')  # list of nodes with age group
    nodes = []
    # selecting the nodes with age group id <= 1 (0-17)
    for i in range(0, len(ages)):
        if ages[i] >= 2 and ages[i] <= 4:
            nodes.append(i)
            
    random.shuffle(nodes)
    n_empresas = 0  # <ctm> n_companies
        
    avg_size = 0    
    while nodes:
        n_empresas += 1
        nodes_to_conect = []
        
        size = random.randint(tamanhos[0], tamanhos[1])
        
        # calculation of the probability of infection, depending on the size, interaction time
        # and average number of people nearby (details in the worksheet)
        prob = (3/(size))*((8*5)/(24*7))*prob_infection         
        
        for i in range(0, size):  # creates the community according to its size
            if nodes:  # while we still have...
                # as the ages are determined randomly, taking the last one on the tbm is random
                random_node = nodes.pop() 
                nodes_to_conect.append(random_node)  # select random node from nodes
            else:
                break
                
        # connect everyone with everyone in the nodes_to_connect group, edge weight=infection prob
        # this is the community
        avg_size += len(nodes_to_conect)
        for i in range(0, len(nodes_to_conect)):
            for j in range(i+1, len(nodes_to_conect)):            
                G.add_edge(nodes_to_conect[i], nodes_to_conect[j], weight=prob, layer=6)
                
    if parameters['verbose']:                
        print('number of companies:', n_empresas, ' - Average size: ', avg_size/(n_empresas))  # <ctm-code> print('Qtde empresas:', n_empresas, ' - Tamanho medio: ', avg_size/(n_empresas))
    
    return G


















