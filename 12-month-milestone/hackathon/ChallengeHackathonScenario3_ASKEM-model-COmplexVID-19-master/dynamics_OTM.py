# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 08:37:33 2020

Dynamics of the epidemic

@author: scabini
"""

# <ctm> UTI = unidade de Tratamento Intensivo --> English: Intensive Care Unit (ICU)

import os
import sys
from multiprocessing import Pool
import networkx as nx
import model
from model import isolate_node, remove_layer, include_layer, change_layer, scale_layer
import numpy as np
import random
from collections import Counter
from comunities import *
import matplotlib.pyplot as plt


def flip_coin(probability):
    return (random.choices([True, False], [probability, 1-probability], k=1))[0]


def simulate(G, parameters, infected_t0=1, days=1, hospital_isolation='hospital - total', home_isolation='home - total'):
    layer_names = ['casas', 'aleatorio', 'trabalho', 'transporte', 'escolas', 'igrejas']  # order used

    layers = {'casas': 1, 'escolas': 2, 'igrejas': 3, 'transporte': 4, 'aleatorio': 5, 'trabalho': 6}
    np.random.seed(parameters['seed'])
    ##### dynamics parameters
    # States:
    # 0 - susceptible
    # possible infection states (4 types) Source: Elnara Negri (USP-SP), ref. 3 papers
    all_states = {'susceptible': 0, 'infected - asymptomatic': 1, 'infected - mild': 2, 'infected - severe': 3,
                  'infected - critical': 4, 'recovered': 5, 'dead': 6}
    infected_states = ['infected - asymptomatic', 'infected - mild', 'infected - severe', 'infected - critical']  # normal list for random choices
    # probability of infection states
    prob_states = [0.3, 0.55, 0.1, 0.05]
    # 30% - asymptomatic
    # 55% - mild/moderate symptoms
    # 10% - strong symptoms   -> hospitalized
    # 5%  - critical symptoms -> UTI (ICU)

    # 5 - recovered
    # 6 - dead
    ###########################
    # ASYMPTOMATIC RECOVERING PROBABILITIES
    as_recovering = 0.09  # after 8 days, there is a chance to recover, 100% more or less on the 18th

    # Evolution probability distributions
    days_start_symptoms = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    sample = np.random.lognormal(mean=1.621, sigma=0.418, size=1000)  # incubation time format (paper reference) ; Portuguese: amostragem = English: sampling
    prob_start_symptoms = ((np.histogram(sample, bins=days_start_symptoms, density=True))[0])
    days_start_symptoms = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

    # death probability distributions
    severe_death_prob = 0.15/5   # 15% mortality split between days 20-24 (5 days)
    critical_death_prob = 0.5/4  # 50% mortality split between days 21-24 (4 days)

    ###########################

    states = nx.get_node_attributes(G, 'condition')
    degree = G.degree()  # <ctm> replaced variable 'grau' with 'degree': Portuguese: grau = Enlglish: degree

    k = [i for i in degree]

    avg_k = np.mean([i[1] for i in k])
    max_k = np.max([i[1] for i in k])

    if parameters['verbose']:
        print('Middle degree : ', avg_k)  # <ctm-code> print('Grau medio : ', avg_k)
        print('Maximum degree : ', max_k)  # <ctm-code> print('Grau maximo : ', max_k)

    i=1

    ######### creating a group of early infected near the middle grade (day 0)
    valid_initial_infected = False
    infected_list = []

    while not valid_initial_infected:
        nodes = list(G.nodes)
        random.seed(parameters['seed'] - i*10000)
        random.shuffle(nodes)
        infected_list = []
        i0=0
        try:
            while i0<infected_t0:
                node = nodes.pop()
                if degree[node] >= avg_k-i and degree[node] <= avg_k+i:
                    states[node] = all_states[(random.choices(infected_states, prob_states, k=1))[0]]
                    infected_list.append(node)
                    i0 += 1
            valid_initial_infected = True
        except:
            i += 1
    nx.set_node_attributes(G, states, 'condition')
    # states = nx.get_node_attributes(G, 'condition')
    ###################################################################

    count = np.zeros((11, days))
    days_infected = nx.get_node_attributes(G, 'days infected')
    days_symptoms = np.ones((G.order()))*-1
    leitos_ocupados = 0
    utis_ocupadas = 0
    for day in range(0, days):
        iii = -1
        for acao in parameters['acoes']:
            iii +=1
            if day == acao:
                for layerr in parameters['layers_tirar'][iii]:  # <ctm> in this context, Portuguese tirar appears closes to English: out
                    G = remove_layer(G, layerr)

                for layerr in parameters['layers_por'][iii]:  # <ctm> in this context, Portuguese por, English could be: per, for, in
                    G = include_layer(G, layerr, parameters)

                changed = len(parameters['layers_por'][iii]) - len(parameters['layers_tirar'][iii])
                if changed > 0:
                    # added more layers than removed, then scale weights on house layer by +20%
                    variable_prob = changed * 0.2
                    G = scale_layer(G, 'casas', variable_prob)  # <ctm> Portuguese: casas -> English: houses
                elif changed < 0: #removed more layers than added
                    # removed more layers than added, then scale weights on house layer 
                    variable_prob = np.abs(changed)*0.2
                    G = scale_layer(G, 'casas', 1/variable_prob)


        dayly_new_cases = {'dead':0, 'recovered':0, 'undiagnosed':0, 'diagnosed':0, 'infected - asymptomatic':0, 'infected - mild':0, 'infected - severe':0, 'infected - critical':0}
        # print(len(G.edges()))
        ############# dynamics in relation to the number of days and types of infection
        remove_infected = []
        # print("lista: ", len(infected_list))
        for i in infected_list: #iterate through all infected nodes
        # for i in range(0, G.order()):
            ##################################################################
            if states[i] == all_states['infected - asymptomatic']:
                days_infected[i] +=1
                if days_infected[i] >= 8:  # asymptomatic can recover after 8 days
                    if flip_coin(as_recovering):
                        remove_infected.append(i)
                        dayly_new_cases['recovered']+=1
                        states[i] = all_states['recovered']
                        G=isolate_node(G, i, how='hospital - total') #use this to remove recovered from the network, since they can't infect other people under this model

            elif states[i] == all_states['infected - mild']:
                days_infected[i] +=1
                if days_infected[i] == days_symptoms[i]:
                    if flip_coin(0.2):
                        dayly_new_cases['diagnosed'] += 1  # whereas only 20% of these are tested
                    else:
                        dayly_new_cases['undiagnosed'] += 1

                    G=isolate_node(G, i, how=home_isolation)  # isolate themselves at home when symptoms start

                elif days_infected[i] >= 10 and days_infected[i] > days_symptoms[i]:
                    # mild can recover evenly between day 10 and day 20
                    if flip_coin(1/(20-days_symptoms[i])):
                        remove_infected.append(i)
                        dayly_new_cases['recovered'] += 1
                        states[i] = all_states['recovered']
                        G = isolate_node(G, i, how='hospital - total') #recovered, so remove from network

                    elif days_infected[i] >= 20:
                        remove_infected.append(i)
                        dayly_new_cases['recovered']+=1
                        states[i] = all_states['recovered']
                        G = isolate_node(G, i, how='hospital - total') #recovered, so remove from network

            elif states[i] == all_states['infected - severe']:
                if days_symptoms[i] < 9:
                    days_symptoms[i] = 9

                days_infected[i] +=1
                if days_infected[i] == days_symptoms[i]:
                    leitos_ocupados = leitos_ocupados+1  # <ctm> Portuguese: "leitos ocupados" -> English:  "occupied beds"
                    dayly_new_cases['diagnosed'] += 1  # diagnosed cases = severe and critical when going to hospital
                    G = isolate_node(G, i, how=hospital_isolation)  # isolate themselves in the hospital

                elif days_infected[i] >= 20:
                    if flip_coin(severe_death_prob):
                        remove_infected.append(i)
                        leitos_ocupados = leitos_ocupados-1
                        dayly_new_cases['dead'] += 1
                        states[i] = all_states['dead']
                        G = isolate_node(G, i, how='hospital - total') #dead, remove from network
                    elif days_infected[i] >= 25:
                        remove_infected.append(i)
                        leitos_ocupados = leitos_ocupados-1
                        dayly_new_cases['recovered'] += 1
                        states[i] = all_states['recovered']
                        G = isolate_node(G, i, how='hospital - total') #recovered, remove from network

            elif states[i] == all_states['infected - critical']:
                if days_symptoms[i] < 8:
                    days_symptoms[i] = 8

                days_infected[i] +=1
                if days_infected[i] == days_symptoms[i]:
                    utis_ocupadas = utis_ocupadas+1  # <ctm> Portuguese "utis_ocupadas" -> English: "busy utils" ("utis" by itself == "useful")
                    dayly_new_cases['diagnosed'] += 1  # diagnosed cases = severe and critical when going to hospital
                    G = isolate_node(G, i, how=hospital_isolation)
                elif days_infected[i] >= 21:
                    if flip_coin(critical_death_prob):
                        remove_infected.append(i)
                        utis_ocupadas = utis_ocupadas-1
                        dayly_new_cases['dead'] += 1
                        states[i] = all_states['dead']
                        G=isolate_node(G, i, how='hospital - total')
                    elif days_infected[i] >= 25:
                        remove_infected.append(i)
                        utis_ocupadas = utis_ocupadas-1
                        dayly_new_cases['recovered'] += 1
                        states[i] = all_states['recovered']
                        G=isolate_node(G, i, how='hospital - total')

        nx.set_node_attributes(G, days_infected, 'days infected')  # initial state of all nodes is susceptible

        ############################################    
        new_states = dict(states)
        # print("remove: ", len(remove_infected))

        infected_list = [x for x in infected_list if x not in remove_infected]
        # print("lista: ", len(infected_list))

        edges = G.edges()
        include_infected = []
        for i in infected_list:
            edges_i = G.edges(i)
            for edge in edges_i:
                link = edges[edge[0], edge[1]]
                ###########standard infection
                if states[edge[0]] >= all_states['infected - asymptomatic'] and states[edge[0]] <= all_states['infected - critical']:  # if it's some kind of infected
                    if all_states['susceptible'] == states[edge[1]]:  # if adjacent is susceptible
                        if flip_coin(link['weight']):  # try to infect
                            estado = (random.choices(infected_states, prob_states, k=1))[0]  # <ctm> Portuguese: "estado" -> English: "state"
                            new_states[edge[1]] = all_states[estado]  # randomly determines a case type
                            days_symptoms[edge[1]] = (random.choices(days_start_symptoms, prob_start_symptoms, k=1))[0]  # determines when symptoms will start (ignores asymptomatic case)
                            dayly_new_cases[estado] += 1  # daily cases
                            include_infected.append(edge[1])

                elif states[edge[1]] >= all_states['infected - asymptomatic'] and states[edge[1]] <= all_states['infected - critical']:  # if it's some kind of infected
                    if all_states['susceptible'] == states[edge[0]]:  # if adjacent is susceptible
                        if flip_coin(link['weight']):  # try to infect
                            estado = (random.choices(infected_states, prob_states, k=1))[0]
                            new_states[edge[0]] = all_states[estado]  # randomly determines a case type
                            days_symptoms[edge[0]] = (random.choices(days_start_symptoms, prob_start_symptoms, k=1))[0]  # determines when symptoms will start (ignores asymptomatic case)
                            dayly_new_cases[estado] += 1  # daily cases
                            include_infected.append(edge[0])

        infected_list = infected_list + include_infected
        if parameters['verbose']:
            print('Day ', day+1, ': ', dayly_new_cases)

        nx.set_node_attributes(G, new_states , 'condition')
        states = nx.get_node_attributes(G, 'condition')

        st = list(states.values())
        count[0][day] = st.count(0)
        count[1][day] = dayly_new_cases['infected - asymptomatic']
        count[2][day] = dayly_new_cases['infected - mild']
        count[3][day] = dayly_new_cases['infected - severe']
        count[4][day] = dayly_new_cases['infected - critical']
        count[5][day] = dayly_new_cases['infected - asymptomatic'] + dayly_new_cases['undiagnosed']
        count[6][day] = dayly_new_cases['diagnosed']
        count[7][day] = dayly_new_cases['dead']
        count[8][day] = dayly_new_cases['recovered']
        count[9][day] = np.max((0,leitos_ocupados))
        count[10][day] = np.max((0,utis_ocupadas))

    return G, count





