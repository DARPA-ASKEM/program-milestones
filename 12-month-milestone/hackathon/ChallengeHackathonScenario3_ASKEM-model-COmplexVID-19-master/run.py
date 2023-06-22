# -*- coding: utf-8 -*-
"""
Created on Thu May  7 16:44:23 2020

<lang:english>
Main code to perform an experiment. Run in parallel using the
     number of colors - 2, each thread is an iteration, final results are the average
     Here it is fixed at 10 iterations, in the article we used 100 (ideal, but expensive)
     The script has several parameters that must be set manually, the
     comments should guide you. Any questions please contact:
         Leonardo Scabini, scabini@ifsc.usp.br
</lang:english>

<lang:portuguese>
Codigo principal para realizar um experimento. Roda em paralelo usando o
    numero de cores - 2, cada thread é uma iteração, resultados finais sao a media
    Aqui esta fixo em 10 iterações, no artigo usamos 100 (ideal, mas custoso) 
    O script tem varios parametros que devem ser setados manualmente, os 
    comentários devem guiá-lo. Qualquer dúvida entre em contato com:
        Leonardo Scabini, scabini@ifsc.usp.br
</lang:portuguese>

@author: scabini
"""

import os
import sys
from multiprocessing import Pool
import networkx as nx
import model 
from model import isolate_node
import numpy as np
import random
from collections import Counter
from comunities import *
# from dynamics import *
from dynamics_OTM import *
import random
import matplotlib.pyplot as plt
import pickle
from time import perf_counter 
import matplotlib.ticker as ticker

age_dist = np.zeros((6))
fam_structure = np.zeros((10))

####################### PARAMETERS #################################################################################
cidade = "Sao carlos"
n = 1000 #251983  # city population
# we test with n up to 250k, if the population is much larger than that, simulate
# with n=100k and change the scale factor
populacao_br = 57 #n  # scale factor to approximate the results according to
# the actual population of the study site *(see paper)

beta_global = 0.2  # article default = 0.3, simulates reduction of contagion rate in existing activities

# age distribution, sum must be 1
age_dist[0] = 0.18  # 0-13  - schools
age_dist[1] = 0.06  # 14-17 - schools
age_dist[2] = 0.11  # 18-24 - company
age_dist[3] = 0.23  # 25-39 - company
age_dist[4] = 0.26  # 40-59 - company
age_dist[5] = 0.16  # 60+

# family size distribution, sum must be 1
fam_structure[0] = 0.12  # 1 people -
fam_structure[1] = 0.22  # 2 people
fam_structure[2] = 0.25  # 3 people
fam_structure[3] = 0.21  # 4 people
fam_structure[4] = 0.11  # 5 people-
fam_structure[5] = 0.05  # 6 people-
fam_structure[6] = 0.02  # 7 people-
fam_structure[7] = 0.01  # 8 people-
fam_structure[8] = 0.006  # 9 people-
fam_structure[9] = 0.004  # 10 people-

qtde_religiao = 0.4  # fraction of the population that goes to church once a week
qtde_transporte = 0.36  # fraction of the population that uses public transport
tempo_transporte = 1.2   # time per day in hours, spent on transport
primeiro_infectado = "February 26" #"March 18"
layers_0 = ['casas', 'aleatorio', 'trabalho', 'transporte', 'escolas', 'igrejas']  # initial layers of the model, in the case of starting in the SP quarantine, 4
## home, random, work, transportation, school, church
acoes = ["March 24"] #, "May 26"]  # vector with the days of interventions, ex 24/3 of SP
layers_tirar = [["transporte", "escolas", "igrejas"], ["trabalho"]]  # list of lists, layers to remove in each action
layers_por = [[], []]  # lists of lists, layers to add on each action

repetitions = 5  # simulation repetitions, final results are the average. In the paper we did 100 repetitions (costly)

layer_names = ['casas', 'aleatorio', 'trabalho', 'transporte', 'escolas', 'igrejas']

# facilities, list with the days of the year (this format can be simulated until December of this year)
ano = {'January 1': 0, 'January 2': 1, 'January 3': 2, 'January 4': 3, 'January 5': 4, 'January 6': 5, 'January 7': 6, 'January 8': 7, 'January 9': 8, 'January 10': 9, 'January 11': 10, 'January 12': 11, 'January 13': 12, 'January 14': 13, 'January 15': 14, 'January 16': 15, 'January 17': 16, 'January 18': 17, 'January 19': 18, 'January 20': 19, 'January 21': 20, 'January 22': 21, 'January 23': 22, 'January 24': 23, 'January 25': 24, 'January 26': 25, 'January 27': 26, 'January 28': 27, 'January 29': 28, 'January 30': 29, 'January 31': 30, 'February 1': 31, 'February 2': 32, 'February 3': 33, 'February 4': 34, 'February 5': 35, 'February 6': 36, 'February 7': 37, 'February 8': 38, 'February 9': 39, 'February 10': 40, 'February 11': 41, 'February 12': 42, 'February 13': 43, 'February 14': 44, 'February 15': 45, 'February 16': 46, 'February 17': 47, 'February 18': 48, 'February 19': 49, 'February 20': 50, 'February 21': 51, 'February 22': 52, 'February 23': 53, 'February 24': 54, 'February 25': 55, 'February 26': 56, 'February 27': 57, 'February 28': 58, 'February 29': 59, 'March 1': 60, 'March 2': 61, 'March 3': 62, 'March 4': 63, 'March 5': 64, 'March 6': 65, 'March 7': 66, 'March 8': 67, 'March 9': 68, 'March 10': 69, 'March 11': 70, 'March 12': 71, 'March 13': 72, 'March 14': 73, 'March 15': 74, 'March 16': 75, 'March 17': 76, 'March 18': 77, 'March 19': 78, 'March 20': 79, 'March 21': 80, 'March 22': 81, 'March 23': 82, 'March 24': 83, 'March 25': 84, 'March 26': 85, 'March 27': 86, 'March 28': 87, 'March 29': 88, 'March 30': 89, 'March 31': 90, 'April 1': 91, 'April 2': 92, 'April 3': 93, 'April 4': 94, 'April 5': 95, 'April 6': 96, 'April 7': 97, 'April 8': 98, 'April 9': 99, 'April 10': 100, 'April 11': 101, 'April 12': 102, 'April 13': 103, 'April 14': 104, 'April 15': 105, 'April 16': 106, 'April 17': 107, 'April 18': 108, 'April 19': 109, 'April 20': 110, 'April 21': 111, 'April 22': 112, 'April 23': 113, 'April 24': 114, 'April 25': 115, 'April 26': 116, 'April 27': 117, 'April 28': 118, 'April 29': 119, 'April 30': 120, 'May 1': 121, 'May 2': 122, 'May 3': 123, 'May 4': 124, 'May 5': 125, 'May 6': 126, 'May 7': 127, 'May 8': 128, 'May 9': 129, 'May 10': 130, 'May 11': 131, 'May 12': 132, 'May 13': 133, 'May 14': 134, 'May 15': 135, 'May 16': 136, 'May 17': 137, 'May 18': 138, 'May 19': 139, 'May 20': 140, 'May 21': 141, 'May 22': 142, 'May 23': 143, 'May 24': 144, 'May 25': 145, 'May 26': 146, 'May 27': 147, 'May 28': 148, 'May 29': 149, 'May 30': 150, 'May 31': 151, 'June 1': 152, 'June 2': 153, 'June 3': 154, 'June 4': 155, 'June 5': 156, 'June 6': 157, 'June 7': 158, 'June 8': 159, 'June 9': 160, 'June 10': 161, 'June 11': 162, 'June 12': 163, 'June 13': 164, 'June 14': 165, 'June 15': 166, 'June 16': 167, 'June 17': 168, 'June 18': 169, 'June 19': 170, 'June 20': 171, 'June 21': 172, 'June 22': 173, 'June 23': 174, 'June 24': 175, 'June 25': 176, 'June 26': 177, 'June 27': 178, 'June 28': 179, 'June 29': 180, 'June 30': 181, 'July 1': 182, 'July 2': 183, 'July 3': 184, 'July 4': 185, 'July 5': 186, 'July 6': 187, 'July 7': 188, 'July 8': 189, 'July 9': 190, 'July 10': 191, 'July 11': 192, 'July 12': 193, 'July 13': 194, 'July 14': 195, 'July 15': 196, 'July 16': 197, 'July 17': 198, 'July 18': 199, 'July 19': 200, 'July 20': 201, 'July 21': 202, 'July 22': 203, 'July 23': 204, 'July 24': 205, 'July 25': 206, 'July 26': 207, 'July 27': 208, 'July 28': 209, 'July 29': 210, 'July 30': 211, 'July 31': 212, 'August 1': 213, 'August 2': 214, 'August 3': 215, 'August 4': 216, 'August 5': 217, 'August 6': 218, 'August 7': 219, 'August 8': 220, 'August 9': 221, 'August 10': 222, 'August 11': 223, 'August 12': 224, 'August 13': 225, 'August 14': 226, 'August 15': 227, 'August 16': 228, 'August 17': 229, 'August 18': 230, 'August 19': 231, 'August 20': 232, 'August 21': 233, 'August 22': 234, 'August 23': 235, 'August 24': 236, 'August 25': 237, 'August 26': 238, 'August 27': 239, 'August 28': 240, 'August 29': 241, 'August 30': 242, 'August 31': 243, 'September 1': 244, 'September 2': 245, 'September 3': 246, 'September 4': 247, 'September 5': 248, 'September 6': 249, 'September 7': 250, 'September 8': 251, 'September 9': 252, 'September 10': 253, 'September 11': 254, 'September 12': 255, 'September 13': 256, 'September 14': 257, 'September 15': 258, 'September 16': 259, 'September 17': 260, 'September 18': 261, 'September 19': 262, 'September 20': 263, 'September 21': 264, 'September 22': 265, 'September 23': 266, 'September 24': 267, 'September 25': 268, 'September 26': 269, 'September 27': 270, 'September 28': 271, 'September 29': 272, 'September 30': 273, 'October 1': 274, 'October 2': 275, 'October 3': 276, 'October 4': 277, 'October 5': 278, 'October 6': 279, 'October 7': 280, 'October 8': 281, 'October 9': 282, 'October 10': 283, 'October 11': 284, 'October 12': 285, 'October 13': 286, 'October 14': 287, 'October 15': 288, 'October 16': 289, 'October 17': 290, 'October 18': 291, 'October 19': 292, 'October 20': 293, 'October 21': 294, 'October 22': 295, 'October 23': 296, 'October 24': 297, 'October 25': 298, 'October 26': 299, 'October 27': 300, 'October 28': 301, 'October 29': 302, 'October 30': 303, 'October 31': 304, 'November 1': 305, 'November 2': 306, 'November 3': 307, 'November 4': 308, 'November 5': 309, 'November 6': 310, 'November 7': 311, 'November 8': 312, 'November 9': 313, 'November 10': 314, 'November 11': 315, 'November 12': 316, 'November 13': 317, 'November 14': 318, 'November 15': 319, 'November 16': 320, 'November 17': 321, 'November 18': 322, 'November 19': 323, 'November 20': 324, 'November 21': 325, 'November 22': 326, 'November 23': 327, 'November 24': 328, 'November 25': 329, 'November 26': 330, 'November 27': 331, 'November 28': 332, 'November 29': 333, 'November 30': 334, 'December 1': 335, 'December 2': 336, 'December 3': 337, 'December 4': 338, 'December 5': 339, 'December 6': 340, 'December 7': 341, 'December 8': 342, 'December 9': 343, 'December 10': 344, 'December 11': 345, 'December 12': 346, 'December 13': 347, 'December 14': 348, 'December 15': 349, 'December 16': 350, 'December 17': 351, 'December 18': 352, 'December 19': 353, 'December 20': 354, 'December 21': 355, 'December 22': 356, 'December 23': 357, 'December 24': 358, 'December 25': 359, 'December 26': 360, 'December 27': 361, 'December 28': 362, 'December 29': 363, 'December 30': 364, 'December 31': 365}
ano_rev = inv_map = {v: k for k, v in ano.items()}
days = 300 #- ano[primeiro_infectado]  # days to simulate, starts counting from the first infected
begin = ano[primeiro_infectado]
acoes1 = [ano[i]-begin for i in acoes]  # <ctm> Portuguese: acoes = English: actions <ctm>



# model parameters
parameters = dict(
    seed=999666,  # fixed seed for random operations
    age_dist=age_dist,  # age distribution
    fam_structure=fam_structure,  # family size distribution
    tempo_transporte=tempo_transporte,  # <ctm> time of transport
    qtde_transporte=qtde_transporte,  # <ctm> transport
    qtde_religiao=qtde_religiao,  # <ctm> religion
    layers_0=layers_0,  # model initial layers
    acoes=acoes1,  # <ctm> intervention dates
    layers_tirar=layers_tirar,  # <ctm> layers remove
    layers_por=layers_por,  # <ctm> layers add
    n_nodes=n,  # number of people used to estimate the % of the epidemic
    # the real probabilities are dynamic and depend on several things, which are
    # defined there within the creation of communities. This value here is "how much
    # to consider of these values: 0->nullifies, 1->original, 0.5->half, 2->double
    ## below betas are multiplied by layer-relevant modifier, within the comunities.py file, when each layer is instantiated/created in the graph
    prob_home=beta_global,  # original> whole family, 3hrs/day. Each layer removed increases home interaction by 25%
    prob_random=beta_global,  # original > 1 next, 1hrs/week, everyone has a chance of having 1 to 10 random connections
    prob_work=beta_global,  # original> 4 next/company size, 6hrs/day, 5 days/week, all 18 to 59 years old.
    prob_transport=beta_global,  # original > 5 next/vehicle size, 1.2hrs/day, 50% population (random)
    prob_school=beta_global,
    # original> 4 next/room size, 5hrs/day, 5 days/week, all population from 0 to 17 years old
    prob_religion=beta_global,  # original> 6 next/church size, 2hrs/week, 40% of the population (random)
    verbose=False  # print or not the information during construction and simulation
)


##################### INFECTION PARAMETERS
infected_t0 = 1  # initial infected
home_isolation = 'home - total'  # 'total' or 'partial', isolate the vertex completely at home (house edges stay), or keep the random ones as well
hospital_isolation = 'hospital - total'  # 'total' or 'partial', isolate the vertex in the hospital entirely or keep the random links
                             #   representing connections with people in the hospital


def fun(data):
    return " ".join([item for var in data for item in var])


def analyze(i):
    parameters['seed'] = i
    np.random.seed(parameters['seed'])
    random.seed(parameters['seed'])
    count=[]
    G = model.createGraph(parameters)
    G, count = simulate(G, parameters, infected_t0=infected_t0, days=days, hospital_isolation=hospital_isolation, home_isolation=home_isolation)
    G = []
    count = np.true_divide(count, n)
    return count


if __name__ == '__main__':
    print("Home network: ", parameters['layers_0'])  # <ctm-code> print("Rede inicial: ", parameters['layers_0'])
    print("Actions days: ", acoes)  # <ctm-code> print("Ações dias: ", acoes)
    print("Layers to insert: ", layers_por)  # <ctm-code> print("Camadas a inserir: ", layers_por)
    print("Layers to be removed: ", layers_tirar)  # <ctm-code> print("Camadas a tirar: ", layers_tirar)

    file = 'experimentos/' + cidade + '_REALISTIC_-' + str(parameters['n_nodes']) + '_reps-'
    file = file + str(repetitions) + '_beta-' + str(beta_global) + '_actions-' + fun(acoes)  # <ctm-code> file = file + str(repetitions) + '_beta-' + str(beta_global) + '_acoes-' + fun(acoes)
    file = file + '_add-(' +  fun(layers_por)  # <ctm-code> file = file + '_por-(' +  fun(layers_por)
    file = file + ')_remove-(' + fun(layers_tirar) +  ').pickle'  # <ctm-code> file = file + ')_tirar-(' + fun(layers_tirar) +  ').pickle'

    print(file)

    exists = os.path.isfile(file)
    if exists:
        with open(file, 'rb') as f:
            count = pickle.load(f)
            f.close()
    else:
        
        index_list = [i for i in range(1, repetitions+1)]
        processes = os.cpu_count()-2
        pool = Pool(processes)
        print('Running with', processes, 'threads')
        
        count=[]
        t1_start = perf_counter() 
        result = pool.map(analyze, iterable=index_list, chunksize=None)
        # result =analyze(1)

        t1_stop = perf_counter() 
        print("Spent ", t1_stop-t1_start, " seconds")
        
        count = np.zeros((11, days, repetitions))
        # count[:,:,0] = result
        i=0 
        for it in result:    
            count[:,:,i] = it  
            i+=1

        with open(file, 'wb') as f:
            pickle.dump(count, f)
    
    final = np.copy(count) 
    final3 = np.copy(count)*populacao_br
    
    final = final[:,:, final3[0,days-1,:]!=n-1]
    final3 = np.copy(final)*populacao_br
    
    count = np.mean(final, axis=2)  # iteration average
    std_mat = np.std(final, axis=2)  # standard deviation
    
    final2 = np.copy(count)  
    count =count*populacao_br
    std_mat = std_mat*populacao_br

    total_casos = np.cumsum(final[5], axis=0)*populacao_br
    std_und = np.std(total_casos, axis=1)
    und_casos = np.mean(total_casos, axis=1)
    
    total_casos = np.cumsum(final[6], axis=0)*populacao_br
    std_diag = np.std(total_casos, axis=1)
    diag_casos = np.mean(total_casos, axis=1)
    
    total_casos = np.cumsum(final[8], axis=0)*populacao_br
    recovered = np.mean(total_casos, axis=1)
    std_recv =  np.std(total_casos, axis=1)
    
    errorspace = 2
    plt.figure(0)
    plt.rcParams.update({'font.size': 15})
    
    eb=plt.errorbar(range(0, days), und_casos, yerr=std_und, lw=2, color='red', label='não diagnosticado', errorevery=errorspace)
    eb[-1][0].set_linewidth(1)
    
    eb=plt.errorbar(range(0, days), diag_casos, yerr=std_diag, lw=2, color='orange', label='diagnosticado', errorevery=errorspace)
    eb[-1][0].set_linewidth(1)
    
    eb=plt.errorbar(range(0, days), recovered, yerr=std_recv, lw=2, alpha=0.65, color='green', label='recuperado', errorevery=errorspace)
    eb[-1][0].set_linewidth(1)
    
    # plt.ylim(0, 5100000)
    # plt.xlim([0,210])
    plt.xlabel('Days since first case')  # <ctm-code> plt.xlabel('Dias desde o primeiro caso')
    plt.ylabel('Total cases')  # <ctm-code> plt.ylabel('Total de casos')
    ax = plt.gca()
    # ax.yaxis.set_major_formatter(ticker.EngFormatter())

    plt.xticks([0, acoes1[0], days],  ("0\n"+primeiro_infectado, str(acoes1[0]) + "\n" + acoes[0], str(days) +"\nDezembro 31"))

    # plt.yticks([])
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()
    # plt.savefig(path2 + "remove2_casos_total.pdf", dpi=500)
    
    
    
    plt.figure(1)
    plt.rcParams.update({'font.size': 15})
    
    # plt.plot((count[9]), label='leito', color='orange')
    # plt.plot((count[10]), label='UTI', color='red')    
    # plt.xticks([0, 34, 64, 95, 125, 156, 186, 217, 247, 278],  ('1º', 'Abr.', 'Maio', 'Jun.', "Jul.", "Ag.", "Set.", "Out.", "Nov.", "Dez."))

    peak = ((np.where(count[9] == np.amax(count[9])))[0] + (np.where(count[10] == np.amax(count[10])))[0])/2
    peak = int(peak[0])
    ICUbeds = np.round(count[10, peak])
    bedsstd = np.round(count[9, peak])
    print("Peak of occupied beds: ", ano_rev[ano[primeiro_infectado] + peak]  ," - UTIs ", ICUbeds, "(+-", np.round(std_mat[10,peak]) ,") - normal ", bedsstd, "(+-", np.round(std_mat[9,peak]) ,")")  # <ctm-code> print("Pico de leitos ocupados: ", ano_rev[ano[primeiro_infectado] + peak]  ," - UTIs ", ICUbeds, "(+-", np.round(std_mat[10,peak]) ,") - normal ", bedsstd, "(+-", np.round(std_mat[9,peak]) ,")")
    # plt.ylim(0, 280000)
    # plt.yticks([0, 43000, 86000, 135000, 215000])
    
    plt.errorbar(range(0, days), (count[9]), yerr=std_mat[9], label='normal', color='orange',zorder=0,errorevery=errorspace)
    plt.errorbar(range(0, days), (count[10]), yerr=std_mat[10], label='UTI/Respirator', color='red',zorder=0,errorevery=errorspace)  # <ctm-code> plt.errorbar(range(0, days), (count[10]), yerr=std_mat[10], label='UTI/Respirador', color='red',zorder=0,errorevery=errorspace)
    plt.grid(zorder=10)
    plt.xlabel('Days since first case')  # <ctm-code> plt.xlabel('Dias desde o primeiro caso')
    plt.ylabel('Occupied beds')  # <ctm-code> plt.ylabel('Leitos ocupados')
    
    plt.xticks([0, peak, days],  ("0\n"+primeiro_infectado, str(peak) + "\n" + ano_rev[begin+peak], str(days) +"\nDezembro 31"))

    # plt.xlabel('Vagas ocupadas')
    # plt.ylabel('Total')
    # plt.xlim([0,210])
    plt.legend(loc='upper left')
    ax = plt.gca()
    # ax.yaxis.set_major_formatter(ticker.EngFormatter())
    plt.tight_layout()
    plt.show()
    # plt.savefig(path2 + "donothing_leitos.pdf", dpi=500)
    
    
    data = count[6]
    plt.figure(2)
    plt.rcParams.update({'font.size': 16})
    # fig, (ax1) = plt.subplots(2)
    
    # ax1.scatter(250,1, label="no new cases")
    plt.bar(range(1, len(data)+1), data, color='blue',width=1.0)

    # plt.ylim(0, 41000)
    # plt.yticks([0,8700, 20000, 30000, 40000])

    peak = (np.where(data == np.amax(data)))[0]
    plt.xticks([0, peak[0], days],  ("0\n"+primeiro_infectado, str(peak[0]) + "\n" + ano_rev[begin+peak[0]], str(days) +"\nDezembro 31"))

    plt.xlabel('Days since first case')  # <ctm-code> plt.xlabel('Dias desde o primeiro caso')
    plt.ylabel('New cases')  # <ctm-code> plt.ylabel('Novos casos')
    plt.grid()
    # plt.legend()
    plt.tight_layout()
    plt.show()
    # plt.savefig(path2 + "casos_daily.pdf", dpi=500)

    
    data = count[7]
    plt.figure(3)
    plt.rcParams.update({'font.size': 16})
    # fig, (ax1) = plt.subplots(3)
    
    # ax1.scatter(281,0, label="no new deaths")
    plt.bar(range(1, len(data)+1), data, color='red', width=1.0)
    
    # plt.ylim(0, 8000)
    peak = (np.where(data == np.amax(data)))[0]
    
    plt.xticks([0, peak[0], days],  ("0\n"+primeiro_infectado, str(peak[0]) + "\n" + ano_rev[begin+peak[0]], str(days) +"\nDezembro 31"))

    plt.xlabel('Days since first case')  # <ctm-code> plt.xlabel('Dias desde o primeiro caso')
    plt.ylabel('New cases')  # <ctm-code> plt.ylabel('Novos casos')
    plt.grid()
    plt.tight_layout()
    plt.show()
    # plt.savefig(path2 + "mortes_daily.pdf", dpi=500)


   
    casos = np.sum(final[6], axis=0)*populacao_br
    casos = casos[casos != 0]
    casos2 = np.sum(final[5], axis=0)*populacao_br
    casos2 = casos2[casos2 != 0]
    mortes = np.sum(final[7], axis=0)*populacao_br
    mortes = mortes[mortes != 0]
    
    print('Statistics according to the total population:')  # <ctm-code> print('Estatistics de acordo com a população total:')
    print('Undiagnosed cases:        ',  np.round(np.mean(casos2)), '(+-', np.round(np.std(casos2)),')' ,' -  (', np.round(max(np.cumsum(final2[5]))*100, decimals=4), '%)')  # <ctm-code> print('Casos não diagnosticados: ',  np.round(np.mean(casos2)), '(+-', np.round(np.std(casos2)),')' ,' -  (', np.round(max(np.cumsum(final2[5]))*100, decimals=4), '%)')
    print('Diagnosed cases:          ',  np.round(np.mean(casos)), '(+-', np.round(np.std(casos)),')' ,' -  (', np.round(max(np.cumsum(final2[6]))*100, decimals=4), '%)')  # <ctm-code> print('Casos diagnosticados:     ',  np.round(np.mean(casos)), '(+-', np.round(np.std(casos)),')' ,' -  (', np.round(max(np.cumsum(final2[6]))*100, decimals=4), '%)')
    print('Deaths:                   ',  np.round(np.mean(mortes)), '(+-', np.round(np.std(mortes)),')' ,' -  (', np.round(max(np.cumsum(final2[7]))*100, decimals=4), '%)')  # <ctm-code> print('Mortes:                   ',  np.round(np.mean(mortes)), '(+-', np.round(np.std(mortes)),')' ,' -  (', np.round(max(np.cumsum(final2[7]))*100, decimals=4), '%)')
    print('Recovered:                ',  (max(np.cumsum(count[8]))), ' (', np.round(max(np.cumsum(final2[8]))*100, decimals=4), '%)')  # <ctm-code> print('Recuperados:              ',  (max(np.cumsum(count[8]))), ' (', np.round(max(np.cumsum(final2[8]))*100, decimals=4), '%)')
    
    
    
    
    
    
    
    
    
    
    
   