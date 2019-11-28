#This module is for temporary methods

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

#Transition System plot
def TS_plot(log1,log2, merg=True):
    '''

    :param log1: a ist of lists : log1 = [['a', 'e','f','d']] * 10 + [['b', 'a', 'c']] * 4
    :param log2: a list of lists: [['c', 'b','a','f','a']] * 10 + [['b', 'a', 'c','b','a']] *1
    :param merg: Plotting together or separate
    :return:
    '''

    #Plotting two logs together
    if merg == True:
        log= log1+log2
        uniq_event =list(set([k for t in log for k in t]))

        #Create adjuncy matrix
        adj_matrix = np.zeros((len(uniq_event),len(uniq_event)))

        #filling the matrix
        for t in log:
            for i in range(len(t)-1):
                row_ind = uniq_event.index(t[i])
                col_ind = uniq_event.index(t[i+1])
                adj_matrix[row_ind][col_ind]=1

    #Creating a directed graph
    G = nx.DiGraph()

    #Adding node
    #G.add_nodes_from(uniq_event)

    # nodes=[]
    # for e in uniq_event:
    #     nodes.append(dot.node(e))


    #adding edges
    temp=[]
    for t in log:
        for i in range(len(t) - 1):
            if([t[i],t[i+1]] not in temp):
                #print "Temmp:",temp
                temp.append((t[i], t[i + 1]))
                #G.add_edge(t[i],t[i+1])

    G.add_edges_from(temp)

    nx.draw(G)
    plt.show()
    return G
    # # Graphviz must be installed in order to save output in pdf format
    # if platform.system() == "Windows":
    #     os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
    # elif platform.system() == "Linux":
    #     if ('graphviz' not in sys.modules):
    #         os.system('sudo apt-get install graphviz')
    #
    # dot.render('test-output/'+time.strftime("%m/%d/%Y, %H%M%S") +".pdf", view=True)

##########################################################################