
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import numpy as np
import seaborn as sns
import Variables as v
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import os
import platform
import sys
from graphviz import Digraph
import time
from scipy import stats
from collections import Counter
import pprint
import Deviation as dv

###############################################################################
#Plotting three D of the three first principle components
def ThreeDPlot(event):
    fig = plt.figure(figsize=(9, 6))
    ax = Axes3D(fig)

    col = ['red'] * len(v.m2_reduced[:, 0]) + ['blue'] * len(v.m_reduced[:, 0])
    zs = np.concatenate((v.m2_reduced[:, 2], v.m_reduced[:, 2]))
    ys = np.concatenate((v.m2_reduced[:, 1], v.m_reduced[:, 1]))
    xs = np.concatenate((v.m2_reduced[:, 0], v.m_reduced[:, 0]))

    df = {'x': list(xs), 'y': list(ys), 'z': list(zs), 'col': col}
    df = pd.DataFrame(df)

    ax.scatter(xs, ys, zs, c=col, marker='o', data=df)
    ax.set_xlabel('First Principle Component')
    ax.set_ylabel('Second Principle Component')
    ax.set_zlabel('Third Principle Component')
    ax.set_title("Distribution for:"+str(event))
    plt.show()


###############################################################
##Temporary
def LDA_Plot():
    obs=100
    mean1= [1,1]
    cov1 = [[1,0],[0,1]]

    mean2 = [-5,-9]
    cov2 = [[9,7],[7,11]]

    x1= np.random.multivariate_normal(mean1,cov1, size=obs)
    x2= np.random.multivariate_normal(mean2,cov2,size=obs)

    X = np.concatenate((np.random.multivariate_normal(mean1,cov1, size=obs), np.random.multivariate_normal(mean2,cov2,size=obs) ))


    clf = LinearDiscriminantAnalysis()

    Y=['red']*obs + ['blue']*obs

    df = {'X': X[0], 'Y': Y}
    print df
    df = pd.DataFrame(df)
    print df

    X_lda= clf.fit_transform(X, Y)
    print "Coeff:", clf.coef_


    plt.plot(X_lda[0:obs, 0], 'gx')
    plt.plot(X_lda[obs:, 0], 'go')
    # plt.plot(X[0:obs,0], X[0:obs,1], 'rx')
    # plt.plot(X[obs:, 0], X[obs:, 1], 'bo')
    # plt.quiver([0, 0], clf.coef_[0])
    # plt.plot(x1, y1,'x')
    # plt.plot(x2,y2,'o')
    plt.axis('equal')
    plt.show()
################################################################
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
    dot = Digraph(comment='TS system')

    nodes=[]
    for e in uniq_event:
        nodes.append(dot.node(e,shape="rect"))


    #adding edges
    temp=[]
    for t in log:
        for i in range(len(t) - 1):
            if([t[i],t[i+1]] not in temp):
                #print "Temmp:",temp
                temp.append([t[i], t[i + 1]])
                dot.edge(t[i],t[i+1])

    # Graphviz must be installed in order to save output in pdf format
    if platform.system() == "Windows":
        os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
    elif platform.system() == "Linux":
        if ('graphviz' not in sys.modules):
            os.system('sudo apt-get install graphviz')

    dot.render('test-output/'+time.strftime("%m/%d/%Y, %H%M%S") +".pdf", view=True)

##########################################################################
def Log_Filter(log, id, deviant_elements=None):
    '''
    This method filters traces that do not contain any significant element
    :param log: a ist of lists : log1 = [['a', 'e','f','d']] * 10 + [['b', 'a', 'c']] * 4
    :param deviant_elements: a list of statistically significant elements. [('a','b'),('a','c'),...]
    :param id: #like {0: '0,1,2,3,4,5,6,7,8,9', 1:'12,13,....}
    :return:
    '''

    deviant_elements = [u.event[0] for u in v.deviant_obj if ((u.pValue < 0.05) or(u.class_proportion['proportion'] == None))]

    deviant_elements_freq = {}
    deviant_elements_f_score={}
    if (log == v.log1):
        for u in v.deviant_obj:
            if (u.pValue < 0.05):
                deviant_elements_freq[u.event[0]] = u.class_proportion['no_matrix1']
                deviant_elements_f_score[u.event[0]] = np.min((u.f1_avg['f1C1'], u.f1_avg['f1C2']))

    else:
        for u in v.deviant_obj:
            if(u.pValue<0.05):
                deviant_elements_freq[u.event[0]] = u.class_proportion['no_matrix2']
                deviant_elements_f_score[u.event[0]] = np.min((u.f1_avg['f1C1'], u.f1_avg['f1C2']))


    total_trace = np.sum([len(id[i].split(',')) for i in id.keys()])

    log_filtered=[]
    kgram = len(deviant_elements[0])
    for j in range(len(log)):
        trace=log[j]
    #for trace in log:
        flag=0
        count=0
        for i in range(len(trace)-kgram+1):
            if ( ( tuple(trace[i:i+kgram]) in deviant_elements) and flag==0):
                print float(len(id[j].split(",")))/ total_trace
                #if ( float(len(id[j].split(",")))/ total_trace > 0.01 ):
                condition = ( (float(len(id[j].split(","))) / total_trace > 0.01) and (deviant_elements_freq[tuple(trace[i:i+kgram])] >20) and(deviant_elements_f_score[tuple(trace[i:i+kgram])] !=0 ) )
                if ( condition ):
                    print "hi",j
                    print "Freq:",deviant_elements_freq[tuple(trace[i:i+kgram])] , tuple(trace[i:i+kgram])
                    print "trace:", trace
                    print "Fscore:", deviant_elements_f_score[tuple(trace[i:i+kgram])]
                    count+=1
                    print "The ratio :",float(count)/len(trace)
                    if(float(count)/len(trace)>0.19):
                        log_filtered.append(trace)
                    #flag=1

    #Computing the Entropy of the obtained graph
    log= log_filtered
    uniq_event =list(set([k for t in log for k in t]))

    #Create adjuncy matrix
    adj_matrix = np.zeros((len(uniq_event),len(uniq_event)))

    #filling the matrix
    for t in log:
        for i in range(len(t)-1):
            row_ind = uniq_event.index(t[i])
            col_ind = uniq_event.index(t[i+1])
            adj_matrix[row_ind][col_ind]=1

    print "The adjucency matrix is:", adj_matrix
    print "The entropy is:", stats.entropy(list(Counter(adj_matrix.flatten()).values()), base=2)
    #Calling TS_plot()
    TS_plot(log_filtered,log_filtered)


    # #----------------------------------------------------------------------
##########################################################################
#This method plot the Directly flow graph of traces containing deviant elements (it plots the the process model at Entorpy <0.7)
def Log_Filter_ControlFlow(log, id, deviant_elements=None):
    '''
    This method filters traces that do not contain any significant element in terms of controlflow
    :param log: a ist of lists : log1 = [['a', 'e','f','d']] * 10 + [['b', 'a', 'c']] * 4
    :param deviant_elements: a list of statistically significant elements. [('a','b'),('a','c'),...]
    :param id: #like {0: '0,1,2,3,4,5,6,7,8,9', 1:'12,13,....}
    :return:
    '''


    # Storing a dictionary for the elements and their information that will be showed on the process model
    v.deviant_element_plot = {}

    #Storing deviant elements names its like  ('Add#penalty', 'Payment')
    #deviant_elements = [u.event[0] for u in v.deviant_obj if ((u.pValue < 0.05) or(u.class_proportion['proportion'] == None))]
    deviant_elements = [u.event[0] for u in v.deviant_obj if (u.pValue < 0.05)]

    #Storing elements that happened in one class only
    one_class_element=[]

    #Storing the frequency and Minimum F-score for each element
    deviant_elements_freq = {}
    deviant_elements_f_score={}
    if (log == v.log1):
        for u in v.deviant_obj:
            if (u.pValue < 0.05):
                deviant_elements_freq[u.event[0]] = u.class_proportion['no_matrix1']
                deviant_elements_f_score[u.event[0]] = np.min((u.f1_avg['f1C1'], u.f1_avg['f1C2']))

            #Storing elements only happened in log1
            elif(u.class_proportion['proportion'] == None):
                one_class_element.append(u.event[0])
                deviant_elements_freq[u.event[0]] = u.class_proportion['no_matrix1']

    else:
        for u in v.deviant_obj:
            if(u.pValue<0.05):
                deviant_elements_freq[u.event[0]] = u.class_proportion['no_matrix2']
                deviant_elements_f_score[u.event[0]] = np.min((u.f1_avg['f1C1'], u.f1_avg['f1C2']))

                # Storing elements only happened in log2
            elif (u.class_proportion['proportion'] == None):
                one_class_element.append(u.event[0])
                deviant_elements_freq[u.event[0]] = u.class_proportion['no_matrix2']


    #The total number of traces.
    total_trace = np.sum([len(id[i].split(',')) for i in id.keys()])





    #Iterating over traces and filter those that do not have any deviant element
    iter=1
    deviation_ratio_per_trace=0.1
    trace_frequecny = 0.001
    entropy_record=[]
    while(iter==1):
        log_filtered=[]
        kgram = len(deviant_elements[0])
        for j in range(len(log)):
            trace=log[j]
        #for trace in log:
            count=0
            for i in range(len(trace)-kgram+1):
                if ( ( tuple(trace[i:i+kgram]) in deviant_elements)):
                    print float(len(id[j].split(",")))/ total_trace
                    condition = ( (float(len(id[j].split(","))) / total_trace > trace_frequecny) and (deviant_elements_freq[tuple(trace[i:i+kgram])] >20) and(deviant_elements_f_score[tuple(trace[i:i+kgram])] >0.05 ) )
                    if ( condition ):
                        print "Iter",j
                        print "Freq:",deviant_elements_freq[tuple(trace[i:i+kgram])] , tuple(trace[i:i+kgram])
                        print "trace:", trace
                        print "F-score:", deviant_elements_f_score[tuple(trace[i:i+kgram])]
                        count+=1
                        print "The ratio :",float(count)/len(trace)

                        # Storing the information of deviant elements for plotting
                        v.deviant_element_plot[tuple(trace[i:i + kgram])] = deviant_elements_freq[tuple(trace[i:i + kgram])]

                        if(float(count)/len(trace)>deviation_ratio_per_trace):
                            if (trace not in log_filtered):
                                print "Trace is added--------------------------------------------"
                                log_filtered.append(trace)



                #Finding traces that contain elements that happened in exclusively one log
                elif ((tuple(trace[i:i + kgram]) in one_class_element)):
                    condition = ((float(len(id[j].split(","))) / total_trace > 0.009) and (deviant_elements_freq[tuple(trace[i:i + kgram])] > 20))
                    if(condition):
                        print "One class Freq:", deviant_elements_freq[tuple(trace[i:i + kgram])], tuple(trace[i:i + kgram])
                        if (trace not in log_filtered):
                            log_filtered.append(trace)


                        #flag=1

        #Computing the Entropy of the obtained graph
        log= log_filtered
        uniq_event =list(set([k for t in log_filtered for k in t]))

        #Create adjacency matrix
        adj_matrix = np.zeros((len(uniq_event),len(uniq_event)))

        #filling the matrix
        for t in log_filtered:
            for i in range(len(t)-1):
                row_ind = uniq_event.index(t[i])
                col_ind = uniq_event.index(t[i+1])
                adj_matrix[row_ind][col_ind]=1
        #From adj matrix the entropy of the graph is computed
        model_entropy=stats.entropy(list(Counter(adj_matrix.flatten()).values()), base=2)
        print "The adj matrix is:", adj_matrix
        print "The entropy is:", model_entropy
        print "The normalized entropy is:", model_entropy/np.sum(adj_matrix)
        print "The Deviation ratio:",deviation_ratio_per_trace

        if(model_entropy <0.8):
            iter =0
        else:
            #Increasing the entropy slowly slowly by looking at its standard deviation (Adaptive increasing)
            entropy_record.append(model_entropy)
            deviation_ratio_per_trace+=0.01 + 0.001*1/(1+np.std(entropy_record))

            #Increasing the frequency threshold to include traces
            trace_frequecny +=0.001



    #Calling TS_plot()
    print "entropy record:",entropy_record
    print "the log filtered:", log_filtered
    #-----------------------------------------------------------------

    TS_plot(log_filtered,log_filtered)

    return log_filtered, v.deviant_element_plot, deviant_elements_freq

###########################################################################################
# #This method plot the Directly flow graph of traces containing deviant elements (it plots the the process model at Entorpy <0.7)
# def Log_Filter_ControlFlow_TimeStamp(log, id, deviant_elements=None):
#     '''
#     This method filters traces that do not contain any significant element in terms of controlflow
#     :param log: a ist of lists : log1 = [['a', 'e','f','d']] * 10 + [['b', 'a', 'c']] * 4
#     :param deviant_elements: a list of statistically significant elements. [('a','b'),('a','c'),...]
#     :param id: #like {0: '0,1,2,3,4,5,6,7,8,9', 1:'12,13,....}
#     :return:
#     '''
#
#
#     # Storing a dictionary for the elements and their information that will be showed on the process model
#     v.deviant_element_plot = {}
#     v.deviant_element_timeStamp_plot={}
#
#     #Storing deviant elements names its like  ('Add#penalty', 'Payment')
#     #deviant_elements = [u.event[0] for u in v.deviant_obj if ((u.pValue < 0.05) or(u.class_proportion['proportion'] == None))]
#     deviant_elements = [u.event[0] for u in v.deviant_obj if (u.pValue < 0.05)]
#     deviant_elements_timestamp = [u.event[0] for u in v.deviant_obj_timestamp if ( (u.class_proportion['no_matrix1']>20) and (u.class_proportion['no_matrix2']>20))]
#
#     #Storing elements that happened in one class only
#     one_class_element=[]
#
#     #Storing the frequency and Minimum F-score for each element
#     deviant_elements_freq = {}
#     deviant_elements_f_score={}
#     if (log == v.log1):
#         for u in v.deviant_obj:
#             if (u.pValue < 0.05):
#                 deviant_elements_freq[u.event[0]] = u.class_proportion['no_matrix1']
#                 deviant_elements_f_score[u.event[0]] = np.min((u.f1_avg['f1C1'], u.f1_avg['f1C2']))
#
#             #Storing elements only happened in log1
#             elif(u.class_proportion['proportion'] == None):
#                 one_class_element.append(u.event[0])
#                 deviant_elements_freq[u.event[0]] = u.class_proportion['no_matrix1']
#
#     else:
#         for u in v.deviant_obj:
#             if(u.pValue<0.05):
#                 deviant_elements_freq[u.event[0]] = u.class_proportion['no_matrix2']
#                 deviant_elements_f_score[u.event[0]] = np.min((u.f1_avg['f1C1'], u.f1_avg['f1C2']))
#
#                 # Storing elements only happened in log2
#             elif (u.class_proportion['proportion'] == None):
#                 one_class_element.append(u.event[0])
#                 deviant_elements_freq[u.event[0]] = u.class_proportion['no_matrix2']
#
#
#     #The total number of traces.
#     total_trace = np.sum([len(id[i].split(',')) for i in id.keys()])
#
#
#
#
#
#     #Iterating over traces and filter those that do not have any deviant element
#     iter=1
#     deviation_ratio_per_trace=0.1
#     trace_frequecny = 0.001
#     entropy_record=[]
#     while(iter==1):
#         log_filtered=[]
#         kgram = len(deviant_elements[0])
#         for j in range(len(log)):
#             trace=log[j]
#         #for trace in log:
#             count=0
#             for i in range(len(trace)-kgram+1):
#                 if ( ( tuple(trace[i:i+kgram]) in deviant_elements)):
#                     print float(len(id[j].split(",")))/ total_trace
#                     condition = ( (float(len(id[j].split(","))) / total_trace > trace_frequecny) and (deviant_elements_freq[tuple(trace[i:i+kgram])] >20) and(deviant_elements_f_score[tuple(trace[i:i+kgram])] >0.05 ) )
#                     if ( condition ):
#                         print "Iter",j
#                         print "Freq:",deviant_elements_freq[tuple(trace[i:i+kgram])] , tuple(trace[i:i+kgram])
#                         print "trace:", trace
#                         print "F-score:", deviant_elements_f_score[tuple(trace[i:i+kgram])]
#                         count+=1
#                         print "The ratio :",float(count)/len(trace)
#
#                         # Storing the information of deviant elements for plotting
#                         v.deviant_element_plot[tuple(trace[i:i + kgram])] = deviant_elements_freq[tuple(trace[i:i + kgram])]
#
#                         if(float(count)/len(trace)>deviation_ratio_per_trace):
#                             if (trace not in log_filtered):
#                                 log_filtered.append(trace)
#
#
#
#                 #Finding traces that contain elements that happened in exclusively one log
#                 elif ((tuple(trace[i:i + kgram]) in one_class_element)):
#                     condition = ((float(len(id[j].split(","))) / total_trace > 0.009) and (deviant_elements_freq[tuple(trace[i:i + kgram])] > 20))
#                     if(condition):
#                         print "One class Freq:", deviant_elements_freq[tuple(trace[i:i + kgram])], tuple(trace[i:i + kgram])
#                         if (trace not in log_filtered):
#                             log_filtered.append(trace)
#
#
#                 #Filtering log that contains elements that are deviant in terms of TimeStamp
#                 elif(tuple(trace[i:i + kgram]) in deviant_elements_timestamp):
#                     condition = ((float(len(id[j].split(","))) / total_trace > trace_frequecny))
#                     if (condition):
#                         count+=1
#                         print "TimeStamp Deviation!!", tuple(trace[i:i + kgram])
#                         if (float(count) / len(trace) > deviation_ratio_per_trace):
#                             v.deviant_element_timeStamp_plot[tuple(trace[i:i + kgram])] = 1
#                             if (trace not in log_filtered):
#                                 log_filtered.append(trace)
#
#
#                         #flag=1
#
#         #Computing the Entropy of the obtained graph
#         log= log_filtered
#         uniq_event =list(set([k for t in log_filtered for k in t]))
#
#         #Create adjacency matrix
#         adj_matrix = np.zeros((len(uniq_event),len(uniq_event)))
#
#         #filling the matrix
#         for t in log_filtered:
#             for i in range(len(t)-1):
#                 row_ind = uniq_event.index(t[i])
#                 col_ind = uniq_event.index(t[i+1])
#                 adj_matrix[row_ind][col_ind]=1
#         #From adj matrix the entropy of the graph is computed
#         model_entropy=stats.entropy(list(Counter(adj_matrix.flatten()).values()), base=2)
#         print "The adj matrix is:", adj_matrix
#         print "The entropy is:", model_entropy
#         print "The Deviation ratio:",deviation_ratio_per_trace
#
#         if(model_entropy <0.9):
#             iter =0
#         else:
#             #Increasing the entropy slowly slowly by looking at its standard deviation (Adaptive increasing)
#             entropy_record.append(model_entropy)
#             deviation_ratio_per_trace+=0.01 + 0.001*1/(1+np.std(entropy_record))
#
#             #Increasing the frequency threshold to include traces
#             trace_frequecny +=0.001
#
#
#
#     #Calling TS_plot()
#     print "entropy record:",entropy_record
#     print "the log filtered:", log_filtered
#     #-----------------------------------------------------------------
#
#     TS_plot(log_filtered,log_filtered)
#
#     return log_filtered, v.deviant_element_plot, deviant_elements_freq, v.deviant_element_timeStamp_plot


###########################################################################################
#Plotting Transition systems for variants
import Main as m
def Deviant_Plot_ControlFlow(log1,log2):
    m.Print_Screen_Off('off')

    #Obtaining the required elements
    ## deviant_element_plot1 = :{('Admission#NC', 'Admission#NC'): 118,
                                 # ('Admission#NC', 'Leucocytes'): 253,
                                 # ('CRP', 'Admission#NC'): 120,
                                 # ('CRP', 'CRP'): 152,
                                 # ('CRP', 'LacticAcid'): 282,
                                 # ('LacticAcid', 'CRP'): 190,
                                 # ('LacticAcid', 'IV#Liquid'): 88,
                                 # ('LacticAcid', 'Leucocytes'): 241,
                                 # ('Leucocytes', 'CRP'): 503,
                                 # ('Leucocytes', 'LacticAcid'): 203,
                                 # ('Leucocytes', 'Leucocytes'): 159,
                                 # ('Leucocytes', 'Release#A'): 149,
                                 # ('Release#A', 'Return#ER'): 189}


    log_filtered1 , deviant_element_plot1, deviant_elements_freq1 = Log_Filter_ControlFlow(v.log1, v.id1)
    log_filtered2 , deviant_element_plot2, deviant_elements_freq2 = Log_Filter_ControlFlow(v.log2, v.id2)
    m.Print_Screen_Off('on')

    print "Deviant part from log1:", pprint.pprint(deviant_element_plot1)
    print "Deviant part from log2:", pprint.pprint(deviant_element_plot2)

    #Finding common edges between to plots
    edge1=[]
    for t in log_filtered1:
        for i in range(len(t) - 1):
            if ([t[i], t[i + 1]] not in edge1):
                print "common edge:", edge1
                edge1.append((t[i], t[i + 1]))
            else:
                edge1.append((t[i], t[i + 1]))

    edge2=[]
    for t in log_filtered2:
        for i in range(len(t) - 1):
            if ([t[i], t[i + 1]] not in edge2):
                print "common edge:", edge2
                edge2.append((t[i], t[i + 1]))
            else:
                edge2.append((t[i], t[i + 1]))

    common_edge = list( set(edge1).intersection(set(edge2)) )





    ############################# Plotting##########################################
    # Plotting the first log
    log = log_filtered1
    uniq_event = list(set([k for t in log for k in t]))

    # Create adjuncy matrix
    adj_matrix = np.zeros((len(uniq_event), len(uniq_event)))

    # filling the matrix
    for t in log:
        for i in range(len(t) - 1):
            row_ind = uniq_event.index(t[i])
            col_ind = uniq_event.index(t[i + 1])
            adj_matrix[row_ind][col_ind] = 1

    # Creating a directed graph
    dot = Digraph(comment='TS system')

    nodes = []
    for e in uniq_event:
        nodes.append(dot.node(e))

    common_deviant_elements = set(deviant_element_plot2.keys()).intersection(set(deviant_element_plot1.keys()))

    # adding edges
    temp = []
    for t in log:
        for i in range(len(t) - 1):
            if ([t[i], t[i + 1]] not in temp):
                print "Temmp:", temp
                #if(tuple([t[i], t[i + 1]]) not in deviant_element_plot1.keys() + deviant_element_plot2.keys()):
                if (tuple([t[i], t[i + 1]]) not in common_deviant_elements):
                #if( tuple([t[i], t[i + 1]]) in common_edge):
                    #if ( ( tuple([t[i], t[i + 1]]) in deviant_element_plot1 ) or (tuple([t[i], t[i + 1]]) in deviant_element_plot2) ):
                    temp.append([t[i], t[i + 1]])
                    dot.edge(t[i], t[i + 1],style='dashed' )
                # else:
                #     temp.append([t[i], t[i + 1]])
                #     dot.edge(t[i], t[i + 1], style='dashed')

                # else:
                #     temp.append([t[i], t[i + 1]])
                #     dot.edge(t[i], t[i + 1],style='dashed')

    #Superimposing other deviations
    for key in common_deviant_elements:
    #for key in set(deviant_element_plot1.keys() + deviant_element_plot2.keys()):
         dot.edge(key[0], key[1], color='red', label=str(deviant_elements_freq1[tuple([key[0], key[1]])]), style='bold')


    # Graphviz must be installed in order to save output in pdf format
    if platform.system() == "Windows":
        os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
    elif platform.system() == "Linux":
        if ('graphviz' not in sys.modules):
            os.system('sudo apt-get install graphviz')

    dot.render('test-output/' + time.strftime("%m/%d/%Y, %H%M%S")  + "log1"+".pdf", view=True)
    # -----------------------------------------------------------------

    # Plotting the second log
    log = log_filtered2
    uniq_event = list(set([k for t in log for k in t]))

    # Create adjuncy matrix
    adj_matrix = np.zeros((len(uniq_event), len(uniq_event)))

    # filling the matrix
    for t in log:
        for i in range(len(t) - 1):
            row_ind = uniq_event.index(t[i])
            col_ind = uniq_event.index(t[i + 1])
            adj_matrix[row_ind][col_ind] = 1

    # Creating a directed graph
    dot = Digraph(comment='TS system')

    nodes = []
    for e in uniq_event:
        nodes.append(dot.node(e))




    # adding edges
    temp = []
    for t in log:
        for i in range(len(t) - 1):
            if ([t[i], t[i + 1]] not in temp):
                print "Temmp:", temp
                #if(tuple([t[i], t[i + 1]]) not in deviant_element_plot1.keys() + deviant_element_plot2.keys()):
                if (tuple([t[i], t[i + 1]]) not in common_deviant_elements):
                #if( tuple([t[i], t[i + 1]]) in common_edge):
                    #if ( ( tuple([t[i], t[i + 1]]) in deviant_element_plot1 ) or (tuple([t[i], t[i + 1]]) in deviant_element_plot2) ):
                    temp.append([t[i], t[i + 1]])
                    dot.edge(t[i], t[i + 1],style='dashed' )
                # else:
                #     temp.append([t[i], t[i + 1]])
                #     dot.edge(t[i], t[i + 1], style='dashed')


    # # adding edges
    # temp = []
    # for t in log:
    #     for i in range(len(t) - 1):
    #         if ([t[i], t[i + 1]] not in temp):
    #             print "Temmp:", temp
    #             if( tuple([t[i], t[i + 1]]) in common_edge):
    #                 if ( ( tuple([t[i], t[i + 1]]) in deviant_element_plot1 ) or (tuple([t[i], t[i + 1]]) in deviant_element_plot2) ):
    #                     temp.append([t[i], t[i + 1]])
    #                     dot.edge(t[i], t[i + 1], color='red', label= str( deviant_elements_freq2[tuple([t[i], t[i + 1]])] ) ,style='bold')
    #                 else:
    #                     temp.append([t[i], t[i + 1]])
    #                     dot.edge(t[i], t[i + 1],style='dashed')
    #
    #             else:
    #                 temp.append([t[i], t[i + 1]])
    #                 dot.edge(t[i], t[i + 1],style='dashed')
    #
    # # Superimposing other deviations
    for key in common_deviant_elements:
    #for key in set(deviant_element_plot1.keys() + deviant_element_plot2.keys()):
            dot.edge(key[0], key[1], color='red', label=str(deviant_elements_freq2[tuple([key[0], key[1]])]),
                     style='bold')

    # Graphviz must be installed in order to save output in pdf format
    if platform.system() == "Windows":
        os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
    elif platform.system() == "Linux":
        if ('graphviz' not in sys.modules):
            os.system('sudo apt-get install graphviz')

    dot.render('test-output/' + time.strftime("%m/%d/%Y, %H%M%S") + "log2"+".pdf", view=True)

    # -----------------------------------------------------------------

###########################################################################################
#Plotting Transition systems for variants +Duration time
def Deviant_Plot_ControlFlow_Time(log1,log2):
    m.Print_Screen_Off('off')

    #Obtaining the required elements
    ## deviant_element_plot1 = :{('Admission#NC', 'Admission#NC'): 118,
                                 # ('Admission#NC', 'Leucocytes'): 253,
                                 # ('CRP', 'Admission#NC'): 120,
                                 # ('CRP', 'CRP'): 152,
                                 # ('CRP', 'LacticAcid'): 282,
                                 # ('LacticAcid', 'CRP'): 190,
                                 # ('LacticAcid', 'IV#Liquid'): 88,
                                 # ('LacticAcid', 'Leucocytes'): 241,
                                 # ('Leucocytes', 'CRP'): 503,
                                 # ('Leucocytes', 'LacticAcid'): 203,
                                 # ('Leucocytes', 'Leucocytes'): 159,
                                 # ('Leucocytes', 'Release#A'): 149,
                                 # ('Release#A', 'Return#ER'): 189}

    #Getting elements that are statistically significant in terms of duration time
    deviant_element_duration_time = [u.event[0] for u in v.deviant_obj_timestamp]


    log_filtered1 , deviant_element_plot1, deviant_elements_freq1 = Log_Filter_ControlFlow(v.log1, v.id1)
    log_filtered2 , deviant_element_plot2, deviant_elements_freq2 = Log_Filter_ControlFlow(v.log2, v.id2)
    m.Print_Screen_Off('on')

    print "Deviant part from log1:", pprint.pprint(deviant_element_plot1)
    print "Deviant part from log2:", pprint.pprint(deviant_element_plot2)

    #Finding common edges between to plots
    edge1=[]
    for t in log_filtered1:
        for i in range(len(t) - 1):
            if ([t[i], t[i + 1]] not in edge1):
                print "common edge:", edge1
                edge1.append((t[i], t[i + 1]))
            else:
                edge1.append((t[i], t[i + 1]))

    edge2=[]
    for t in log_filtered2:
        for i in range(len(t) - 1):
            if ([t[i], t[i + 1]] not in edge2):
                print "common edge:", edge2
                edge2.append((t[i], t[i + 1]))
            else:
                edge2.append((t[i], t[i + 1]))

    common_edge = list( set(edge1).intersection(set(edge2)) )





    ############################# Plotting##########################################
    # Plotting the first log
    log = log_filtered1
    uniq_event = list(set([k for t in log for k in t]))

    # Create adjuncy matrix
    adj_matrix = np.zeros((len(uniq_event), len(uniq_event)))

    # filling the matrix
    for t in log:
        for i in range(len(t) - 1):
            row_ind = uniq_event.index(t[i])
            col_ind = uniq_event.index(t[i + 1])
            adj_matrix[row_ind][col_ind] = 1

    # Creating a directed graph
    dot = Digraph(comment='TS system')

    nodes = []
    for e in uniq_event:
        nodes.append(dot.node(e,shape="rect"))

    common_deviant_elements = set(deviant_element_plot2.keys()).intersection(set(deviant_element_plot1.keys()))

    # adding edges
    temp = []
    for t in log:
        for i in range(len(t) - 1):
            if ([t[i], t[i + 1]] not in temp):
                print "Temmp:", temp
                #if(tuple([t[i], t[i + 1]]) not in deviant_element_plot1.keys() + deviant_element_plot2.keys()):
                if (tuple([t[i], t[i + 1]]) not in common_deviant_elements):
                    if (tuple([t[i], t[i + 1]]) in deviant_element_duration_time):
                #if( tuple([t[i], t[i + 1]]) in common_edge):
                    #if ( ( tuple([t[i], t[i + 1]]) in deviant_element_plot1 ) or (tuple([t[i], t[i + 1]]) in deviant_element_plot2) ):
                        temp.append([t[i], t[i + 1]])
                        dot.edge(t[i], t[i + 1],style='dashed', color = 'black' )
                    else:
                        temp.append([t[i], t[i + 1]])
                        dot.edge(t[i], t[i + 1], style='solid')


    #Superimposing other deviations
    for key in common_deviant_elements:
        if (key in deviant_element_duration_time):
            dot.edge(key[0], key[1], color='red', label=str(deviant_elements_freq1[tuple([key[0], key[1]])]), style='dashed')
        else:
            dot.edge(key[0], key[1], color='red', label=str(deviant_elements_freq1[tuple([key[0], key[1]])]),
                     style='bold')


    # Graphviz must be installed in order to save output in pdf format
    if platform.system() == "Windows":
        os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
    elif platform.system() == "Linux":
        if ('graphviz' not in sys.modules):
            os.system('sudo apt-get install graphviz')

    dot.render('test-output(log1)/' + time.strftime("%m/%d/%Y, %H%M%S")  + "log1"+".pdf", view=True)
    # -----------------------------------------------------------------

    # Plotting the second log
    log = log_filtered2
    uniq_event = list(set([k for t in log for k in t]))

    # Create adjuncy matrix
    adj_matrix = np.zeros((len(uniq_event), len(uniq_event)))

    # filling the matrix
    for t in log:
        for i in range(len(t) - 1):
            row_ind = uniq_event.index(t[i])
            col_ind = uniq_event.index(t[i + 1])
            adj_matrix[row_ind][col_ind] = 1

    # Creating a directed graph
    dot = Digraph(comment='TS system')

    nodes = []
    for e in uniq_event:
        nodes.append(dot.node(e,shape="rect"))




    # adding edges
    temp = []
    for t in log:
        for i in range(len(t) - 1):
            if ([t[i], t[i + 1]] not in temp):
                print "Temmp:", temp
                #if(tuple([t[i], t[i + 1]]) not in deviant_element_plot1.keys() + deviant_element_plot2.keys()):
                if (tuple([t[i], t[i + 1]]) not in common_deviant_elements):
                    if (tuple([t[i], t[i + 1]]) in deviant_element_duration_time):
                #if( tuple([t[i], t[i + 1]]) in common_edge):
                    #if ( ( tuple([t[i], t[i + 1]]) in deviant_element_plot1 ) or (tuple([t[i], t[i + 1]]) in deviant_element_plot2) ):
                        temp.append([t[i], t[i + 1]])
                        dot.edge(t[i], t[i + 1],style='dashed', color = 'black' )
                    else:
                        temp.append([t[i], t[i + 1]])
                        dot.edge(t[i], t[i + 1], style='solid')



    # # adding edges
    # temp = []
    # for t in log:
    #     for i in range(len(t) - 1):
    #         if ([t[i], t[i + 1]] not in temp):
    #             print "Temmp:", temp
    #             if( tuple([t[i], t[i + 1]]) in common_edge):
    #                 if ( ( tuple([t[i], t[i + 1]]) in deviant_element_plot1 ) or (tuple([t[i], t[i + 1]]) in deviant_element_plot2) ):
    #                     temp.append([t[i], t[i + 1]])
    #                     dot.edge(t[i], t[i + 1], color='red', label= str( deviant_elements_freq2[tuple([t[i], t[i + 1]])] ) ,style='bold')
    #                 else:
    #                     temp.append([t[i], t[i + 1]])
    #                     dot.edge(t[i], t[i + 1],style='dashed')
    #
    #             else:
    #                 temp.append([t[i], t[i + 1]])
    #                 dot.edge(t[i], t[i + 1],style='dashed')
    #
    # # Superimposing other deviations
    for key in common_deviant_elements:
        if (key in deviant_element_duration_time):
            dot.edge(key[0], key[1], color='red', label=str(deviant_elements_freq2[tuple([key[0], key[1]])]),
                     style='dashed')
        else:
            dot.edge(key[0], key[1], color='red', label=str(deviant_elements_freq2[tuple([key[0], key[1]])]),
                     style='bold')

    # Graphviz must be installed in order to save output in pdf format
    if platform.system() == "Windows":
        os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
    elif platform.system() == "Linux":
        if ('graphviz' not in sys.modules):
            os.system('sudo apt-get install graphviz')

    dot.render('test-output(log2)/' + time.strftime("%m/%d/%Y, %H%M%S") + "log2"+".pdf", view=True)

    # -----------------------------------------------------------------






##########################################################################################
############################################################################################
# #Plotting Transition systems for variants
# import Main as m
# def Deviant_Plot_ControlFlow_TimeStamp(log1,log2):
#     m.Print_Screen_Off('off')
#
#     #Obtaining the required elements
#     ## deviant_element_plot1 = :{('Admission#NC', 'Admission#NC'): 118,
#                                  # ('Admission#NC', 'Leucocytes'): 253,
#                                  # ('CRP', 'Admission#NC'): 120,
#                                  # ('CRP', 'CRP'): 152,
#                                  # ('CRP', 'LacticAcid'): 282,
#                                  # ('LacticAcid', 'CRP'): 190,
#                                  # ('LacticAcid', 'IV#Liquid'): 88,
#                                  # ('LacticAcid', 'Leucocytes'): 241,
#                                  # ('Leucocytes', 'CRP'): 503,
#                                  # ('Leucocytes', 'LacticAcid'): 203,
#                                  # ('Leucocytes', 'Leucocytes'): 159,
#                                  # ('Leucocytes', 'Release#A'): 149,
#                                  # ('Release#A', 'Return#ER'): 189}
#
#
#     log_filtered1 , deviant_element_plot1, deviant_elements_freq1, deviant_elements_timestamp_plot = Log_Filter_ControlFlow_TimeStamp(v.log1, v.id1)
#     log_filtered2 , deviant_element_plot2, deviant_elements_freq2, deviant_elements_timestamp_plot  = Log_Filter_ControlFlow_TimeStamp(v.log2, v.id2)
#     m.Print_Screen_Off('on')
#
#     print "Deviant part from log1:", pprint.pprint(deviant_element_plot1)
#     print "Deviant part from log2:", pprint.pprint(deviant_element_plot2)
#     print "Deviant parts Timestamp:", pprint.pprint(deviant_elements_timestamp_plot)
#
#     #Finding common edges between to plots
#     edge1=[]
#     for t in log_filtered1:
#         for i in range(len(t) - 1):
#             if ([t[i], t[i + 1]] not in edge1):
#                 print "common edge:", edge1
#                 edge1.append((t[i], t[i + 1]))
#             else:
#                 edge1.append((t[i], t[i + 1]))
#
#     edge2=[]
#     for t in log_filtered2:
#         for i in range(len(t) - 1):
#             if ([t[i], t[i + 1]] not in edge2):
#                 print "common edge:", edge2
#                 edge2.append((t[i], t[i + 1]))
#             else:
#                 edge2.append((t[i], t[i + 1]))
#
#     common_edge = list( set(edge1).intersection(set(edge2)) )
#
#
#
#
#
#     ############################# Plotting##########################################
#     # Plotting the first log
#     log = log_filtered1
#     uniq_event = list(set([k for t in log for k in t]))
#
#     # Create adjuncy matrix
#     adj_matrix = np.zeros((len(uniq_event), len(uniq_event)))
#
#     # filling the matrix
#     for t in log:
#         for i in range(len(t) - 1):
#             row_ind = uniq_event.index(t[i])
#             col_ind = uniq_event.index(t[i + 1])
#             adj_matrix[row_ind][col_ind] = 1
#
#     # Creating a directed graph
#     dot = Digraph(comment='TS system')
#
#     nodes = []
#     for e in uniq_event:
#         nodes.append(dot.node(e))
#
#     common_deviant_elements = set(deviant_element_plot2.keys()).intersection(set(deviant_element_plot1.keys()))
#
#     # adding edges
#     temp = []
#     for t in log:
#         for i in range(len(t) - 1):
#             if ([t[i], t[i + 1]] not in temp):
#                 print "Temmp:", temp
#                 #if(tuple([t[i], t[i + 1]]) not in deviant_element_plot1.keys() + deviant_element_plot2.keys()):
#                 if ((tuple([t[i], t[i + 1]]) not in ( common_deviant_elements.union( set( deviant_elements_timestamp_plot.keys() ) ) ) ) ):
#                 #if( tuple([t[i], t[i + 1]]) in common_edge):
#                     #if ( ( tuple([t[i], t[i + 1]]) in deviant_element_plot1 ) or (tuple([t[i], t[i + 1]]) in deviant_element_plot2) ):
#                     temp.append([t[i], t[i + 1]])
#                     dot.edge(t[i], t[i + 1],style='dashed' )
#                 # else:
#                 #     temp.append([t[i], t[i + 1]])
#                 #     dot.edge(t[i], t[i + 1], style='dashed')
#
#                 # else:
#                 #     temp.append([t[i], t[i + 1]])
#                 #     dot.edge(t[i], t[i + 1],style='dashed')
#
#     #Superimposing other deviations
#     for key in common_deviant_elements:
#
#     #for key in set(deviant_element_plot1.keys() + deviant_element_plot2.keys()):
#         if(key in common_deviant_elements.intersection(set(deviant_elements_timestamp_plot.keys()))):
#             dot.edge(key[0], key[1], color='red', label=str(deviant_elements_freq1[tuple([key[0], key[1]])]), style='dashed')
#
#         elif(key in common_deviant_elements.difference(set(deviant_elements_timestamp_plot.keys()))):
#             dot.edge(key[0], key[1], color='red', label=str(deviant_elements_freq1[tuple([key[0], key[1]])]), style='bold')
#
#         elif (key in set(deviant_elements_timestamp_plot.keys()).difference(common_deviant_elements)):
#             dot.edge(key[0], key[1], color='green', style='bold')
#
#
#
#     # Graphviz must be installed in order to save output in pdf format
#     if platform.system() == "Windows":
#         os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
#     elif platform.system() == "Linux":
#         if ('graphviz' not in sys.modules):
#             os.system('sudo apt-get install graphviz')
#
#     dot.render('test-output/' + time.strftime("%m/%d/%Y, %H%M%S")  + "log1"+".pdf", view=True)
#     # -----------------------------------------------------------------
#
#     # Plotting the second log
#     log = log_filtered2
#     uniq_event = list(set([k for t in log for k in t]))
#
#     # Create adjuncy matrix
#     adj_matrix = np.zeros((len(uniq_event), len(uniq_event)))
#
#     # filling the matrix
#     for t in log:
#         for i in range(len(t) - 1):
#             row_ind = uniq_event.index(t[i])
#             col_ind = uniq_event.index(t[i + 1])
#             adj_matrix[row_ind][col_ind] = 1
#
#     # Creating a directed graph
#     dot = Digraph(comment='TS system')
#
#     nodes = []
#     for e in uniq_event:
#         nodes.append(dot.node(e))
#
#
#
#
#     # adding edges
#     temp = []
#     for t in log:
#         for i in range(len(t) - 1):
#             if ([t[i], t[i + 1]] not in temp):
#                 print "Temmp:", temp
#                 #if(tuple([t[i], t[i + 1]]) not in deviant_element_plot1.keys() + deviant_element_plot2.keys()):
#                 if ((tuple([t[i], t[i + 1]]) not in ( common_deviant_elements.union( set( deviant_elements_timestamp_plot.keys() ) ) ) ) ):
#                 #if( tuple([t[i], t[i + 1]]) in common_edge):
#                     #if ( ( tuple([t[i], t[i + 1]]) in deviant_element_plot1 ) or (tuple([t[i], t[i + 1]]) in deviant_element_plot2) ):
#                     temp.append([t[i], t[i + 1]])
#                     dot.edge(t[i], t[i + 1],style='dashed' )
#                 # else:
#                 #     temp.append([t[i], t[i + 1]])
#                 #     dot.edge(t[i], t[i + 1], style='dashed')
#
#
#     # # adding edges
#     # temp = []
#     # for t in log:
#     #     for i in range(len(t) - 1):
#     #         if ([t[i], t[i + 1]] not in temp):
#     #             print "Temmp:", temp
#     #             if( tuple([t[i], t[i + 1]]) in common_edge):
#     #                 if ( ( tuple([t[i], t[i + 1]]) in deviant_element_plot1 ) or (tuple([t[i], t[i + 1]]) in deviant_element_plot2) ):
#     #                     temp.append([t[i], t[i + 1]])
#     #                     dot.edge(t[i], t[i + 1], color='red', label= str( deviant_elements_freq2[tuple([t[i], t[i + 1]])] ) ,style='bold')
#     #                 else:
#     #                     temp.append([t[i], t[i + 1]])
#     #                     dot.edge(t[i], t[i + 1],style='dashed')
#     #
#     #             else:
#     #                 temp.append([t[i], t[i + 1]])
#     #                 dot.edge(t[i], t[i + 1],style='dashed')
#     #
#     #Superimposing other deviations
#     for key in common_deviant_elements:
#
#     #for key in set(deviant_element_plot1.keys() + deviant_element_plot2.keys()):
#         if(key in common_deviant_elements.intersection(set(deviant_elements_timestamp_plot.keys()))):
#             dot.edge(key[0], key[1], color='red', label=str(deviant_elements_freq2[tuple([key[0], key[1]])]), style='dashed')
#
#         elif(key in common_deviant_elements.difference(set(deviant_elements_timestamp_plot.keys()))):
#             dot.edge(key[0], key[1], color='red', label=str(deviant_elements_freq2[tuple([key[0], key[1]])]), style='bold')
#
#         elif (key in set(deviant_elements_timestamp_plot.keys()).difference(common_deviant_elements)):
#             dot.edge(key[0], key[1], color='green', style='bold')
#
#
#
#     # Graphviz must be installed in order to save output in pdf format
#     if platform.system() == "Windows":
#         os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
#     elif platform.system() == "Linux":
#         if ('graphviz' not in sys.modules):
#             os.system('sudo apt-get install graphviz')
#
#     dot.render('test-output/' + time.strftime("%m/%d/%Y, %H%M%S") + "log2"+".pdf", view=True)
#
#     # -----------------------------------------------------------------

##############################################################################################



#This is a temporary function (for plotting histograms)
def Event_Position_Plot(log, id, log2 ,id2, event = ('Payment', 'Payment')):
    '''

    :param log: the log input, imported by Readin_log_all()
    :param id:  the log id
    :param event:
    :return:
    '''
    temp=[]
    for i in id:
        for j in range(len(log[i])):
            if (log[i][j] ==  ('Payment', 'Send#for#Credit#Collection')):
                temp+=[j]*len(id[i].split(','))




    temp2=[]
    for i in id2:
        for j in range(len(log2[i])):
            if (log2[i][j] ==  ('Payment', 'Send#for#Credit#Collection')):
                temp2+=[j]*len(id2[i].split(','))


    plt.style.use('ggplot')
    #plt.hist(temp, bins=20,  color="steelblue", label="Fine >=50")
    plt.hist(temp2, bins=20,  color="gray", label="Fine<50")
    # sns.distplot(temp, norm_hist=True,  color="skyblue", label="Sepal Length")
    # sns.distplot(temp2, norm_hist= True, color="gray", label="Sepal Length")
    plt.xticks(np.arange(20))
    plt.xlabel('Position of (Payment, Send for Credit Collection)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()
    return temp


























