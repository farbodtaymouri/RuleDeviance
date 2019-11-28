import Preparation as pr
import operator
import itertools
import numpy as np
import pandas as pd
import os
import sys
import time
import tqdm
import pprint
import Reading_Wriitng as rw
import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt
import DimReduction as dr
import Visualization as vs
import Variables as v
import Deviation as dv



def Print_Screen_Off(exp):
    if exp =='off':
        # -----------This psrt of code is to prevent printing on the screen!!!--------
        if (v.sys_stdout_old == None):
            v.sys_stdout_old =  sys.stdout

        sys.stdout = open(os.devnull, "w")
        f = open(os.devnull, 'w')
        sys.stdout = f
    else:
        sys.stdout = v.sys_stdout_old
# ----------------------------------------------------------------------------


def Main(log_path, destination):

    #Print_Screen_Off('on')


    #List of lists as input
    '''log = [['a', 'b', 'c', 'a', 'a', 'c'], ['a', 'b', 'b', 'a', 'a', 'c'],
           ['a', 'b', 'b', 'a', 'c', 'c']]'''

    log = [['d','c','e','f'], ['c','e','d','f'], ['e','d','f','c'], ['e','d','c','f'],  ['e','c','d','f'],['c','f'] ]

    id = {0: '0', 1:'1', 2:'2', 3:'3'}


    ##log, id = rw.Reading_log_all(log_path, destination)


    #Normalizing the log (equating the length)
    obj1 = pr.Seq2Vec()
    log = obj1.Log_Normalization(log)


    matrix=[]

    for trace in log:
        obj =pr.Seq2Vec()    #Creating an object for the trace
        obj.unique_Events_log = obj1.unique_Events_log   #Giving the unique name of event log
        obj.Start(trace)     #Starting Sequence to vector
        #Flattening the dictionary of coefficients to make a matrix row
        print "the cf:",obj.dic_WaveletCfVector.values()
        matrix.append(reduce(operator.concat, obj.dic_WaveletCfVector.values()) )

    #print "The matrix is:", matrix
    print np.matrix(matrix).shape
    plt.figure(figsize=(10, 5))
    matrix = np.array(matrix)
    dend = shc.dendrogram(shc.linkage(matrix, method='ward'), p=30, truncate_mode = 'lastp')
    plt.savefig('dendogram.png')

    # create clusters
    '''hc = AgglomerativeClustering(n_clusters=2, affinity='euclidean', linkage='ward')
    y= hc.fit_predict(matrix)
    for row in range(len(matrix)):
        print "id:", id[row], ", cluster:", y[row]'''

    #from sklearn.decomposition import PCA
    #pca = PCA(n_components=2)
    #principalComponents = pca.fit_transform(db)

    return matrix
#############################################################################

###This has been created on 25/08/2019
def Main2(log_path, destination):

    #Print_Screen_Off('on')


    #List of lists as input
    '''log = [['a', 'b', 'c', 'a', 'a', 'c'], ['a', 'b', 'b', 'a', 'a', 'c'],
           ['a', 'b', 'b', 'a', 'c', 'c']]'''

    # log2 = [['d','c','e','f'], ['c','e','d','f'], ['c','d','f','e'], ['c','d','e','f'],  ['c','e','d','f'],['c','f'] ]*10
    #
    # #log2 = [['f','c','d','e'], ['c','f']]
    # log = [['c', 'e', 'd', 'f'],['d', 'c', 'e', 'f'], ['c', 'd', 'e', 'f'] , ['c', 'd', 'f', 'e'], ['c', 'e', 'd', 'f']]*10

    log = [['a','b']]*6 + [['b','a','c']]*30
    log2 = [['a', 'c']] * 6 + [['b', 'a', 'c']] * 30

    # log =[['a','b','c', 'a','b']]*3
    # log2=[['b','a'], ['a','b','c']]*15



    # log, id = rw.Reading_log_all(log_path, destination)
    # log = log[0:len(log)/2]
    # log2 = log[(len(log)/2 +1):]

    # log, id = rw.Reading_log_all('C:/Users/taymouri/Desktop/model temp/trash2/BPIC15_1.xes', destination)
    # log2, id = rw.Reading_log_all('C:/Users/taymouri/Desktop/model temp/trash2/BPIC15_2.xes', destination)
    # log=log[1:50]
    # log2 = log2[1:50]





    #Normalizing the log (equating the length)
    obj1 = pr.Seq2Vec()
    log_tot = obj1.Log_Normalization(log + log2)


    matrix=[]
    matrix_log2 = []

    #for trace in log:
    for i in range(len(log_tot)):
        if (i <= (len(log)-1)):
            trace = log_tot[i]
            obj =pr.Seq2Vec()    #Creating an object for the trace
            obj.unique_Events_log = obj1.unique_Events_log   #Giving the unique name of event log
            obj.Start(trace)     #Starting Sequence to vector
            #Flattening the dictionary of coefficients to make a matrix row
            print "the cf:",obj.dic_WaveletCfVector.values()
            matrix.append(reduce(operator.concat, obj.dic_WaveletCfVector.values()) )
            #del obj

        else:
            trace = log_tot[i]
            obj = pr.Seq2Vec()  # Creating an object for the trace
            obj.unique_Events_log = obj1.unique_Events_log  # Giving the unique name of event log
            obj.Start(trace)  # Starting Sequence to vector
            # Flattening the dictionary of coefficients to make a matrix row
            print "the cf:", obj.dic_WaveletCfVector.values()
            matrix_log2.append(reduce(operator.concat, obj.dic_WaveletCfVector.values()))
            #del obj











    #print "The matrix is:", matrix
    print "Design matrix for variants 1:", np.matrix(matrix).shape
    plt.figure(figsize=(10, 5))
    matrix = np.array(matrix)
    dend = shc.dendrogram(shc.linkage(matrix, method='ward'), p=30, truncate_mode = 'lastp')
    plt.savefig('dendogram.png')

    # print "The matrix is:", matrix
    print "Design matrix for variants 2:", np.matrix(matrix_log2).shape
    plt.figure(figsize=(10, 5))
    matrix_log2 = np.array(matrix_log2)
    dend = shc.dendrogram(shc.linkage(matrix_log2, method='ward'), p=30, truncate_mode='lastp')
    plt.savefig('dendogram.png')



    #-------------
    print "Dimensionality_Reduction(m) has been started"
    # matrix = matrix[:, range(0, 5) + range(10, 15)]
    # matrix_log2 = matrix_log2[:, range(0, 5) + range(10, 15)]
    print "log1 matrix:\n", matrix
    print "log2 matrix:\n", matrix_log2
    dr.Ttset_Variants(matrix, matrix_log2)
    #-----------

    #Selecting the best feature to discriminate the classes
    dr.Feature_Selection(matrix,matrix_log2, obj)

    return matrix, obj

############################################################################
#This version of Main is created on 11/09/2109
def Main3(log_path,d):

    Print_Screen_Off('off')


    # log1 = [['a', 'e','f','d']] * 10 + [['b', 'a', 'c']] * 40
    # log2 = [['c', 'b','a','f','a']] * 10 + [['b', 'a', 'c']] *10

    log1 = [['c', 'b','a','f','a']] * 20  + [['b', 'a', 'c']] * 10
    log2 = [['a', 'e','f','d']] * 10+ [['b', 'a', 'c']] *10


    log1, id1 = dv.Log_dic(log1)
    log2, id2 = dv.Log_dic(log2)



    #Plotlog
    #viz.TS_plot(log1,log2)

    # log =[['a','b','c', 'a','b']]*3
    # log2=[['b','a'], ['a','b','c']]*15

    # log, id = rw.Reading_log_all(log_path, destination)
    # log = log[0:len(log)/2]
    # log2 = log[(len(log)/2 +1):]

    path1= 'C:/Users/taymouri/Desktop/model temp/trash2/Road_teraffic_fineamount_Above and equal50.xes'
    path2 ='C:/Users/taymouri/Desktop/model temp/trash2/Road_teraffic_fineamount_Below50.xes'

    # path1= 'C:/Users/taymouri/Desktop/model temp/trash2/Sepsis(age above 70).xes'
    # path2= 'C:/Users/taymouri/Desktop/model temp/trash2/Sepsis(age under 35).xes'

    # path1= 'C:/Users/taymouri/Desktop/model temp/trash2/logs/BPIC13_incidents_orgline_A2.xes'
    # path2= 'C:/Users/taymouri/Desktop/model temp/trash2/logs/BPIC13_incidents_orgline_C.xes'

    # path1 = 'C:/Users/taymouri/Desktop/model temp/trash2/logs/BPIC15_1.xes'
    # path2 = 'C:/Users/taymouri/Desktop/model temp/trash2/logs/BPIC15_2.xes'

    #
    # path1= 'C:/Users/taymouri/Desktop/model temp/trash2/BPI_2013_OrgCountry_br.xes'
    # path2= 'C:/Users/taymouri/Desktop/model temp/trash2/BPI_2013_OrgCountry_US.xes'


    log1, id1 = rw.Reading_log_all(path1, d)
    log2, id2 = rw.Reading_log_all(path2, d)

    # log1, id1 = rw.Reading_log_all('C:/Users/taymouri/Desktop/model temp/trash2/Sepsis(age above 70).xes', d)
    # log2, id2 = rw.Reading_log_all('C:/Users/taymouri/Desktop/model temp/trash2/Sepsis(age under 35).xes', d)

    # log1 = log1[0:100]
    # log2 = log2[0:100]


    #Storing the id of unique traces and unique logs
    v.id1 = id1
    v.id2 = id2
    v.log1=log1[:]
    v.log2=log2[:]

    # #Log analysing
    # log1,log2, id1, id2 = dr.Log_Stat()
    # #return


    #return
    #Computing Ngram
    kgram=2
    log1 = dr.Ngram_Compute(log1,kgram)
    log2 = dr.Ngram_Compute(log2,kgram)
    print "Bigram1:", pprint.pprint(log1)
    print "Bigram2:",pprint.pprint(log2)


    # #--------------------------------
    #The method that reads log, only extracts unique ones and the number of occurrence, therefore, we take it back to the normal mode
    # temp= [[log1[i]] * len(id1[i].split(',')) for i in range(len(log1))]
    # log1 =  [case for l in temp for case in l]
    #
    # temp = [[log2[i]] * len(id2[i].split(',')) for i in range(len(log2))]
    # log2 = [case for l in temp for case in l]

    #------------------------------


    #Creating a dictionary of feature-value. See the definition of the function
    log_featureValue_dic, log2_featureValue_dic = dr.FeatureValue_Representation(log1, log2)

    # Ranking the features
    sorted_events = dr.Feature_Ranking(log_featureValue_dic, log2_featureValue_dic) #[('d', 0.178), ('f', 0.17), ('e', 0.0)]

    # Applying T-test to the sorted events
    ttest_pvalues={}
    v.stat_test_result = {}
    v.deviant_obj=[]
    for i in tqdm.tqdm(range(len(sorted_events))):
        event = sorted_events[i]
    #for event in sorted_events:
        event_matrix1, event_matrix2 = dr.Event_Matrxi_Creation([event[0]], log_featureValue_dic,log2_featureValue_dic)
        #dr.LDA_Compute([event[0]],event_matrix1,event_matrix2)

        #Creating an object for each event to save corresponding information
        evObj = dv.Deviant_Element([event[0]])
        v.deviant_obj.append(evObj)
        dr.SVM_Classifier(evObj, event_matrix1, event_matrix2)
        # dr.SVM_Classifier([event[0]],event_matrix1,event_matrix2)

        ##dr.LDA_Compute2(evObj, event_matrix1, event_matrix2)

        #Finding at which level (coefficient) the behavior of event is different between classes
        evObj.Level_difference(event_matrix1,event_matrix2)

        ##dr.Regression_Logistic([event[0]],event_matrix1,event_matrix2)

        # if( (np.sum(np.absolute(event_matrix1)) == 0 ) or (np.sum(np.absolute(event_matrix2)) == 0 ) ):
        #     continue
        # ttest_pvalues[event[0]] = dr.Ttset_Variants(np.array(event_matrix1), np.array(event_matrix2))
        # viz.ThreeDPlot(event[0])

    Print_Screen_Off('on')
    print("Computing duration time is started")
    Print_Screen_Off('off')

    #Computing Duration time differences
    dv.Deviation_Timestamp(path1,path2)

    Print_Screen_Off('on')
    print("Plotting fingerprints just started!")
    Print_Screen_Off('off')
    vs.Deviant_Plot_ControlFlow_Time('log','log')

    print "Sorted vales:\n"
    pprint.pprint(sorted_events)
    print "The pvalues are:\n",
    pprint.pprint(ttest_pvalues)



