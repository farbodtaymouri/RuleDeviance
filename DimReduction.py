import numpy as np
from sklearn.decomposition import PCA
from numpy.linalg import inv
from scipy.spatial import distance
from scipy import stats
import Preparation as pr
import pprint
import Reading_Wriitng as rw
import seaborn as sns
import matplotlib.pyplot as plt
import operator
import Variables as v
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from nltk.util import ngrams
import collections
from scipy.stats import chisquare, chi2_contingency, t
from scipy.special import expit
from sklearn.metrics.cluster import contingency_matrix
import pandas as pd
from sklearn.svm import SVC
from sklearn.metrics import f1_score
import gc
import copy
import tqdm
import sys
import os
from sklearn import preprocessing



def Dimensionality_Reduction(m):
    '''

    :param m: Design matrix, where each row shows an observed trace, and columns show Wavelet coefficients
    :return:
    '''
    #If the dimension is greater than 3 we use PCA
    if(m.shape[1]>5):
        # Computing PCA with 100% of total variation in the dataset
        # Perturbing the dataset to avoid of getting singular matrix
        np.random.seed(100)
        m = m + np.random.uniform(size=m.shape) / 1000
        #Computing PCA with .9 of total variation in the dataset
        pc = PCA(0.9)
        m_reduced = pc.fit_transform(m)
        print "The dimension before reduction:", m.shape
        print "The dimension of the reduced matrix:", m_reduced.shape


        #Computing covariance matrix of the reduced dataset (numpy.cov(), considers the features are presented in rows)
        cov_reduced = np.cov(m_reduced.T)
        #Computing the inverse of covariance matrix
        try:
            inv_cov_reduced = np.round( inv(cov_reduced),3)
        except np.linalg.LinAlgError :
            print "The result of PCA has only one column! No covariance matrix"
            inv_cov_reduced = None

        #Compute mean across the columns
        mean_col = np.round( np.mean(m_reduced, axis= 0), 3 )

        v.m_reduced = m_reduced

        #When we want to transform a new vector x into this dimension we use pc.transform(x.reshape([1,-1])
    else:

        # Computing PCA with 100% of total variation in the dataset
        # Perturbing the dataset to avoid of getting singular matrix
        np.random.seed(100)
        m = m + np.random.uniform(size=m.shape) / 10000
        pc = PCA()
        m_reduced = pc.fit_transform(m)
        print "The dimension before reduction:", m.shape
        print "The dimension of the reduced matrix:", m_reduced.shape

        # Computing covariance matrix of the reduced dataset (numpy.cov(), considers the features are presented in rows)
        cov_reduced = np.cov(m_reduced.T)
        # # Computing the inverse of covariance matrix
        # inv_cov_reduced = np.round(inv(cov_reduced), 3)

        # Computing the inverse of covariance matrix
        try:
            inv_cov_reduced = np.round(inv(cov_reduced), 3)
        except np.linalg.LinAlgError:
            print "The result of PCA has only one column! No covariance matrix"
            inv_cov_reduced = None



        v.m_reduced = m_reduced

        # Compute mean across the columns
        mean_col = np.round(np.mean(m_reduced, axis=0), 3)



    return inv_cov_reduced, mean_col, pc

#####################################################################################

def Compute_distance(inv_cov_reduced, mean_col, pc, matrix):
    '''

    :param inv_cov_reduced: Inverse covariance matrix, a numpy 2D array
    :param mean_col: Mean of the dataset, a 1D numpy array
    :param pc:  Object of PCA
    :param log:  The event log foe which we want to compute the distance (it must be transformed into multidimensional space)
    :return:
    '''
    distance_values = []

    #Transform the log into the space we want to compare
    matrix_transform = pc.transform(matrix)
    v.m2_reduced = matrix_transform


    #Computing Mahalanobis distance for every single trace
    #If the result of PCA transformation results in one column then no mahalonobis distance is computed
    if (inv_cov_reduced is None):
        #Flatenning the the since it looks like log_transform =[ [1.23],[2.3],[1.76],..]
        distance_values = [dist for element in matrix_transform  for dist in element]
    else:
        for trace in matrix_transform:
            distance_values.append(distance.mahalanobis(trace, mean_col, inv_cov_reduced))

    return distance_values

#####################################################################################

#Apply t-test between two event logs, represented as two design matrices.
def  Ttset_Variants(matrix1, matrix2):
    if ( (matrix1.shape[0] <=5) or (matrix2.shape[0] <=5)):
        raise Exception("The number of observations is not enough for t-test!", matrix1.shape[0], matrix2.shape[0] )
    #Compute the distances from vectors of matrix 2 to matrix1
    inv_cov_reduced, mean_col, pc = Dimensionality_Reduction(matrix1)
    di1 = Compute_distance(inv_cov_reduced, mean_col, pc, matrix2)

    #Computing the distances from vectors of matrix 1 to itself
    di2 = Compute_distance(inv_cov_reduced, mean_col, pc, matrix1)

    # Compute the distances from vectors of matrix 1 to matrix2
    # inv_cov_reduced, mean_col, pc = Dimensionality_Reduction(matrix2)
    # di2 = Compute_distance(inv_cov_reduced, mean_col, pc, matrix1)

    print"distance 1:\n", di1
    print"distance 2:\n", di2
    print "Covariance matrix 1:\n", inv_cov_reduced
    #Apply t-test
    t = stats.ttest_ind(di1, di2, equal_var=False)
    print "The P-value of t-test:", np.round( t.pvalue, 9)
    return np.round( t.pvalue, 9)

######################################################################################
#Computing which features better can discriminate the classes
def Feature_Selection(matrix1, matrix2, obj):

    mean_gloabal = np.mean(np.concatenate((matrix1,matrix2)))

    #Augmenting matrices with the class labels
    # matrix1 = [matrix1, np.repeat(0, matrix1.shape[0]).T]
    # matrix2 = [matrix2, np.repeat(1, matrix2.shape[0]).T]


    # #Creating the column names
    # # For example obj.dic_WaveletCfVector =
    # # {'c': [0.25, 0.25, 0.5, 0.0],
    # #  'd': [0.0, 0.0, 0.0, 0.0],
    # #  'e': [0.0, 0.0, 0.0, 0.0],
    # #  'f': [0.25, 0.25, -0.5, 0.0]}
    # colnames =[]
    # for event in obj.dic_WaveletCfVector.keys():
    #     for i in range(len(obj.dic_WaveletCfVector[event])):
    #         colnames.append(event+"*"+str(i))
    #
    # print("Matrix1:",matrix1)
    # print("Matrix2:", matrix2)
    # print ("The column names is:", colnames)
    #
    #
    # #A vector that shows mean values for each column
    # mean_matrix1 = np.mean(matrix1,axis=0)
    # print("mean matrix1:", mean_matrix1)
    # mean_matrix2 = np.mean(matrix2, axis=0)
    # var_matrix1 = np.var(matrix1, axis=0) +np.random.uniform(low=0,high=0.002, size=matrix1.shape[1])
    # var_matrix2 = np.var(matrix2, axis=0) +np.random.uniform(low=0,high=0.002, size=matrix2.shape[1])
    portion1 = float(matrix1.shape[0]) /(matrix1.shape[0]+matrix2.shape[0])
    portion2 = float(matrix2.shape[0]) / (matrix1.shape[0] + matrix2.shape[0])
    #
    #
    fisher_score = {}
    # for i in range(len(colnames)):
    #     print "The variable under consideration is:", colnames[i]
    #     numerator = 0
    #     denominator = 0
    #
    #     numerator+= np.power( (mean_matrix1[i] - mean_gloabal), 2) *portion1; print"The numerator for class 1:", numerator
    #     numerator += np.power((mean_matrix2[i] - mean_gloabal), 2) * portion2; print"The numerator for class 2:", numerator
    #
    #
    #     denominator += var_matrix1[i]*portion1; print"The dinomenator for class 1:", denominator
    #     denominator += var_matrix2[i]*portion2; print"The dinomenator for class 2:", denominator
    #
    #     fisher_score[colnames[i]] = float(numerator)/denominator

    k=0
    dim1 = matrix1.shape
    dim2 =matrix2.shape
    matrix1_tall = []
    for event in obj.dic_WaveletCfVector.keys():
        numerator=0
        denominator=0
        step = len(obj.dic_WaveletCfVector[event])
        print"event:",event
        print "step:",step
        print "Dim:",dim1
        print "matrix:", matrix1[:, range(step*k, (k+1)*step)]
        temp1= matrix1[:, range(step*k, (k+1)*step)].reshape((dim1[0]*step, 1),order="F")
        temp2 = matrix2[:, range(step * k, (k + 1) * step)].reshape((dim2[0] * step, 1), order="F")
        print "temp1:",temp1
        print "temp2:",temp2
        print "temp1 Con temp2:", np.concatenate((temp1,temp2))
        mean_event = np.mean(np.concatenate((temp1, temp2)))

        numerator = np.power((np.mean(temp1) - mean_event), 2)*portion1  +  np.power((np.mean(temp2) - mean_event), 2)* portion2
        denominator = np.var(temp1) * portion1  + np.var(temp2)*portion2

        fisher_score[event] = float(numerator) / denominator
        # matrix1_tall.append(temp1)
        # matrix1_tall = np.array(matrix1_tall).flatten()
        k+=1




    print ("The tall matrix", fisher_score)
    print("sorted dic is:"), sorted(fisher_score,key=fisher_score.__getitem__)

#####################################################################################

#Implementing the feature:value representation for sparse dataset
def FeatureValue_Representation(log, log2):
    '''

    :param log: A list of list like [['a', 'b']] * 6 + [['b', 'a', 'c','d']] * 3
    :return: A nested dictionary , where the main key is the case ID, other keys are events with non zero list
    {9: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
 10: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
 11: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
 12: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
 13: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
 14: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
 15: {'a': [0.25, 0.25, -0.5, 0.0],
      'b': [0.25, 0.25, 0.5, 0.0],
      'c': [0.25, -0.25, 0.0, 0.5]},
 16: {'a': [0.25, 0.25, -0.5, 0.0],
      'b': [0.25, 0.25, 0.5, 0.0],
      'c': [0.25, -0.25, 0.0, 0.5]},
 17: {'a': [0.25, 0.25, -0.5, 0.0],
      'b': [0.25, 0.25, 0.5, 0.0],
      'c': [0.25, -0.25, 0.0, 0.5]}}
    '''


    # log = [['a', 'b']] * 6 + [['b', 'a', 'c','d']] * 3
    # log2 = [['a', 'c']] * 6 + [['b', 'a', 'c']] * 3
    # #
    # log, id = rw.Reading_log_all('C:/Users/taymouri/Desktop/model temp/trash2/Road_teraffic_fineamount_Above and equal50.xes', d)
    # log2, id = rw.Reading_log_all('C:/Users/taymouri/Desktop/model temp/trash2/Road_teraffic_fineamount_Below50.xes', d)
    # print "Reading log is done!!"


    # log = Log_Stat(log, cut='yes')
    # log2 = Log_Stat(log2, cut='yes')



    # Normalizing the log (equating the length)
    obj1 = pr.Seq2Vec()
    log_tot = obj1.Log_Normalization(log + log2)


    log_featureValue_dic={}
    log2_featureValue_dic={}

    # for trace in log:
    for i in tqdm.tqdm(range(len(log_tot))):
        if (i <= (len(log) - 1)):
            trace = log_tot[i]
            obj = pr.Seq2Vec()  # Creating an object for the trace
            obj.unique_Events_log = obj1.unique_Events_log  # Giving the unique name of event log
            obj.Start(trace)  # Starting Sequence to vector
            # Flattening the dictionary of coefficients to make a matrix row
            print "the cf:", obj.dic_WaveletCfVector.values()
            # obj.dic_WaveletCfVector =
            # {'a': [0.25, 0.25, -0.5, 0.0],
            #  'b': [0.25, 0.25, 0.5, 0.0],
            #  'c': [0.25, -0.25, 0.0, 0.5],
            #  'd': [0.0, 0.0, 0.0, 0.0]}

            #We only select those key, foe which the list is not totally zero
            nonzero_keys = [key  for key in obj.dic_WaveletCfVector  if np.sum(np.absolute(obj.dic_WaveletCfVector[key])) !=0]
            log1_id =v.id1.keys()
            log_featureValue_dic[log1_id[i]] = dict((k, obj.dic_WaveletCfVector[k]) for k in nonzero_keys)
            #log_featureValue_dic[i] = dict((k, obj.dic_WaveletCfVector[k]) for k in nonzero_keys)


            #del obj

        else:
            trace = log_tot[i]
            obj = pr.Seq2Vec()  # Creating an object for the trace
            obj.unique_Events_log = obj1.unique_Events_log  # Giving the unique name of event log
            obj.Start(trace)  # Starting Sequence to vector
            # Flattening the dictionary of coefficients to make a matrix row
            print "the cf:", obj.dic_WaveletCfVector.values()

            # We only select those key, foe which the list is not totally zero
            nonzero_keys = [key for key in obj.dic_WaveletCfVector if np.sum(np.absolute(obj.dic_WaveletCfVector[key])) != 0]
            log2_id = v.id2.keys()
            log2_featureValue_dic[log2_id[i-len(log)]] = dict((k, obj.dic_WaveletCfVector[k]) for k in nonzero_keys)
            #log2_featureValue_dic[i] = dict((k, obj.dic_WaveletCfVector[k]) for k in nonzero_keys)

    # print"The computed matrix is:\n", pprint.pprint(obj.matrix_transform)
    # print"The computed inverse matrix is:\n", pprint.pprint(obj.matrix_transform_inv)
    print "The Matrix is:\n", pprint.pprint(obj.matrix_transform)
    print "The inverse Matrix is:\n", pprint.pprint(obj.matrix_transform_inv)
    print"the log dictionary one is:\n", pprint.pprint(log_featureValue_dic)
    print "the log dictionary two is:\n", pprint.pprint(log2_featureValue_dic)
    return log_featureValue_dic, log2_featureValue_dic

#################################################################################
#This function only provides statistics for the log (such as maximum lenght and distribution)
def Log_Stat():
    '''

    :param log: log = [['a', 'b']] * 6 + [['b', 'a', 'c','d']]
    :param cut: A string 'yes' or 'no'
    :return: The truncated log
    '''


    lengths=[]
    for i in range(len(v.log1)):
        lengths += len(v.id1[i].split(",")) * [len(v.log1[i])]

    #lengths = [len(e) for e in v.log1]
    print "The number of traces 1 is:", len(lengths)
    print "The maximum length of the trace is:", np.max(lengths)
    print "The average length of the trace is:", np.mean(lengths)
    print "The standard deviation of the trace is:", np.std(lengths)
    print "The minimum length of the trace is:", np.min(lengths)
    #sns.distplot(lengths)
    #plt.show()

    log1=[]
    id1={}
    k=0
    for i in range(len(v.log1)):
        if len(v.log1[i]) <= np.mean(lengths) + 2 * np.std(lengths):
            log1.append(v.log1[i])
            id1[k] = v.id1[i]
            k+=1



    #-------------------------------------------------
    lengths = []
    for i in range(len(v.log2)):
        lengths += len(v.id2[i].split(",")) * [len(v.log2[i])]

    #lengths = [len(e) for e in v.log2]
    print "The number of traces 1 is:", len(lengths)
    print "The maximum length of the trace is:", np.max(lengths)
    print "The average length of the trace is:", np.mean(lengths)
    print "The standard deviation of the trace is:", np.std(lengths)
    print "The minimum length of the trace is:", np.min(lengths)
    #sns.distplot(lengths)
    #plt.show()

    log2=[]
    id2={}
    k=0
    for i in range(len(v.log2)):
        if len(v.log2[i]) <= np.mean(lengths) + 2 * np.std(lengths):
            log2.append(v.log2[i])
            id2[k] = v.id2[i]
            k+=1


    v.log1 = copy.deepcopy(log1)
    v.id1 = copy.deepcopy(id1)

    v.log2 =copy.deepcopy( log2)
    v.id2 = copy.deepcopy(id2)

    # if(cut == 'yes'):
    #     log = [e for e in log if len(e) <= np.mean(lengths)+2*np.std(lengths) ]
    #
    #
    # return log
    gc.collect()
    return log1,log2, id1, id2

################################################################################
#This function ranks the events in terms of discrimination (Fisher score)
def Feature_Ranking(log_featureValue_dic, log2_featureValue_dic):

    '''
    :param log_featureValue_dic:
            A nested dictionary , where the main key is the case ID, other keys are events with non zero list
            {9: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
         10: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
         11: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
         12: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
         13: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
         14: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
         15: {'a': [0.25, 0.25, -0.5, 0.0],
              'b': [0.25, 0.25, 0.5, 0.0],
              'c': [0.25, -0.25, 0.0, 0.5]},
         16: {'a': [0.25, 0.25, -0.5, 0.0],
              'b': [0.25, 0.25, 0.5, 0.0],
              'c': [0.25, -0.25, 0.0, 0.5]},
         17: {'a': [0.25, 0.25, -0.5, 0.0],
              'b': [0.25, 0.25, 0.5, 0.0],
              'c': [0.25, -0.25, 0.0, 0.5]}}
    :param log2_featureValue_dic: Similar to the above
    :return: An ordered list like [('c', 0.014702191780821921), ('a', 0.01254339888714726), ('b', 0.004671108599976933), ('d', 0.0)]
    '''

    #Finding unique events
    events_unique = set([e for trace in log2_featureValue_dic for e in log2_featureValue_dic[trace]])
    events_unique = events_unique.union(set([e for trace in log_featureValue_dic for e in log_featureValue_dic[trace]]))
    print "Unique events are:",events_unique


    #Computing scores for each event
    portion1 = float(len(log_featureValue_dic)) / (float(len(log_featureValue_dic))+float(len(log2_featureValue_dic)) )
    portion2 = float(len(log2_featureValue_dic)) / (float(len(log_featureValue_dic)) + float(len(log2_featureValue_dic)))


    fscores={}
    for e in events_unique:

        el1,el2=[],[]
        for key in log_featureValue_dic:
            if(e in log_featureValue_dic[key]):
                el1 += log_featureValue_dic[key][e]
            else:
                pass

        for key in log2_featureValue_dic:
            if(e in log2_featureValue_dic[key]):
                el2 += log2_featureValue_dic[key][e]
            else:
                pass


        #This is for the case that an event only happened in one event log
        if(len(el1)==0):
            el1 = [0]*len(el2)
            tot_mean = np.mean(el2+ el1)
            mean1 = np.mean(el1)
            mean2 = np.mean(el2)
            fscores[e] = ( portion1*np.power((mean1 - tot_mean),2) + portion2*np.power((mean2-tot_mean) ,2) /  portion1*np.var(el1) +portion2*np.var(el2)   )
            fscores[e] = np.power(mean2,2) /  np.var(el2)



        elif (len(el2)==0):
            el2 = [0] * len(el1)
            mean1 = np.mean(el1)
            mean2 = np.mean(el2)
            tot_mean = np.mean(el2+ el1)
            fscores[e] = (portion1 * np.power((mean1 - tot_mean), 2) + portion2 * np.power((mean2 - tot_mean),
                                                                                     2) / portion1 * np.var(
                el1) + portion2 * np.var(el2))
            fscores[e] = np.power(mean1, 2) / np.var(el1)


        else:
            mean1 = np.mean(el1)
            mean2 = np.mean(el2)
            tot_mean = np.mean(el2+ el1)
            fscores[e] = (portion1 * np.power((mean1 - tot_mean), 2) + portion2 * np.power((mean2 - tot_mean),
                                                                                     2) / portion1 * np.var(el1) + portion2 * np.var(el2))


    sorted_events = sorted(fscores.items(), key=operator.itemgetter(1),reverse=True)
    print("sorted dic is:"),  sorted_events
    return  sorted_events
################################################################################
#This function given a set of events create two matrices (related to event log1 and log2)
def Event_Matrxi_Creation(event_list, log_featureValue_dic, log2_featureValue_dic ):
    '''

    :param event_list: A list of events like ['a','b']
    :param log_featureValue_dic: log_featureValue_dic:
            A nested dictionary , where the main key is the case ID, other keys are events with non zero list
            {9: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
         10: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
         11: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
         12: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
         13: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
         14: {'a': [0.25, 0.25, 0.5, 0.0], 'c': [0.25, 0.25, -0.5, 0.0]},
         15: {'a': [0.25, 0.25, -0.5, 0.0],
              'b': [0.25, 0.25, 0.5, 0.0],
              'c': [0.25, -0.25, 0.0, 0.5]},
         16: {'a': [0.25, 0.25, -0.5, 0.0],
              'b': [0.25, 0.25, 0.5, 0.0],
              'c': [0.25, -0.25, 0.0, 0.5]},
         17: {'a': [0.25, 0.25, -0.5, 0.0],
              'b': [0.25, 0.25, 0.5, 0.0],
              'c': [0.25, -0.25, 0.0, 0.5]}}
    :param log2_featureValue_dic: Similar to the above
    :return: two matrix. like :
                            [[0.2, 0.133, 0.167, 0.5, 0.0],
                            [0.2, 0.133, 0.167, 0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0]]
    '''

    #dimension of each element vector
    t_id = log_featureValue_dic.keys()[0]
    el =log_featureValue_dic[t_id].keys()[0]
    event_dim =  len(log_featureValue_dic[t_id][el])

    #Iterate over the first event log
    event_matrix1=[]
    for trace in log_featureValue_dic:
        #Fisrt finding how many times it happened
        freq= len(v.id1[trace].split(","))
        print"Tjhe freq:",freq
        temp=[]
        for e in event_list:
            if(e in log_featureValue_dic[trace]):
                temp+=log_featureValue_dic[trace][e]
            # else:
            #     temp+=[0]*event_dim
        [event_matrix1.append(t) for t in [temp] * freq]
        #event_matrix1.append(temp)

    #Iterate over the second event log
    event_matrix2 = []
    for trace in log2_featureValue_dic:
        # Fisrt finding how many times it happened
        freq = len(v.id2[trace].split(","))
        print"Tjhe freq22:", freq
        temp = []
        for e in event_list:
            if (e in log2_featureValue_dic[trace]):
                temp += log2_featureValue_dic[trace][e]
            # else:
            #     temp += [0] * event_dim
        [event_matrix2.append(t) for t in [temp] * freq]
        #event_matrix2.append(temp)


    print"The matrix according to elements:", event_list
    pprint.pprint(event_matrix1)
    print"---------------------"
    pprint.pprint(event_matrix2)

    return event_matrix1, event_matrix2

###############################################################################
#This function computes LDA vector
def LDA_Compute(events,matrix1, matrix2):
    '''

    :param matrix1: like :[[0.2, 0.133, 0.167, 0.5, 0.0],
                            [0.2, 0.133, 0.167, 0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0]]
    :param matrix2: Similar oto the above matrix
    :return:
    '''

    # #-------------------------------------------------------
    # ##### AS aTEST, Just considering rows that are non-zero
    # zero_row_index = np.where(np.all(np.array(matrix1) == 0, axis=1))[0]
    # matrix1 = [matrix1[i] for i in range(len(matrix1)) if i not in zero_row_index]
    # #Regularzring matrix1
    # print "Matrix1 after removing zero rows:\n", matrix1
    #
    # #for matrix 2
    # zero_row_index = np.where(np.all(np.array(matrix2) == 0, axis=1))[0]
    # matrix2 = [matrix2[i] for i in range(len(matrix2)) if i not in zero_row_index]
    # print "Matrix2 after removing zero rows:\n", matrix2
    #
    # #When the event happens only in one class
    # if(matrix2 == []):
    #     print("The element:", str(events), ", happened in only class 1");return
    # elif (matrix1 ==[]):
    #     print("The element:", str(events), ", happened in only class 2");return
    #
    # #----------------------------------------------------

    X = np.concatenate((np.array(matrix1), np.array(matrix2)))
    print "Before smoothing:",X
    #X = preprocessing.normalize(X)
    X=np.log(X+10)
    #Smoothing via moving average
    #See https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    # for i in range(X.shape[1]):
    #     N=X.shape[0]
    #     X[:,i] = np.convolve( X[:,i], np.ones((N,))/N,'same' )
    Y = [0]*len(matrix1) + [1]*len(matrix2)
    prob_class1 = len(matrix1) / float(X.shape[0])
    prob_class2 = len(matrix2) / float(X.shape[0])





    clf = LinearDiscriminantAnalysis(tol=0.01, priors=[prob_class1,prob_class2])

    X_lda = clf.fit_transform(X, Y)
    print"The initial X is:\n",X
    print "The LDA component is:\n",X_lda
    print "The coeff is:", clf.coef_
    #print "The probabilities are:", clf.predict_proba(X)
    print "The predictions using expit:", expit(X_lda)
    print "The predicted class labels directly:", clf.predict(X)
    u=np.zeros(len(X_lda))
    #u[np.where(np.array(expit(X_lda))>0.5)[0]] =1
    #u[np.where(np.array(expit(clf.decision_function(X))) > 0.5)[0]] = 1
    u[np.where(np.array(expit(X_lda)) > 0.5)[0]] = 1
    u= np.array(clf.predict(X))
    print "The predicted labels u:",list(u)
    print "the correct labels Y:",Y



    # #Computing KL Divergence
    # prob_class1 = len(matrix1) / float(X.shape[0])
    # prob_class2 = len(matrix2) / float(X.shape[0])
    # print "The priors are:", prob_class1,prob_class2
    # probabilities= np.array(clf.predict_proba(X))
    # probabilities[0:len(matrix1),0] = probabilities[0:len(matrix1),0]/prob_class1
    # probabilities[0:len(matrix1), 1] = probabilities[0:len(matrix1), 1] / prob_class2
    # probabilities[len(matrix1):X.shape[0], 0] = probabilities[len(matrix1):X.shape[0], 0] / prob_class1
    # probabilities[len(matrix1):X.shape[0], 1] = probabilities[len(matrix1):X.shape[0], 1] / prob_class2
    #
    # print "The probabilities for KL are:", probabilities
    #Measuring goodness of fit with chi square
    ##target_predicted = collections.Counter(clf.predict(X)); print "Target predicted:",target_predicted
    target_predicted = collections.Counter(u) ; print "Target predicted:",target_predicted
    target_true = collections.Counter(Y) ; print "Target true values:",target_true

    print "Contingency matrix is:", contingency_matrix(Y, u, eps=None, sparse=False)
    # CM_Chi2_Test(contingency_matrix(Y, u, eps=None, sparse=False))
    #Sometime the contingency matrix provided by python comes up with 2*1 matrix instead of 2*2, this is why I developed it manually
    z = [[0,0],[0,0]]
    for i in range(2):
        for j in range(2):
            for x in zip(Y,u):
                if((i,j) == x):
                    z[i][j]+=1
    print "The handmade contingency matrix is:", z
    CM_Chi2_Test(z)
    Model_Test2(Y,u)

    #print "Cont table result:", chi2_contingency(contingency_matrix(Y, u, eps=None, sparse=False))

    keys = set(target_predicted.keys() + target_true.keys())
    obs=[]
    trues=[]
    for  cls in keys:
        if cls in target_predicted.keys():
            obs.append(target_predicted[cls])
        else:
            obs.append(0)

        if cls in target_true.keys():
            trues.append(target_true[cls])
        else:
            trues.append(0)

    #print "the observed and true are:",obs,trues
    print "the obs:",obs
    print "The rue:",trues
    #print "The result of chi test:",chisquare(obs, f_exp=trues)
    #CM_Chi2_Test([trues,obs])
    #print "NEW:", chisquare(obs,trues)
    #print"Neww2:",CM_Chi2_Test([trues,obs])





    #This is for jittering
    y_axis = np.random.uniform(size=X.shape[0])/100
    x_axis = np.random.uniform(size=X.shape[0])/10000
    col = ['c1']*len(X_lda[0:len(matrix1), 0]) + ['c2']*len(X_lda[len(matrix1):,0])

    df1 = pd.DataFrame({'x': X_lda[:,0], 'y':y_axis, 'col':col})
    sns.set_style("darkgrid")
    sns.lmplot(x='x', y='y', data=df1, fit_reg=False, hue='col', legend=False, palette={'c1':'red', 'c2':'blue'})

    # plt.plot(X_lda[0:len(matrix1), 0]+x_axis[0:len(matrix1)], y_axis[0:len(matrix1)] , 'gx')
    # plt.plot(X_lda[len(matrix1):,0]+x_axis[len(matrix1):], y_axis[len(matrix1):],'ro')
    plt.title("LDA projection for:"+str(events))
    plt.show()

#############################################################################


###############################################################################
#This function computes LDA vector
def LDA_Compute2(evObj,matrix1, matrix2):

    # sys.stdout = open(os.devnull, "w")
    # f = open(os.devnull, 'w')
    # sys.stdout = f
    '''

    :param matrix1: like :[[0.2, 0.133, 0.167, 0.5, 0.0],
                            [0.2, 0.133, 0.167, 0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0]]
    :param matrix2: Similar oto the above matrix
    :return:
    '''

    # #-------------------------------------------------------
    # ##### AS aTEST, Just considering rows that are non-zero
    # zero_row_index = np.where(np.all(np.array(matrix1) == 0, axis=1))[0]
    # matrix1 = [matrix1[i] for i in range(len(matrix1)) if i not in zero_row_index]
    # #Regularzring matrix1
    # print "Matrix1 after removing zero rows:\n", matrix1
    #
    # #for matrix 2
    # zero_row_index = np.where(np.all(np.array(matrix2) == 0, axis=1))[0]
    # matrix2 = [matrix2[i] for i in range(len(matrix2)) if i not in zero_row_index]
    # print "Matrix2 after removing zero rows:\n", matrix2
    #
    # #When the event happens only in one class
    # if(matrix2 == []):
    #     print("The element:", str(events), ", happened in only class 1");return
    # elif (matrix1 ==[]):
    #     print("The element:", str(events), ", happened in only class 2");return
    #
    # #----------------------------------------------------

    # Some features only happened in one class, therefore no need to run a classifier for them
    # like matrix1=[[],...[]]
    if (np.sum([0 if t == [] else 1 for t in matrix1]) == 0):
        evObj.pValue = 0
        evObj.class_proportion = {"matrix1": 0, "matrix2:": 1}
        return
    if (np.sum([0 if t == [] else 1 for t in matrix2]) == 0):
        evObj.pValue = 0
        evObj.class_proportion = {"matrix1": 1, "matrix2:": 0}
        return
    else:
        matrix1= [t for t in matrix1 if t != []]
        matrix2 = [t for t in matrix2 if t != []]


    X = np.concatenate((np.array(matrix1), np.array(matrix2)))
    print "Before smoothing:",X
    #X = preprocessing.normalize(X)
    #X=np.log(X+10)
    #Smoothing via moving average
    #See https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    # for i in range(X.shape[1]):
    #     N=X.shape[0]
    #     X[:,i] = np.convolve( X[:,i], np.ones((N,))/N,'same' )
    Y = [0]*len(matrix1) + [1]*len(matrix2)
    prob_class1 = len(matrix1) / float(X.shape[0])
    prob_class2 = len(matrix2) / float(X.shape[0])





    clf = LinearDiscriminantAnalysis(tol=0.01, priors=[prob_class1,prob_class2])

    X_lda = clf.fit_transform(X, Y)
    #print"The initial X is:\n",X
    #print "The LDA component is:\n",X_lda
    print "The coeff is:", clf.coef_
    #print "The probabilities are:", clf.predict_proba(X)
    #print "The predictions using expit:", expit(X_lda)
    #print "The predicted class labels directly:", clf.predict(X)
    u=np.zeros(len(X_lda))
    #u[np.where(np.array(expit(X_lda))>0.5)[0]] =1
    #u[np.where(np.array(expit(clf.decision_function(X))) > 0.5)[0]] = 1
    u[np.where(np.array(expit(X_lda)) > 0.5)[0]] = 1
    u= np.array(clf.predict(X))
    print "The predicted labels u:",list(u)
    print "the correct labels Y:",Y



    # #Computing KL Divergence
    # prob_class1 = len(matrix1) / float(X.shape[0])
    # prob_class2 = len(matrix2) / float(X.shape[0])
    # print "The priors are:", prob_class1,prob_class2
    # probabilities= np.array(clf.predict_proba(X))
    # probabilities[0:len(matrix1),0] = probabilities[0:len(matrix1),0]/prob_class1
    # probabilities[0:len(matrix1), 1] = probabilities[0:len(matrix1), 1] / prob_class2
    # probabilities[len(matrix1):X.shape[0], 0] = probabilities[len(matrix1):X.shape[0], 0] / prob_class1
    # probabilities[len(matrix1):X.shape[0], 1] = probabilities[len(matrix1):X.shape[0], 1] / prob_class2
    #
    # print "The probabilities for KL are:", probabilities
    #Measuring goodness of fit with chi square
    ##target_predicted = collections.Counter(clf.predict(X)); print "Target predicted:",target_predicted
    target_predicted = collections.Counter(u) ; print "Target predicted:",target_predicted
    target_true = collections.Counter(Y) ; print "Target true values:",target_true

    print "Contingency matrix is:", contingency_matrix(Y, u, eps=None, sparse=False)
    # CM_Chi2_Test(contingency_matrix(Y, u, eps=None, sparse=False))
    #Sometime the contingency matrix provided by python comes up with 2*1 matrix instead of 2*2, this is why I developed it manually
    # z = [[0,0],[0,0]]
    # for i in range(2):
    #     for j in range(2):
    #         for x in zip(Y,u):
    #             if((i,j) == x):
    #                 z[i][j]+=1
    # print "The handmade contingency matrix is:", z
    # CM_Chi2_Test(z)
    # Model_Test2(Y,u)

    #print "Cont table result:", chi2_contingency(contingency_matrix(Y, u, eps=None, sparse=False))

    keys = set(target_predicted.keys() + target_true.keys())
    obs=[]
    trues=[]
    for  cls in keys:
        if cls in target_predicted.keys():
            obs.append(target_predicted[cls])
        else:
            obs.append(0)

        if cls in target_true.keys():
            trues.append(target_true[cls])
        else:
            trues.append(0)

    #print "the observed and true are:",obs,trues
    print "the obs:",obs
    print "The rue:",trues
    #print "The result of chi test:",chisquare(obs, f_exp=trues)
    #CM_Chi2_Test([trues,obs])
    #print "NEW:", chisquare(obs,trues)
    #print"Neww2:",CM_Chi2_Test([trues,obs])

    # Storing the corresponding designe matrix
    evObj.design_matrix = np.concatenate((X, np.array([Y]).T), axis=1)
    evObj.class_proportion = np.average(Y)
    Model_Test2(evObj, Y, clf.predict(X))



    # #This is for jittering
    # y_axis = np.random.uniform(size=X.shape[0])/100
    # x_axis = np.random.uniform(size=X.shape[0])/10000
    # col = ['c1']*len(X_lda[0:len(matrix1), 0]) + ['c2']*len(X_lda[len(matrix1):,0])
    #
    # df1 = pd.DataFrame({'x': X_lda[:,0], 'y':y_axis, 'col':col})
    # sns.set_style("darkgrid")
    # sns.lmplot(x='x', y='y', data=df1, fit_reg=False, hue='col', legend=False, palette={'c1':'red', 'c2':'blue'})
    #
    # # plt.plot(X_lda[0:len(matrix1), 0]+x_axis[0:len(matrix1)], y_axis[0:len(matrix1)] , 'gx')
    # # plt.plot(X_lda[len(matrix1):,0]+x_axis[len(matrix1):], y_axis[len(matrix1):],'ro')
    # plt.title("LDA projection for:"+str(events))
    # plt.show()

#############################################################################
#Generating the Ngram of the given log
def Ngram_Compute(log, k=2):
    '''

    :param log: [['a', 'b','e','f','d']] * 4 + [['b', 'a', 'c']] * 3
    :param k: A parameter that shows unigram or bigram
    :return: [[('a', 'b'), ('b', 'e'), ('e', 'f'), ('f', 'd')],
             [('a', 'b'), ('b', 'e'), ('e', 'f'), ('f', 'd')],
             [('a', 'b'), ('b', 'e'), ('e', 'f'), ('f', 'd')],
             [('a', 'b'), ('b', 'e'), ('e', 'f'), ('f', 'd')],
             [('b', 'a'), ('a', 'c')],
             [('b', 'a'), ('a', 'c')],
             [('b', 'a'), ('a', 'c')]]
    '''

    temp=[]
    for trace in log:
        temp.append(list(ngrams(trace,k)))

    return temp

#############################################################################
def CM_Chi2_Test(table):
    '''
    :param table: like: [[2 3]
                    [0 6]]
    :return:
    '''

    print"The fisher Exact test:", stats.fisher_exact(table)

    # table = np.array(table)
    # expected_table = np.zeros(table.shape)
    # expected_table[0,0] = np.sum(table[0,:])/float(np.sum(table))
    # expected_table[1, 1] = np.sum(table[1, :]) / float(np.sum(table))
    #
    # temp=0
    # temp+=np.power((expected_table[0,0] -table[0,0]), 2) /expected_table[0,0]
    # temp += np.power((expected_table[1, 1] - table[1, 1]), 2) / expected_table[1, 1]
    # print "Temp:",temp
    # observed_values = [expected_table[0,0],expected_table[1, 1] ,expected_table[1, 0]  ]
    # expected_values =[table[0, 0], table[1,1], table[1,0]]
    # print"gggggggggggggggggg:", chisquare(observed_values, f_exp=expected_values)

    # print "ffff:",list(expected_table.flatten()), table.flatten()
    # print "The result of chi test With handmade conf matrix:", chisquare(list(expected_table.flatten()), f_exp=list(table.flatten()))

    #Compute the expected values for each cell
    # table = np.array(table)
    # expected_table = np.zeros(table.shape)
    # expected_table[0,:] = table[0,:]*np.sum(table,axis=1)[0]/float(np.sum(table))
    # expected_table[1, :] = table[1, :] * np.sum(table, axis=1)[1] / float(np.sum(table))

    #print "The result of chi test NEWWWWWWW:", chi2_contingency(table.flatten(), f_exp=expected_table.flatten())



############################################################################
#Implementing logisitc regression
def Regression_Logistic(events, matrix1,matrix2):

    '''

    :param matrix1: like :[[0.2, 0.133, 0.167, 0.5, 0.0],
                            [0.2, 0.133, 0.167, 0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0]]
    :param matrix2: Similar oto the above matrix
    :return:
    '''

    X = np.concatenate((np.array(matrix1), np.array(matrix2)))
    # regularzing
    # X = X + np.random.uniform(size=X.shape)/100000
    Y = [0] * len(matrix1) + [1] * len(matrix2)
    prob_class1 = len(matrix1) / float(X.shape[0])
    prob_class2 = len(matrix2) / float(X.shape[0])

    # logreg = LogisticRegression(penalty='l1', class_weight='balanced')
    # logreg.fit(X, Y)
    # print "The score of logisitc regression:", logreg.predict(X)
    # # logit_model = sm.Logit(Y, X)
    # # result = logit_model.fit()
    # # print(result.summary2())

#######################################################################
def SVM_Classifier(evObj, matrix1,matrix2):
    '''

    :param matrix1: like :[[0.2, 0.133, 0.167, 0.5, 0.0],
                            [0.2, 0.133, 0.167, 0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0]]
    :param matrix2: Similar oto the above matrix
    :return:
        '''

    # sys.stdout = open(os.devnull, "w")
    # f = open(os.devnull, 'w')
    # sys.stdout = f

    # Some features only happened in one class, therefore no need to run a classifier for them
    # like matrix1=[[],...[]]
    if (np.sum([0 if t == [] else 1 for t in matrix1]) == 0):
        evObj.pValue = 1
        evObj.class_proportion = {"no_matrix1": 0, "no_matrix2": len([t for t in matrix2 if t != []]), "proportion": None}
        return
    if (np.sum([0 if t == [] else 1 for t in matrix2]) == 0):
        evObj.pValue = 1
        evObj.class_proportion = {"no_matrix1": len([t for t in matrix1 if t != []]), "no_matrix2": 0, "proportion": None}
        return
    else:
        matrix1= [t for t in matrix1 if t != []]
        matrix2 = [t for t in matrix2 if t != []]


    X = np.concatenate((np.array(matrix1), np.array(matrix2)))
    # regularzing
    # X = X + np.random.uniform(size=X.shape)/100000
    Y = [1] * len(matrix1) + [0] * len(matrix2)
    prob_class1 = len(matrix1) / float(X.shape[0])
    prob_class2 = len(matrix2) / float(X.shape[0])

    clf = SVC(C=100000000000, kernel='rbf',
    random_state=None, shrinking=True, tol=0.1, verbose=False)
    clf.fit(X, Y)


    print"The initial X is:\n", list(X)
    print "The SVM component is:\n", clf.decision_function(X)
    #print "The prob is:", clf.densify


    u = np.zeros(len(Y))
    # u[np.where(np.array(expit(X_lda))>0.5)[0]] =1
    u[np.where(np.array(expit(clf.decision_function(X))) > 0.5)[0]] = 1
    print "The predicted labels u:", clf.predict(X)
    print "the correct labels Y:", Y

    # #Computing KL Divergence
    # prob_class1 = len(matrix1) / float(X.shape[0])
    # prob_class2 = len(matrix2) / float(X.shape[0])
    # print "The priors are:", prob_class1,prob_class2
    # probabilities= np.array(clf.predict_proba(X))
    # probabilities[0:len(matrix1),0] = probabilities[0:len(matrix1),0]/prob_class1
    # probabilities[0:len(matrix1), 1] = probabilities[0:len(matrix1), 1] / prob_class2
    # probabilities[len(matrix1):X.shape[0], 0] = probabilities[len(matrix1):X.shape[0], 0] / prob_class1
    # probabilities[len(matrix1):X.shape[0], 1] = probabilities[len(matrix1):X.shape[0], 1] / prob_class2
    #
    # print "The probabilities for KL are:", probabilities
    # Measuring goodness of fit with chi square
    ##target_predicted = collections.Counter(clf.predict(X)); print "Target predicted:",target_predicted
    target_predicted = collections.Counter(u);
    print "Target predicted:", target_predicted
    target_true = collections.Counter(Y);
    print "Target true values:", target_true

    print "Contingency matrix is:", contingency_matrix(Y, u, eps=None, sparse=False)
    # CM_Chi2_Test(contingency_matrix(Y, u, eps=None, sparse=False))
    # Sometime the contingency matrix provided by python comes up with 2*1 matrix instead of 2*2, this is why I developed it manually
    z = [[0, 0], [0, 0]]
    for i in range(2):
        for j in range(2):
            for x in zip(Y, u):
                if ((i, j) == x):
                    z[i][j] += 1
    print "The handmade contingency matrix is:", z

    #Storing the corresponding designe matrix
    evObj.design_matrix= np.concatenate((X, np.array([Y]).T), axis=1)
    print {"no_matrix1":len(matrix1), "no_matrix2":len(matrix2), "proportion": np.average(Y)}
    evObj.class_proportion = {"no_matrix1": len(matrix1), "no_matrix2": len(matrix2), "proportion": np.average(Y)}
    Model_Test2(evObj,Y, clf.predict(X))

####################################################################################

#https://pingouin-stats.org/generated/pingouin.multivariate_normality.html
def multivariate_normality(X, alpha=.05):
    """Henze-Zirkler multivariate normality test.

    Parameters
    ----------
    X : np.array
        Data matrix of shape (n_samples, n_features).
    alpha : float
        Significance level.

    Returns
    -------
    normal : boolean
        True if X comes from a multivariate normal distribution.
    p : float
        P-value.

    See Also
    --------
    normality : Test the univariate normality of one or more variables.
    homoscedasticity : Test equality of variance.
    sphericity : Mauchly's test for sphericity.

    Notes
    -----
    The Henze-Zirkler test has a good overall power against alternatives
    to normality and is feasable for any dimension and any sample size.

    Adapted to Python from a Matlab code by Antonio Trujillo-Ortiz and
    tested against the R package MVN.

    Rows with missing values are automatically removed using the
    :py:func:`remove_na` function.

    References
    ----------
    .. [1] Henze, N., & Zirkler, B. (1990). A class of invariant consistent
       tests for multivariate normality. Communications in Statistics-Theory
       and Methods, 19(10), 3595-3617.

    .. [2] Trujillo-Ortiz, A., R. Hernandez-Walls, K. Barba-Rojo and L.
       Cupul-Magana. (2007). HZmvntest: Henze-Zirkler's Multivariate
       Normality Test. A MATLAB file.

    Examples
    --------
    # >>> import pingouin as pg
    # >>> data = pg.read_dataset('multivariate')
    # >>> X = data[['Fever', 'Pressure', 'Aches']]
    # >>> normal, p = pg.multivariate_normality(X, alpha=.05)
    # >>> print(normal, round(p, 3))
    # True 0.717
    """
    from scipy.stats import lognorm

    # Check input and remove missing values
    X = np.asarray(X)
    assert X.ndim == 2, 'X must be of shape (n_samples, n_features).'
    X = X[~np.isnan(X).any(axis=1)]
    n, p = X.shape
    assert n >= 3, 'X must have at least 3 rows.'
    assert p >= 2, 'X must have at least two columns.'

    # Covariance matrix
    S = np.cov(X, rowvar=False, bias=True)
    S_inv = np.linalg.pinv(S)
    difT = X - X.mean(0)

    # Squared-Mahalanobis distances
    Dj = np.diag(np.linalg.multi_dot([difT, S_inv, difT.T]))
    Y = np.linalg.multi_dot([X, S_inv, X.T])
    Djk = -2 * Y.T + np.repeat(np.diag(Y.T), n).reshape(n, -1) + \
        np.tile(np.diag(Y.T), (n, 1))

    # Smoothing parameter
    b = 1 / (np.sqrt(2)) * ((2 * p + 1) / 4)**(1 / (p + 4)) * \
        (n**(1 / (p + 4)))

    hz = n * 4
    if np.linalg.matrix_rank(S) == p:
        hz = n * (1 / (n**2) * np.sum(np.sum(np.exp(-(b**2) / 2 * Djk))) - 2
                  * ((1 + (b**2))**(-p / 2)) * (1 / n)
                  * (np.sum(np.exp(-((b**2) / (2 * (1 + (b**2)))) * Dj)))
                  + ((1 + (2 * (b**2)))**(-p / 2)))

    wb = (1 + b**2) * (1 + 3 * b**2)
    a = 1 + 2 * b**2
    # Mean and variance
    mu = 1 - a**(-p / 2) * (1 + p * b**2 / a + (p * (p + 2)
                                                * (b**4)) / (2 * a**2))
    si2 = 2 * (1 + 4 * b**2)**(-p / 2) + 2 * a**(-p) * \
        (1 + (2 * p * b**4) / a**2 + (3 * p * (p + 2) * b**8) / (4 * a**4)) \
        - 4 * wb**(-p / 2) * (1 + (3 * p * b**4) / (2 * wb)
                              + (p * (p + 2) * b**8) / (2 * wb**2))

    # Lognormal mean and variance
    pmu = np.log(np.sqrt(mu**4 / (si2 + mu**2)))
    psi = np.sqrt(np.log((si2 + mu**2) / mu**2))

    # P-value
    pval = lognorm.sf(hz, psi, scale=np.exp(pmu))
    normal = True if pval > alpha else False
    return normal, pval

####################################################################
#Implementing t-test for the performance of the classifier
#https://www3.nd.edu/~rjohns15/cse40647.sp14/www/content/lectures/28%20-%20Classifier%20Comparisons.pdf
def Model_Test(events,conf_matrix):
    '''
        :param table: like: [[2 3]
                        [0 6]]
        :return:
    '''

    #Computing the error (precision) of the current classifier
    conf_matrix = np.array(conf_matrix)
    errA = 1 - float(np.sum(np.diag(conf_matrix))) / np.sum(conf_matrix)
    # try:
    #     errA = 1 - float(conf_matrix[1,1]) / np.sum(conf_matrix[:,1])
    # except ZeroDivisionError:
    #     errA = 1
    #     print ("The precision is zero and the classifier can not discriminate classes"); return

    #We assume that the features do not provide any information for the classification, indeed no differences
    errH=.9

    #computing variance of differences
    vard = errA*(1-errA)/ np.sum(conf_matrix) + errH*(1-errH)/ np.sum(conf_matrix)

    #Computing the t-value at specified significant level.
    alpha= 0.95
    df = 1*np.sum(conf_matrix) -1
    c_value = stats.t.ppf((1.0 - alpha)/2, df)
    #c_value = stats.norm.ppf((1.0 - alpha)/2)


    #computing the interval
    print"The classifier eror:",errA
    print "The variance is:", vard
    print "The c-value is:", c_value
    print"The terms:",c_value*np.sqrt(vard)
    upper_value = errH-errA + c_value*np.sqrt(vard)
    lower_value = errH-errA - c_value*np.sqrt(vard)

    print "The interval is:", [lower_value,upper_value]
    if lower_value <= 0 <= upper_value:
        print("No statistical significant")

    v.stat_test_result[str(events)] = [upper_value, lower_value]
    print "The resulrs:\n",
    pprint.pprint(v.stat_test_result)
###################################################################
def Model_Test2(evObj, correct_label, pred_label):


    f1C1 = f1_score(correct_label,pred_label)
    print "f1c1:",f1C1

    correct_label = 1-np.array(correct_label)
    pred_label = 1-np.array(pred_label)

    f1C2= f1_score(correct_label,pred_label)
    print "f1c2:",f1C2

    #Weighted F1 score for the observed data
    prop1 = float(evObj.class_proportion['no_matrix1'])/(evObj.class_proportion['no_matrix1'] +evObj.class_proportion['no_matrix2'])
    prop2 = float(evObj.class_proportion['no_matrix2']) / (
                evObj.class_proportion['no_matrix1'] + evObj.class_proportion['no_matrix2'])
    f1_avg = np.mean([prop1*f1C1, prop2*f1C2])
    #f1_avg=(f1C1+f1C2)/2
    errA = f1_avg


    #Considering the F1 score for the hypothetical (null) classifier
    #errH=0.5
    prop = evObj.class_proportion['proportion']
    #Weithed F1 score
    errH = np.mean([(1-prop)*2*prop/(2*prop+1), prop*2/(2+prop)])
    # computing variance of differences
    vard = errA * (1 - errA) / len(correct_label) + errH * (1 - errH) / len(correct_label)

    # Computing the t-value at specified significant level.
    alpha = 0.95
    df = 1 * len(correct_label) - 1
    c_value = stats.t.ppf((1.0 - alpha) / 2, df)
    # c_value = stats.norm.ppf((1.0 - alpha)/2)

    # computing the interval
    print "The level of H0:", errH
    print"The classifier eror (f1_avg):", errA
    print "The variance is:", vard
    print "The c-value is:", c_value
    print"The terms:", c_value * np.sqrt(vard)
    upper_value = errH - errA + c_value * np.sqrt(vard)
    lower_value = errH - errA - c_value * np.sqrt(vard)


    print "The interval is:", [lower_value, upper_value]
    #print "The corresponding p-value is:", float(errH - errA)/np.sqrt(vard),df, ",and, ", t.cdf(float(errH - errA)/np.sqrt(vard),df)
    print "The corresponding statistic and p-value is:", float(errH - errA) / np.sqrt(vard), df, ",and, ", 1-t.sf(float(errH - errA) / np.sqrt(vard), df)
    if lower_value <= 0 <= upper_value:
        print("No statistical significant")

    evObj.pValue = 1-t.sf(float(errH - errA) / np.sqrt(vard), df)
    evObj.interval = [lower_value, upper_value]
    evObj.f1_avg = {'f1C1': f1C1, 'f1C2':f1C2, 'f1_avg data':f1_avg, 'f1_H0':errH }
    v.stat_test_result[str(evObj.event)] = [upper_value,lower_value]

    #
    # print "The resulrs:\n",
    # pprint.pprint(v.stat_test_result)
    #
    # print "The following events follow H0:\n"
    # tr1= v.stat_test_result
    # pprint.pprint([ key for key in tr1 if tr1[key][0]*tr1[key][1] <0])
    #
    # print "The following events follow H1:\n"
    # pprint.pprint([key for key in tr1 if tr1[key][0] * tr1[key][1] > 0])
