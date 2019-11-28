import numpy as np
from scipy import stats
import Preparation as pr
import pprint
import operator
import Variables as v
from nltk.util import ngrams
import collections
from scipy.stats import chisquare, chi2_contingency, t
from scipy.special import expit
from sklearn.metrics.cluster import contingency_matrix
from sklearn.svm import SVC
from sklearn.metrics import f1_score
import gc
import copy
import tqdm


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


    target_predicted = collections.Counter(u);
    print "Target predicted:", target_predicted
    target_true = collections.Counter(Y);
    print "Target true values:", target_true

    print "Contingency matrix is:", contingency_matrix(Y, u, eps=None, sparse=False)
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
#Implementing t-test for the performance of the classifier
#https://www3.nd.edu/~rjohns15/cse40647.sp14/www/content/lectures/28%20-%20Classifier%20Comparisons.pdf
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
