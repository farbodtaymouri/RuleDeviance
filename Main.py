import os
import sys
import time
import tqdm
import pprint
import Reading_Wriitng as rw
import DimReduction as dr
import Visualization as vs
import Variables as v
import Deviation as dv



def Print_Screen_Off(exp):
    if exp =='off':
        # -----------This part of code is to prevent printing on the screen!!!--------
        if (v.sys_stdout_old == None):
            v.sys_stdout_old =  sys.stdout

        sys.stdout = open(os.devnull, "w")
        f = open(os.devnull, 'w')
        sys.stdout = f
    else:
        sys.stdout = v.sys_stdout_old
# ----------------------------------------------------------------------------

############################################################################
#This version of Main is created on 11/09/2109
def Main(path1, path2):
    '''

    :param path1: Location of the first event log
    :param path2: Location of the second event log
    :return:
    '''

    Print_Screen_Off('off')


    # log1 = [['c', 'b','a','f','a']] * 20  + [['b', 'a', 'c']] * 10
    # log2 = [['a', 'e','f','d']] * 10+ [['b', 'a', 'c']] *10
    #
    #
    # log1, id1 = dv.Log_dic(log1)
    # log2, id2 = dv.Log_dic(log2)




    path1= 'C:/Users/taymouri/Desktop/model temp/trash2/logs/BPIC13_incidents_orgline_A2.xes'
    path2 ='C:/Users/taymouri/Desktop/model temp/trash2/logs/BPIC13_incidents_orgline_C.xes'



    d=" "
    log1, id1 = rw.Reading_log_all(path1, d)
    log2, id2 = rw.Reading_log_all(path2, d)




    #Storing the id of unique traces and unique logs
    v.id1 = id1
    v.id2 = id2
    v.log1=log1[:]
    v.log2=log2[:]

    # #Log analysing
    # log1,log2, id1, id2 = dr.Log_Stat()
    # #return



    #Computing Ngram
    kgram=2
    log1 = dr.Ngram_Compute(log1,kgram)
    log2 = dr.Ngram_Compute(log2,kgram)
    print "Bigram1:", pprint.pprint(log1)
    print "Bigram2:",pprint.pprint(log2)





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


        #Finding at which level (coefficient) the behavior of event is different between classes
        #evObj.Level_difference(event_matrix1,event_matrix2)


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



