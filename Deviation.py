import numpy as np
from scipy.stats import ttest_ind
import xml.etree.ElementTree as ET
import re
from dateutil import parser
import DimReduction as dr
import Variables as v


class Deviant_Element:

    def __init__(self, event):
        self.event = event  #This shows the statistically significant element like ['a','b'], it is kind of k-gram
        self.level = {}    #This shoes the level of deviation (it is the wavelet coefficient that deviation happens + t-statistics), like {1:-2,2:0.5,..}
        self.pValue=''
        self.interval=''
        self.design_matrix=''    #It shows the corresponding matrix of coefficients, where the last column is the labels
        self.f1_avg = None
        self.class_proportion = None


    def Level_difference(self,matrix1,matrix2):
        '''
        Computing at which level (i.e., coefficient) the behavior are statistically different
        This procedure only applies to localize the deviation. To check that the behavior is totally different SVM already used
        :param matrix1: like :[[0.2, 0.133, 0.167, 0.5, 0.0],
                            [0.2, 0.133, 0.167, 0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0],
                             [0.2, 0.133, 0.167, -0.5, 0.0]]
            matrix2: Similar oto the above matrix
        :return: a list shows the corresponding levels [0,2,5]
        '''

        if(self.pValue > 0.05):
            #There is nothing to be localized
            return

        # Some features only happened in one class, therefore no need to run a classifier for them
        # like matrix1=[[],...[]]
        matrix1 = [t for t in matrix1 if t != []]
        matrix2 = [t for t in matrix2 if t != []]


        matrix1 = np.array(matrix1); matrix2=np.array(matrix2)
        for col in range(matrix1.shape[1]):
            #This situation t-test does not work since one of the columns are zero or constant
            if( (np.std(matrix1[:,col]) == 0) or (np.std(matrix2[:,col]) == 0)):
                if(np.mean(matrix1[:,col]) != np.mean(matrix2[:,col]) ):
                    self.level[col] = np.mean(matrix1[:,col])  - np.mean(matrix2[:,col])

            else:

                #Checking whether the mean of coefficients are different across the same column
                if(ttest_ind(matrix1[:,col], matrix2[:,col],equal_var=False)[1] < 0.05):

                    #Storing the t-statistics for the significant column
                    self.level[col]= ttest_ind(matrix1[:,col], matrix2[:,col],equal_var=False)[0]





#----------------------------------------------------------------------

#This a function that returns unique traces + corresponding ids
def Log_dic(log):
    '''

    :param log: = [['c', 'b','a','f','a']] * 10  + [['b', 'a', 'c']] * 40
    :return: {0: '0,1,2,3,4,5,6,7,8,9',
              1: '10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49'}
    '''
    dic_id={}
    temp =[]
    k=0
    for i in range(len(log)):
        if log[i] not in temp:
            dic_id[k]=str(i)
            temp.append(log[i])
            k+=1
        else:
            dic_id[temp.index(log[i])]=dic_id[temp.index(log[i])] +','+str(i)

    return temp, dic_id

#----------------------------------------------------------------------

def Deviation_Timestamp(path,path2):
    def Read_XES(path):
        add = path
        tree = ET.parse(add)
        root = tree.getroot()

        # --New Added 30 August 2016-------------------------------------------
        # Some XML files have Namesapce and it is declered at the root like and represented by xmlns like: <pnml xmlns="http://www.pnml.org/version-2009/grammar/pnml">.
        # So in order to acces the tags inside the XML file, the name space must be mentioned!, The general form is like "{name_sapce}tag",
        # For example for reading places tags, the code is like below:
        #   for p in root.iter("{ http://www.pnml.org/version-2009/grammar/pnml }place"):
        #       print p
        # ------
        # First we need to extract the namespace, namely the value of xmlns
        # root.tag='{http://www.pnml.org/version-2009/grammar/pnml}pnml' but we need '{http://www.pnml.org/version-2009/grammar/pnml}'
        # For this issue we use regular expression library (re) of python
        m = re.match('\{.*\}', root.tag)
        # checking whether m is empty or no
        if (m):
            namespace = m.group(0)
        else:
            namespace = ''
            # ------------------------------------------------------------

        temp2 = []
        temp = []
        log = []
        case_name = []
        event_time_stamp = []
        print "We are before for"
        for t in root.iter(namespace + "trace"):
            # ----
            # Reading the case number, like "instance_294"
            for name in t.iter(namespace + "string"):
                # print name.attrib['key']
                if (name.attrib['key'] == 'concept:name'):
                    try:
                        temp_name = name.attrib['value']
                        temp_name = temp_name.split("_")[1]
                        case_name.append(int(temp_name))
                    except IndexError:
                        # When the case number is like A1001
                        temp_name = temp_name.split("_")[0]
                        case_name.append(temp_name)

                    break
            # print "The caee_name:",case_name
            # ---------------------
            for e in t.iter(namespace + "event"):
                for r in e.iter(namespace + "string"):
                    # print r.attrib
                    if (r.attrib['key'] == 'concept:name'):
                        # print r.attrib['value']
                        temp.append(r.attrib['value'])

                # New 14/10/2019
                # Searcing to the time of events
                for r in e.iter(namespace + 'date'):
                    if (r.attrib['key'] == "time:timestamp"):
                        temp2.append(r.attrib['value'])
            event_time_stamp.append(temp2)
            temp2 = []

            # print "the temp is:", temp
            log.append(temp)
            temp = []

            # print t.attrib

        print "We are after for"
        print "The len log is:", len(log)
        print "the len case_name:", len(case_name)
        ###-------
        # re-arranging the position of the log, starting from zero
        # Before arranging case_name[0]='297' and log[0]= case 297
        # After  arranging, case_name[0]='0' and log[0]=case 0
        # Example:

        #    >>> list1 = [3,2,4,1, 1]
        #    >>> list2 = ['three', 'two', 'four', 'one', 'one2']
        #    >>> list1, list2 = zip(*sorted(zip(list1, list2)))
        #    >>> list1
        #        (1, 1, 2, 3, 4)
        #    >>> list2
        #        ('one', 'one2', 'two', 'three', 'four')

        print "The len log is:", log
        print "the len case_name:", case_name
        ##case_name, log=zip(*sorted(zip(case_name,log)))

        # ---------------------------------------

        ##printting the mean length of the traces
        men = 0
        for i in range(len(log)):
            men += len(log[i])
        print "the mean lenght of traces:", men / len(log)

        #return log, case_name
        return log, case_name, event_time_stamp

    ##########################################
    #Reading the log and timestamps
    log1, case_name1, event_time_stamp1 = Read_XES(path)
    log2, case_name2, event_time_stamp2 = Read_XES(path2)


    '''
    log is like this: [ ['ER Registration',
                         'ER Triage',
                         'ER Sepsis Triage',
                         'IV Liquid',
                         'CRP',
                         'Leucocytes',
                         'LacticAcid',
                         'IV Antibiotics',
                         'Admission NC',
                         'Release A',
                         'Return ER'],.....]
                         
    event_time_stamp is like [ ['2014-10-12T07:29:02.000+11:00',
                                 '2014-10-12T07:38:13.000+11:00',
                                 '2014-10-12T07:38:37.000+11:00',
                                 '2014-10-12T07:43:13.000+11:00',
                                 '2014-10-12T07:51:00.000+11:00',
                                 '2014-10-12T07:51:00.000+11:00',
                                 '2014-10-12T07:51:00.000+11:00',
                                 '2014-10-12T08:00:27.000+11:00',],......]
    '''


    #Computing bigrams
    log1 = dr.Ngram_Compute(log1,2)
    event_time_stamp1 = dr.Ngram_Compute(event_time_stamp1,2)

    log2 = dr.Ngram_Compute(log2, 2)
    event_time_stamp2 = dr.Ngram_Compute(event_time_stamp2, 2)

    # --------------------------------
    #Computing difference time for the bigram

    #First log
    time_difference1=[]
    for i in range(len(log1)):
        temp=[]
        for k in list(event_time_stamp1[i]):
            delta = parser.parse(k[1]) - parser.parse(k[0])
            temp.append(delta.total_seconds())
        time_difference1.append(temp)

    #Second log
    time_difference2=[]
    for i in range(len(log2)):
        temp=[]
        for k in list(event_time_stamp2[i]):
            delta = parser.parse(k[1]) - parser.parse(k[0])
            temp.append(delta.total_seconds())
        time_difference2.append(temp)

    #--------------------------------
    #Creating a dictionary of the form {('ER Registration', 'ER Triage'):[12,34,21,45], ('IV Liquid', 'CRP'):[4,55,12],......}

    event_diffTime_dic1={}
    for i in range(len(log1)):
        for j in range(len(log1[i])):
            if(log1[i][j] in event_diffTime_dic1):
                event_diffTime_dic1[log1[i][j]].append(time_difference1[i][j])
            else:
                event_diffTime_dic1[log1[i][j]] = [time_difference1[i][j]]


    #Creating a dictionary of the form {('ER Registration', 'ER Triage'):[12,34,21,45], ('IV Liquid', 'CRP'):[4,55,12],......}
    event_diffTime_dic2={}
    for i in range(len(log2)):
        for j in range(len(log2[i])):
            if(log2[i][j] in event_diffTime_dic2):
                event_diffTime_dic2[log2[i][j]].append(time_difference2[i][j])
            else:
                event_diffTime_dic2[log2[i][j]] = [time_difference2[i][j]]

    # --------------------------------
    deviant_obj_timestamp=[]
    #Computing T-test between common elements
    for key in set(event_diffTime_dic1.keys()).intersection(set(event_diffTime_dic2.keys())):
        pValue = ttest_ind(event_diffTime_dic2[key], event_diffTime_dic1[key])[1]
        if(pValue <0.05):
            #Creating an object to save deviant elements (timstamp deviation)
            #obj = Deviant_Element(key)
            obj = Deviant_Element([tuple([e for e in key])])
            #obj = Deviant_Element([tuple([e.replace(' ', '#') for e in key])] )      #It would be like [('Payment', 'Payment')]
            obj.pValue = pValue
            obj.avgTimestamp = {'log1': np.average( event_diffTime_dic1[key]), 'log2':np.average( event_diffTime_dic2[key]) }
            obj.class_proportion = {"no_matrix1": len(event_diffTime_dic1[key]), "no_matrix2": len(event_diffTime_dic2[key]),
                                      "proportion": float(len(event_diffTime_dic1[key]))/len(event_diffTime_dic2[key])}

            deviant_obj_timestamp.append(obj)


    #Storing the variables
    v.deviant_obj_timestamp = deviant_obj_timestamp


    #return log1, case_name1, event_time_stamp1, time_difference1, event_diffTime_dic1


