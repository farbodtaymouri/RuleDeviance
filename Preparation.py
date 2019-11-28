import math
import numpy as np
from numpy.linalg import inv
from collections import Counter
import operator




def WaveletCompute(s):
    '''

    :param s: a vector list of 0,1 of power 2
    :return: Returning Haar wavelet coefficients
    '''

    #check the length of input is power of two
    p = math.log(len(s) ,2)
    if not ((math.log(len(s) ,2)  - int(math.log(len(s) ,2)) ) == 0):
        raise ValueError('The length of input is not power of 2!')


    #Storing coefficients
    coeffs=[]
    coeffs.append(np.average(s))

    #Traversing each level
    for  level in range(1,int(p)+1):
        #Size of each chunk
        chunk_length = int(len(s)/math.pow(2,level))
        print("chunk lenth:", chunk_length)

        temp = []
        for i in range(int(math.pow(2,level))):
            print('ssss',chunk_length*i,(i+1)*chunk_length )
            chunk = s[chunk_length*i : (i+1)*chunk_length]
            print("the chanuck is:", chunk)
            temp.append(chunk)

            if(len(temp) == 2):
                coeffs.append ( round(( np.average(temp[0]) - np.average(temp[1]) )/2 ,3))
                temp=[]

    #Returning coefficients
    return coeffs


class Transformation:
#This class given a sequence like s=[1,2,3,4,5,6], provides the transformation matrix that is requried to
#transfomr the sequence into a vector
    def __init__(self):
        self.matrix_transform = []
        self.matrix_transform_inv = []


    #-----------------------------------------------------------------
    #Matrix form of wavelet
    def WaveletCf(self, s):

        '''matrix_transform = np.matrix([[1,1,1,1,1],
                                     [1,1,1, -float(3)/2, -float(3)/2 ],
                                     [1,1,-2,0,0],
                                     [1,-1,0,0,0],
                                     [0,0,0,1,-1]])

        s=[1,2,3,4,5]


        wavelet_coefficients = np.dot(s, inv(matrix_transform))
        input_recover = np.round( np.dot(wavelet_coefficients,matrix_transform), 3)

        print"the original input is:", s
        print "The coefficients are:", wavelet_coefficients
        print "The recovered input is:", input_recover'''


        #First computing the transform matrix (which provides class variables 'matrix_transform')
        self.__Matrix_tr_Compute(s)

        '''wavelet_coefficients = np.round(np.dot(s, inv(self.matrix_transform)),3)
        input_recover = np.round(np.dot(wavelet_coefficients, self.matrix_transform), 3)

        print"the original input is:", s
        print "The coefficients are:", wavelet_coefficients
        print "The recovered input is:", input_recover'''





    #--------------------------------------------------

    #This method creates the matrix transformation given an input string 's'
    def __Matrix_tr_Compute(self,s):
        's is a list like s=[1,2,3,4,5]'
        #matrix=[]
        #matrix.append([1] * len(s))
        #Initializing
        self.matrix_transform=[]
        self.matrix_transform.append(np.array([1] * len(s),dtype='f'))


        lower=0
        upper=len(s)-1



        self.__Recursive(s, lower, upper)
        self.matrix_transform_inv = np.round(inv(self.matrix_transform),3)


    def __Recursive(self, s, lower, upper):

        #When we reach two consecutive elements
        if( (upper - lower) ==1):
            temp=np.array([0] * len(s),dtype='f')
            temp[upper] =-1; temp[lower]=1
            #print "The lower, upper and the vector is:", lower, upper, temp
            self.matrix_transform.append(temp)
            #return temp

        elif ((upper - lower) > 1):
            temp = np.array([0] * len(s),dtype='f')
            middle = int(math.floor( (upper + lower)/float(2)))
            temp[lower:(middle+1)] =1
            temp[(middle+1):upper+1] = -1*(float(middle - lower +1)/(upper-middle))

            #print "The lower, middle, upper and the vector is:", lower, middle, upper, temp, float(middle - lower +1)/(upper-middle)
            self.matrix_transform.append(temp)

            self.__Recursive(s, lower, middle)
            self.__Recursive(s, middle+1, upper)

            #return temp
        else:
            pass
#-----------------------------------------------------------

class Seq2Vec (Transformation):

    def __init__(self):
        self.dic_BinaryVector = {}  # Its like {'a': [0,0,0,1,0],....}
        self.dic_WaveletCfVector= {}  # Its like {'a': [.5,-1,0,1,0],....}
        self.unique_Events_log = []   #Its a list contains all the unique events of the log





    def Start(self,s):

        # Computing a dictionary of binary vectors
        self.__BinaryVector(s)

        #Computing coefficient vectors
        #computing the matrix transform
        self.WaveletCf(s)

        #computing coefficients
        for e in self.dic_BinaryVector:
            self.dic_WaveletCfVector[e] = np.round(np.dot(self.dic_BinaryVector[e], inv(self.matrix_transform)), 3).tolist()


    #------------------------------------------------------------
    def __BinaryVector(self, s):
        self.dic_BinaryVector ={}


        max_length = len(s)

        #Finding the number of unique elements in the input trace
        ##events = Counter(s).keys()
        events = self.unique_Events_log


        #Finding all the indexes of each unique event
        for e in events:
            #This case happens when we work with Kgram(k>1)
            if(len(e)>1):
                temp = np.array([0] * max_length)
                for  i in range(len(s)):
                    if(s[i]==e):
                        temp[i] = 1


                self.dic_BinaryVector[e] = list(temp)

            else:
                inds = np.where(np.array(s) == e)[0]       #It returns like (array([1, 2, 3], dtype=int64),)
                temp = np.array([0]*max_length)
                temp[inds] =1
                self.dic_BinaryVector[e] = list(temp)

        print "The dic is:", self.dic_BinaryVector

    #------------------------------------------

    def Log_Normalization(self,log):
        #Flatting the log
        temp = reduce(operator.concat, log)

        # Finding the number of unique elements in the input trace
        unique_events = Counter(temp).keys()
        self.unique_Events_log = unique_events

        # Finding the maximum length
        max_length = max([len(trace) for trace in log])


        #Adding slack elements to the traces to meke them equal length "<T>"
        for i in range(len(log)):
            if(len(log[i]) < max_length):
                log[i]+= ['<T>'] * (max_length - len(log[i]))

        '''print "The unique elemensts are:", unique_events
        print "The maximum len is:", max_length
        print "The augmented log is:", log'''

        return log


#-----------------------------------------------------------











'''m=[[1,-1,0,0,0,0,0,0],
   [0,0,1,-1,0,0,0,0],
   [0,0,0,0,1,-1,0,0],
   [0,0,0,0,0,0,1,-1],
   [1,1,-1,-1,0,0,0,0],
   [0,0,0,0,1,1,-1,-1],
   [1,1,1,1,-1,-1,-1,-1],
   [1,1,1,1,1,1,1,1]]
cf = [1, -0.5, -1, 0.5, 2.25, -0.25, -0.25, 5]'''


