#This variable is just for system printing issue. Related to Print_Screen_Off(exp)
sys_stdout_old=None

#----------------------------
#Storing a dictionary for the elements and their information that will be showed on the process model
#This variable belongs to Visualization.Log_Filter2()
deviant_element_plot={}
deviant_element_timeStamp_plot={}

#----------------------------
#This variable stores objects that contain events which are deviants in terms of timestamps.
#See Deviation.Deviation_Timestamp()
deviant_obj_timestamp=[]

#---------------------------

m_reduced=[]
m2_reduced=[]

stat_test_result={}


#It contains the list of deviant objects. Each object has several attributes. See Deviation.py
deviant_obj=[]
#----------------------------
#The following variables will be used for waiting time analysis
#The logs are row, which means that they contains duplicate traces.
log1_row=[]
log1_timestamps=[]

log2_row=[]
log2_timestamps=[]

#----------------------------
#Dictionaries that store unique traces with frequnecies
#like {0: '0,1,2,3,4,5,6,7,8,9',
#       1: '10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49'}
id1={}
id2={}

#Storing unique logs
log1=[]
log2=[]