This libaray is desgined to extract mutual fingerprints from process variant. It is based on Python 2.7.

#Instruction:
1) Download the whole files in a directory
2) Open Python or Ipython terminal in that directory
3) Import Main as m
4) Specify the paths of process variants by variables path1 and path2
5) Call m.Main(path1,path2)
6) Wait until it comes up with result. The finger prints are pop up as pdf files and are saved in the same directory of codes
########################################################################
7) If you want to access the deviant parts for further analysis you can have all the results in module Variables.py
   For example, to print the discriminatory edges you can run the following codes:

	import Variables as v
	for u in v.deviant_obj:
    	...:     if u.event[0] in v.deviant_element_plot:
    	...:         print u.event
    	...:         print u.f1_avg
    	...:         print u.class_proportion
    	...:         print u.pValue
    	...:         print "-----------------------------"

#######################################################################
If you got any problems, feel free to send me a massage. I try to answer you as early as possible!

Enjoy!