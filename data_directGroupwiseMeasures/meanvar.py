# -*- coding: utf-8 -*-
import numpy as np
import pylab as py
from struct import *

#Variables______________________________________________________________________
gM_direct_length = 18
geneNum = 5000
gNgN = geneNum*geneNum
annotSizesNum = 5
gN_div_aN = geneNum/annotSizesNum

annots = np.arange(annotSizesNum) 
geneProducts = np.arange(1,geneNum + 1)

sim = np.array([np.zeros(5000*5000),np.zeros(5000*5000),np.zeros(5000*5000),np.zeros(5000*5000),np.zeros(5000*5000),np.zeros(5000*5000),
                np.zeros(5000*5000),np.zeros(5000*5000),np.zeros(5000*5000),np.zeros(5000*5000),np.zeros(5000*5000),np.zeros(5000*5000),
                np.zeros(5000*5000),np.zeros(5000*5000),np.zeros(5000*5000),np.zeros(5000*5000),np.zeros(5000*5000),np.zeros(5000*5000)])

#Average similarities of gene products with annots of a certain size to gene products of size 1,10,50,100,1000___________________________              
sim_avg_annots1 = np.array([np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                            np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                            np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5)])    
sim_avg_annots10 = np.array([np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                             np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                             np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5)])  
sim_avg_annots50 = np.array([np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                             np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                             np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5)])                     
sim_avg_annots100 = np.array([np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                              np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                              np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5)])  
sim_avg_annots1000 = np.array([np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                               np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                               np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5)])  
							   
#Variance in similarities of gene products with annots of a certain size to gene products of size 1,10,50,100,1000___________________________               
sim_var_annots1 = np.array([np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                            np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                            np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5)])    
sim_var_annots10 = np.array([np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                             np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                             np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5)])  
sim_var_annots50 = np.array([np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                             np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                             np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5)])                     
sim_var_annots100 = np.array([np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                              np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                              np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5)])  
sim_var_annots1000 = np.array([np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                               np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),
                               np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5)])  
                         
#Open binary data files with similarities.___________________________________________________________________________________________________________________________
data = np.array([open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_1_SIM_FRAMEWORK_DAG_SET_BADER_2003.bin", "rb"),    
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_2_SIM_FRAMEWORK_DAG_SET_BRAUN_BLANQUET_1932.bin", "rb"), 
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_3_SIM_FRAMEWORK_DAG_SET_DICE_1945.bin", "rb"),  
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_4_SIM_FRAMEWORK_DAG_SET_JACCARD_1901.bin", "rb"),   
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_5_SIM_FRAMEWORK_DAG_SET_KNAPPE_2004.bin", "rb"),  
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_6_SIM_FRAMEWORK_DAG_SET_KORBEL_2002.bin", "rb"),   
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_7_SIM_FRAMEWORK_DAG_SET_MARYLAND_BRIDGE_2003.bin", "rb"),  
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_8_SIM_FRAMEWORK_DAG_SET_OCHIAI_1957.bin", "rb"),  
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_9_SIM_FRAMEWORK_DAG_SET_SIMPSON_1960.bin", "rb"),   
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_10_SIM_FRAMEWORK_DAG_SET_SOKAL_SNEATH_1963.bin", "rb"),   
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_11_SIM_FRAMEWORK_DAG_SET_TVERSKY_1977.bin", "rb"),  
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_12_SIM_GROUPWISE_DAG_ALI_DEANE.bin", "rb"),  
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_13_SIM_GROUPWISE_DAG_NTO.bin", "rb"),  
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_14_SIM_GROUPWISE_DAG_NTO_MAX.bin", "rb"),   
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_15_SIM_GROUPWISE_DAG_UI.bin", "rb"), 
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_16_SIM_GROUPWISE_DAG_GIC.bin", "rb"), 
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_17_SIM_GROUPWISE_DAG_GIC.bin", "rb"),
                 open("C:/Users/#ballin/Desktop/individualMeasures_data/directGroupwiseMeasures/similarities_18_SIM_GROUPWISE_DAG_GIC.bin", "rb")])

#Fill similarities array________________________________________
for i in xrange(0,gM_direct_length):
    data[i].seek(0)
    sim[i] = np.fromfile(data[i], dtype='>f8', count=-1, sep='')  
    data[i].close()

#Average similarity to annotations of size 1____________________________________
#_______________________________________________________________________________
for i in xrange (0, gM_direct_length): 
    arr = np.zeros(1000)   
    for j in xrange(0,1000):
        sim_avg_annots1[i][0] += (1.0/1000.0)*np.sum(sim[i][j*geneNum : j*geneNum + 1000])
		arr[j] = (1.0/1000.0)*np.sum(sim[i][j*geneNum : j*geneNum + 1000])
	sim_avg_annots1[i][0] = (1.0/1000.0)*sim_avg_annots1[i][0] 
	for j in xrange(0,1000):
        arr[j] = (arr[j] - sim_avg_annots1[i][0])*(arr[j] - sim_avg_annots1[i][0])
    sim_var_annots1[i][0] = (1.0/1000.0)*np.sum(arr)
    print i    
for i in xrange (0, gM_direct_length): 
    arr = np.zeros(1000)   
    for j in xrange(1000,2000):
        sim_avg_annots1[i][1] += (1.0/1000.0)*np.sum(sim[i][j*geneNum : j*geneNum + 1000])
		arr[j - 1000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum : j*geneNum + 1000])
    sim_avg_annots1[i][1] = (1.0/1000.0)*sim_avg_annots1[i][1] 
	for j in xrange(0,1000):
        arr[j] = (arr[j] - sim_avg_annots1[i][1])*(arr[j] - sim_avg_annots1[i][1])
    sim_var_annots1[i][1] = (1.0/1000.0)*np.sum(arr)
    print i
for i in xrange (0, gM_direct_length): 
    arr = np.zeros(1000) 
    for j in xrange(2000,3000):
        sim_avg_annots1[i][2] += (1.0/1000.0)*np.sum(sim[i][j*geneNum : j*geneNum + 1000])
		arr[j - 2000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum : j*geneNum + 1000])
    sim_avg_annots1[i][2] = (1.0/1000.0)*sim_avg_annots1[i][2] 
	for j in xrange(0,1000):
        arr[j] = (arr[j] - sim_avg_annots1[i][2])*(arr[j] - sim_avg_annots1[i][2])
    sim_var_annots1[i][2] = (1.0/1000.0)*np.sum(arr)
    print i
for i in xrange (0, gM_direct_length):   
	arr = np.zeros(1000) 
    for j in xrange(3000,4000):
        sim_avg_annots1[i][3] += (1.0/1000.0)*np.sum(sim[i][j*geneNum : j*geneNum + 1000])
		arr[j - 3000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum : j*geneNum + 1000])
    sim_avg_annots1[i][3] = (1.0/1000.0)*sim_avg_annots1[i][3] 
	for j in xrange(0,1000):
        arr[j] = (arr[j] - sim_avg_annots1[i][3])*(arr[j] - sim_avg_annots1[i][3])
    sim_var_annots1[i][3] = (1.0/1000.0)*np.sum(arr)
    print i   
for i in xrange (0, gM_direct_length):  
	arr = np.zeros(1000) 
    for j in xrange(4000,geneNum):
        sim_avg_annots1[i][4] += (1.0/1000.0)*np.sum(sim[i][j*geneNum : j*geneNum + 1000])
		arr[j - 4000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum : j*geneNum + 1000])
    sim_avg_annots1[i][4] = (1.0/1000.0)*sim_avg_annots1[i][4] 
	for j in xrange(0,1000):
        arr[j] = (arr[j] - sim_avg_annots1[i][4])*(arr[j] - sim_avg_annots1[i][4])
    sim_var_annots1[i][4] = (1.0/1000.0)*np.sum(arr)	
    print i
    
#Average similarity to annotations of size 10___________________________________
#_______________________________________________________________________________
for i in xrange (0, gM_direct_length):  
	arr = np.zeros(1000)
    for j in xrange(0,1000):
        sim_avg_annots10[i][0] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 1000: j*geneNum + 2000])
		arr[j] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 1000: j*geneNum + 2000])
    sim_avg_annots10[i][0] = (1.0/1000.0)*sim_avg_annots10[i][0] 
	for j in xrange(0,1000):
        arr[j] = (arr[j] - sim_avg_annots10[i][0])*(arr[j] - sim_avg_annots10[i][0])
    sim_var_annots10[i][0] = (1.0/1000.0)*np.sum(arr)		
    print i    
for i in xrange (0, gM_direct_length):  
    arr = np.zeros(1000)
    for j in xrange(1000,2000):
        sim_avg_annots10[i][1] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 1000: j*geneNum + 2000])
		arr[j - 1000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 1000: j*geneNum + 2000])
    sim_avg_annots10[i][1] = (1.0/1000.0)*sim_avg_annots10[i][1] 
	for j in xrange(0,1000):
        arr[j] = (arr[j] - sim_avg_annots10[i][1])*(arr[j] - sim_avg_annots10[i][1])
    sim_var_annots10[i][1] = (1.0/1000.0)*np.sum(arr)	
    print i
for i in xrange (0, gM_direct_length):  
	arr = np.zeros(1000)
    for j in xrange(2000,3000):
        sim_avg_annots10[i][2] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 1000: j*geneNum + 2000])
		arr[j - 2000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 1000: j*geneNum + 2000])
    sim_avg_annots10[i][2] = (1.0/1000.0)*sim_avg_annots10[i][2] 
	for j in xrange(0,1000):
        arr[j] = (arr[j] - sim_avg_annots10[i][2])*(arr[j] - sim_avg_annots10[i][2])
    sim_var_annots10[i][2] = (1.0/1000.0)*np.sum(arr)	
    print i
for i in xrange (0, gM_direct_length):    
	arr = np.zeros(1000)
    for j in xrange(3000,4000):
        sim_avg_annots10[i][3] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 1000: j*geneNum + 2000])
		arr[j - 3000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 1000: j*geneNum + 2000])
    sim_avg_annots10[i][3] = (1.0/1000.0)*sim_avg_annots10[i][3] 
	for j in xrange(0,1000):
        arr[j] = (arr[j] - sim_avg_annots10[i][3])*(arr[j] - sim_avg_annots10[i][3])
    sim_var_annots10[i][3] = (1.0/1000.0)*np.sum(arr)	
    print i   
for i in xrange (0, gM_direct_length): 
	arr = np.zeros(1000)
    for j in xrange(4000,geneNum):
        sim_avg_annots10[i][4] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 1000: j*geneNum + 2000])
		arr[j - 4000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 1000: j*geneNum + 2000])
    sim_avg_annots10[i][4] = (1.0/1000.0)*sim_avg_annots10[i][4] 
	for j in xrange(0,1000):
        arr[j] = (arr[j] - sim_avg_annots10[i][4])*(arr[j] - sim_avg_annots10[i][4])
    sim_var_annots10[i][4] = (1.0/1000.0)*np.sum(arr)	
    print i
    
#Average similarity to annotations of size 50___________________________________
#_______________________________________________________________________________
for i in xrange (0, gM_direct_length): 
    arr = np.zeros(1000)   
    for j in xrange(0,1000):
        sim_avg_annots50[i][0] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 2000: j*geneNum + 3000])
		arr[j] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 2000: j*geneNum + 3000])
    sim_avg_annots50[i][0] = (1.0/1000.0)*sim_avg_annots50[i][0] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots50[i][0])*(arr[j] - sim_avg_annots50[i][0])
	sim_var_annots50[i][0] = (1.0/1000.0)*np.sum(arr)	
    print i    
for i in xrange (0, gM_direct_length):
	arr = np.zeros(1000)
    for j in xrange(1000,2000):
        sim_avg_annots50[i][1] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 2000: j*geneNum + 3000])
		arr[j - 1000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 2000: j*geneNum + 3000])
    sim_avg_annots50[i][1] = (1.0/1000.0)*sim_avg_annots50[i][1] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots50[i][1])*(arr[j] - sim_avg_annots50[i][1])
	sim_var_annots50[i][1] = (1.0/1000.0)*np.sum(arr)	
    print i
for i in xrange (0, gM_direct_length):
	arr = np.zeros(1000)
    for j in xrange(2000,3000):
        sim_avg_annots50[i][2] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 2000: j*geneNum + 3000])
		arr[j - 2000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 2000: j*geneNum + 3000])
    sim_avg_annots50[i][2] = (1.0/1000.0)*sim_avg_annots50[i][2] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots50[i][2)*(arr[j] - sim_avg_annots50[i][2])
	sim_var_annots50[i][2] = (1.0/1000.0)*np.sum(arr)	
    print i
for i in xrange (0, gM_direct_length):
	arr = np.zeros(1000)
    for j in xrange(3000,4000):
        sim_avg_annots50[i][3] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 2000: j*geneNum + 3000])
		arr[j - 3000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 2000: j*geneNum + 3000])
    sim_avg_annots50[i][3] = (1.0/1000.0)*sim_avg_annots50[i][3] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots50[i][3])*(arr[j] - sim_avg_annots50[i][3])
	sim_var_annots50[i][3] = (1.0/1000.0)*np.sum(arr)	
    print i   
for i in xrange (0, gM_direct_length):   
	arr = np.zeros(1000)
    for j in xrange(4000,geneNum):
        sim_avg_annots50[i][4] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 2000: j*geneNum + 3000])
		arr[j - 4000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 2000: j*geneNum + 3000])
    sim_avg_annots50[i][4] = (1.0/1000.0)*sim_avg_annots50[i][4] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots50[i][4])*(arr[j] - sim_avg_annots50[i][4])
	sim_var_annots50[i][4] = (1.0/1000.0)*np.sum(arr)	
    print i
    
#Average similarity to annotations of size 100___________________________________
#_______________________________________________________________________________
for i in xrange (0, gM_direct_length): 
	arr = np.zeros(1000)
    for j in xrange(0,1000):
        sim_avg_annots100[i][0] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 3000: j*geneNum + 4000])
		arr[j] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 3000: j*geneNum + 4000])
    sim_avg_annots100[i][0] = (1.0/1000.0)*sim_avg_annots100[i][0] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots100[i][0])*(arr[j] - sim_avg_annots100[i][0])
	sim_var_annots100[i][0] = (1.0/1000.0)*np.sum(arr)	
    print i    
for i in xrange (0, gM_direct_length):  
	arr = np.zeros(1000)
    for j in xrange(1000,2000):
        sim_avg_annots100[i][1] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 3000: j*geneNum + 4000])
		arr[j - 1000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 3000: j*geneNum + 4000])
    sim_avg_annots100[i][1] = (1.0/1000.0)*sim_avg_annots100[i][1] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots100[i][1])*(arr[j] - sim_avg_annots100[i][1])
	sim_var_annots100[i][1] = (1.0/1000.0)*np.sum(arr)	
    print i
for i in xrange (0, gM_direct_length):   
	arr = np.zeros(1000)
    for j in xrange(2000,3000):
        sim_avg_annots100[i][2] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 3000: j*geneNum + 4000])
		arr[j - 2000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 3000: j*geneNum + 4000])
    sim_avg_annots100[i][2] = (1.0/1000.0)*sim_avg_annots100[i][2] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots100[i][2])*(arr[j] - sim_avg_annots100[i][2])
	sim_var_annots5100[i][2] = (1.0/1000.0)*np.sum(arr)	
    print i
for i in xrange (0, gM_direct_length): 
	arr = np.zeros(1000)
    for j in xrange(3000,4000):
        sim_avg_annots100[i][3] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 3000: j*geneNum + 4000])
		arr[j - 3000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 3000: j*geneNum + 4000])
    sim_avg_annots100[i][3] = (1.0/1000.0)*sim_avg_annots100[i][3] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots100[i][3])*(arr[j] - sim_avg_annots100[i][3])
	sim_var_annots100[i][3] = (1.0/1000.0)*np.sum(arr)	
    print i   
for i in xrange (0, gM_direct_length):  
	arr = np.zeros(1000)
    for j in xrange(4000,geneNum):
        sim_avg_annots100[i][4] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 3000: j*geneNum + 4000])
		arr[j - 4000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 3000: j*geneNum + 4000])
    sim_avg_annots100[i][4] = (1.0/1000.0)*sim_avg_annots100[i][4] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots100[i][4])*(arr[j] - sim_avg_annots100[i][4])
	sim_var_annots100[i][4] = (1.0/1000.0)*np.sum(arr)	
    print i
    
#Average similarity to annotations of size 1000___________________________________
#_______________________________________________________________________________
for i in xrange (0, gM_direct_length):  
	arr = np.zeros(1000)
    for j in xrange(0,1000):
        sim_avg_annots1000[i][0] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 4000: j*geneNum + geneNum])
		arr[j] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 4000: j*geneNum + geneNum])
    sim_avg_annots1000[i][0] = (1.0/1000.0)*sim_avg_annots1000[i][0] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots1000[i][0])*(arr[j] - sim_avg_annots1000[i][0])
	sim_var_annots1000[i][0] = (1.0/1000.0)*np.sum(arr)	
    print i    
for i in xrange (0, gM_direct_length): 
	arr = np.zeros(1000)
    for j in xrange(1000,2000):
        sim_avg_annots1000[i][1] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 4000: j*geneNum + geneNum])
		arr[j - 1000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 4000: j*geneNum + geneNum])
    sim_avg_annots1000[i][1] = (1.0/1000.0)*sim_avg_annots1000[i][1] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots1000[i][1])*(arr[j] - sim_avg_annots1000[i][1])
	sim_var_annots1000[i][1] = (1.0/1000.0)*np.sum(arr)	
    print i
for i in xrange (0, gM_direct_length): 
	arr = np.zeros(1000)
    for j in xrange(2000,3000):
        sim_avg_annots1000[i][2] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 4000: j*geneNum + geneNum])
		arr[j - 2000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 4000: j*geneNum + geneNum])
    sim_avg_annots1000[i][2] = (1.0/1000.0)*sim_avg_annots1000[i][2] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots1000[i][2])*(arr[j] - sim_avg_annots1000[i][2])
	sim_var_annots1000[i][2] = (1.0/1000.0)*np.sum(arr)	
    print i
for i in xrange (0, gM_direct_length): 
	arr = np.zeros(1000)
    for j in xrange(3000,4000):
        sim_avg_annots1000[i][3] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 4000: j*geneNum + geneNum])
		arr[j - 3000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 4000: j*geneNum + geneNum])
    sim_avg_annots1000[i][3] = (1.0/1000.0)*sim_avg_annots1000[i][3] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots1000[i][3])*(arr[j] - sim_avg_annots1000[i][3])
	sim_var_annots1000[i][3] = (1.0/1000.0)*np.sum(arr)	
    print i   
for i in xrange (0, gM_direct_length):   
	arr = np.zeros(1000)
    for j in xrange(4000,geneNum):
        sim_avg_annots1000[i][4] += (1.0/1000.0)*np.sum(sim[i][j*geneNum + 4000: j*geneNum + geneNum])
		arr[j - 4000] = (1.0/1000.0)*np.sum(sim[i][j*geneNum + 4000: j*geneNum + geneNum])
    sim_avg_annots1000[i][4] = (1.0/1000.0)*sim_avg_annots1000[i][4] 
	for j in xrange(0,1000):
		arr[j] = (arr[j] - sim_avg_annots1000[i][4])*(arr[j] - sim_avg_annots1000[i][4])
	sim_var_annots1000[i][4] = (1.0/1000.0)*np.sum(arr)	
    print i
    
#Figure_________________________________________________________________________
f = py.figure(figsize=(13,20), dpi=80)

f.suptitle('Direct Groupwise Similarities')

subplots = np.array([f.add_subplot(6,3,1), f.add_subplot(6,3,2), f.add_subplot(6,3,3), f.add_subplot(6,3,4), f.add_subplot(6,3,5), f.add_subplot(6,3,6),
                     f.add_subplot(6,3,7), f.add_subplot(6,3,8), f.add_subplot(6,3,9), f.add_subplot(6,3,10), f.add_subplot(6,3,11), f.add_subplot(6,3,12),
                     f.add_subplot(6,3,13), f.add_subplot(6,3,14), f.add_subplot(6,3,15), f.add_subplot(6,3,16), f.add_subplot(6,3,17), f.add_subplot(6,3,18)])
subplots[0].set_title('SIM_FRAMEWORK_DAG_SET_BADER_2003', y=1.08, fontsize=7)
subplots[1].set_title('SIM_FRAMEWORK_DAG_SET_BRAUN_BLANQUET_1932', y=1.08, fontsize=7)
subplots[2].set_title('SIM_FRAMEWORK_DAG_SET_DICE_1945', y=1.08, fontsize=7)
subplots[3].set_title('SIM_FRAMEWORK_DAG_SET_JACCARD_1901', y=1.08, fontsize=7)
subplots[4].set_title('SIM_FRAMEWORK_DAG_SET_KNAPPE_2004', y=1.08, fontsize=7)
subplots[5].set_title('SIM_FRAMEWORK_DAG_SET_KORBEL_2002', y=1.08, fontsize=7)
subplots[6].set_title('SIM_FRAMEWORK_DAG_SET_MARYLAND_BRIDGE_2003', y=1.08, fontsize=7)
subplots[7].set_title('SIM_FRAMEWORK_DAG_SET_OCHIAI_1957', y=1.08, fontsize=7)
subplots[8].set_title('SIM_FRAMEWORK_DAG_SET_SIMPSON_1960', y=1.08, fontsize=7)
subplots[9].set_title('SIM_FRAMEWORK_DAG_SET_SOKAL_SNEATH_1963', y=1.08, fontsize=7)
subplots[10].set_title('SIM_FRAMEWORK_DAG_SET_TVERSKY_1977', y=1.08, fontsize=7)   
subplots[11].set_title('SIM_GROUPWISE_DAG_ALI_DEANE', y=1.08, fontsize=7)
subplots[12].set_title('SIM_GROUPWISE_DAG_NTO', y=1.08, fontsize=7)
subplots[13].set_title('SIM_GROUPWISE_DAG_NTO_MAX', y=1.08, fontsize=7)
subplots[14].set_title('SIM_GROUPWISE_DAG_UI', y=1.08, fontsize=7) 
subplots[15].set_title('SIM_GROUPWISE_DAG_GIC + ICI_SANCHEZ_2011 ', y=1.08, fontsize=7)
subplots[16].set_title('SIM_GROUPWISE_DAG_GIC + ICI_RESNIK_1995 ', y=1.08, fontsize=7)
subplots[17].set_title('SIM_GROUPWISE_DAG_GIC + ICI_PROB_OCCURENCE_PROPAGATED', y=1.08, fontsize=7) 

for i in xrange(0,gM_direct_length):
    subplots[i].set_xlabel('Annotation Size', labelpad = 7, fontsize = 7)
    subplots[i].set_ylabel('Similarity', labelpad = 7, fontsize = 7)
    #subplots[i].xticks(np.array([1000,2000,3000,4000,5000]))
    #subplots[i].xlim(0,5001)
    subplots[i].set_ylim(0,1.0)
    subplots[i].set_yticklabels([0.0,'','','','',1.0], fontsize=6)
    #subplots[i].set_xticklabels(['','GP1000','GP2000','GP3000','GP4000','GP5000'], fontsize=6)
'''  
ax.set_xticks(ind+width)
xtickNames = ax.set_xticklabels( ['Group'+str(i) for i in range(1,6)])
plt.setp(xtickNames, rotation=45, fontsize=10)
'''

width = 0.15
#for each direct groupwise measures:
for i in xrange (0, gM_direct_length): 
    #subplots[i].bar(annots, mean_annots[i][:], width, color='black',yerr=var_annots[i][:], error_kw=dict(elinewidth=2,ecolor='red'))
    subplots[i].bar(annots, sim_avg_annots1[i][:], width, color='black', yerr = sim_var_annots1[i][:], label = "Average sim to GPs with annotation size 1")
    subplots[i].bar(annots + width, sim_avg_annots10[i][:], width, color='blue', yerr = sim_var_annots10[i][:], label = "Average sim to GPs with annotation size 10")
    subplots[i].bar(annots + 2.0*width, sim_avg_annots50[i][:], width, color='red', yerr = sim_var_annots50[i][:], label = "Average sim to GPs with annotation size 50")
    subplots[i].bar(annots + 3.0*width, sim_avg_annots100[i][:], width, color='green', yerr = sim_var_annots100[i][:], label = "Average sim to GPs with annotation size 100")
    subplots[i].bar(annots + 4.0*width, sim_avg_annots1000[i][:], width, color='yellow', yerr = sim_var_annots1000[i][:], label = "Average sim to GPs with annotation size 1000")
    subplots[i].set_xlim(-width,5.0 + width)
    subplots[i].set_xticks(annots+0.35)
    subplots[i].set_xticklabels(['1','10','50','100','1000'], fontsize=6)
    for tic in subplots[i].xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
    #subplots[i].set_xticklabels([str(i) for i in range(1,6)])
    #subplots[i].plot(geneProducts, mean[i][:], '-', linewidth = 0.5, label = "mean similarity")
    #subplots[i].plot(geneProducts, var[i][:], '-', linewidth = 0.5, label = "variance in similarity") 
    if i==0:
       subplots[i].legend(loc = 2, fontsize=6)
        
f.subplots_adjust(wspace=0.5, hspace=0.5)


py.savefig("C:\Users\#ballin\Desktop\individualMeasures_data/directGroupwiseMeasures/meanvargraphs.png")

print sim_avg_annots1

py.show()



    
