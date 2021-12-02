#!/usr/bin/python2.7.3
# -*- coding: latin-1 -*-

#==============================================================================
#title           :K-S_test.py
#description     :It is the kolmogorov-Smirnov test on the zeta_vir and zeta_tild
#author          :Dyaa Chbib
#date            :2016_06_01
#version         :0.1
#python_version  :2.7.3
#==============================================================================

import time
print time.ctime()

from pylab import *
import pylab as pl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from random import randint
from numpy import *
import scipy.integrate as si
import math
from random import gauss
import random
import os
import sympy
from scipy.optimize import curve_fit
import scipy.optimize as optimization
import scipy.stats as chi
from scipy.optimize import leastsq
import asciitable as asc
from scipy.integrate import quad
from scipy.integrate import dblquad
from scipy.integrate import tplquad
from scipy import stats
from scipy.stats import norm
#from decimal import *

import Methods
from Methods import Functions
from Cosmological_Model import Model
from HubbleDiagram import ApparentMagnitude
from GradientMethod import Gradient
from Savedata import Save
from Plot import Plotting

if socket.gethostname()=='pc-dchbib' or socket.gethostname()=='cpt31':
	try:
		current_path1 = os.getcwd()
		os.chdir("/media/Elements/")
		current_path = os.getcwd()
		if current_path == '/media/Elements':
			Name_of_machine = '/media/Elements/'
			os.chdir(current_path1)
		else: 
		     print "except OSError ??"
	except OSError: 
		Name_of_machine = '/media/Elements__/'
elif socket.gethostname()=='port-yassine':
	Name_of_machine = '/media/manal/Elements/'
	
elif socket.gethostname()=='mchbib-HP-Pavilion-Notebook' or socket.gethostname()=='mchbib':
	Name_of_machine = '/media/mchbib/Elements/'	

else:
	print "You should determine the host name of machine", "\n"

#path_of_sample = Name_of_machine+'Zip_file/Python/Simulations/POO/OO/NullCorrelation/Mixte_Plot/Image_Simulation/OneHundredSimulation/Second_OneHundredSimulation_2000_Objects/SimulationData_0.71-0.3_NbOfObjects740_limitMag26_NbOfSimulation0Beta2.05.dat'
#path_of_sample = Name_of_machine+'Zip_file/Python/Simulations/POO/OO/NullCorrelation/Mixte_Plot/Image_Simulation/OneHundredSimulation/Second_OneHundredSimulation_2000_Objects/SimulationData_0.71-0.3_NbOfObjects740_limitMag26_NbOfSimulation0Beta2.05.dat'
#path_of_sample = Name_of_machine+'Zip_file/Python/Simulations/POO/OO/NullCorrelation/Mixte_Plot/Image_Simulation/OneHundredSimulation/Second_OneHundredSimulation_2000_Objects/SimulationData_0.71-0.3_NbOfObjects740_limitMag26_NbOfSimulation0Beta2.05.dat'
#path_of_sample = Name_of_machine+'Zip_file/Python/Simulations/POO/OO/NullCorrelation/Mixte_Plot/Image_Simulation/OneHundredSimulation/Second_OneHundredSimulation_2000_Objects/SimulationData_0.71-0.3_NbOfObjects740_limitMag26_NbOfSimulation0Beta2.07.dat'
#path_of_sample = Name_of_machine+'Zip_file/Python/Simulations/POO/OO/SNIa/NullCorrelation/NullCorrelation_Simulation/NullCorrelation_Simulation_lifetime_code/NullCorrelation_Simulation_Like_QSOs_Code/Simulation_Data/SimulatedData_0.701-0.3_alfa1.2_bbeta2.1_s01_c00_sigmas0.12_sigmac0.09_M0-19.8_NbOfObjects500_limitMag26_NbOfSimulation0_Beta2.05.dat'
#path_of_sample = Name_of_machine+'Zip_file/Python/Simulations/POO/OO/SNIa/NullCorrelation/NullCorrelation_Simulation/NullCorrelation_Simulation_lifetime_code/NullCorrelation_Simulation_Like_QSOs_Code/Simulation_Data/SimulatedData_0.9-0.2_alfa0.14_bbeta3.1_s01_c00_sigmas0.12_sigmac0.09_M0-19.01_NbOfObjects500_limitMag26_NbOfSimulation0_Beta2.21.dat'
#path_of_sample = Name_of_machine+'Zip_file/Python/Simulations/POO/OO/SNIa/NullCorrelation/NullCorrelation_Simulation/NullCorrelation_Simulation_lifetime_code/NullCorrelation_Simulation_Like_QSOs_Code/Simulation_Data/SimulatedData_0.85-0.15_alfa0.14_bbeta3.1_s01_c00_sigmas0.12_sigmac0.09_M0-19.01_NbOfObjects500_limitMag26_NbOfSimulation0_Beta2.19.dat'
#path_of_sample = Name_of_machine+'Zip_file/Python/Simulations/POO/OO/SNIa/NullCorrelation/NullCorrelation_Simulation/NullCorrelation_Simulation_lifetime_code/NullCorrelation_Simulation_Like_QSOs_Code/Simulation_Data/SimulatedData_0.85-0.15_alfa0.14_bbeta3.1_s01_c00_sigmas0.12_sigmac0.09_M0-19.01_NbOfObjects500_limitMag26_NbOfSimulation0_Beta2.19.dat'
#path_of_sample = Name_of_machine+'Zip_file/Python/Simulations/POO/OO/SNIa/NullCorrelation/NullCorrelation_Simulation/NullCorrelation_Simulation_lifetime_code/NullCorrelation_Simulation_Like_QSOs_Code/Simulation_Data/SimulatedData_0.701-0.3_alfa0.14_bbeta3.1_s01_c00_sigmas0.08_sigmac0.06_M0-19.01_NbOfObjects500_limitMag26_NbOfSimulation0_Beta2.15.dat'
path_of_sample = Name_of_machine+'Zip_file/Python/Simulations/POO/OO/SNIa/NullCorrelation/NullCorrelation_Simulation/NullCorrelation_Simulation_lifetime_code/NullCorrelation_Simulation_Like_QSOs_Code/Simulation_Data/SimulatedData_0.9-0.2_alfa0.14_bbeta3.1_s01_c00_sigmas0.07_sigmac0.06_M0-19.01_NbOfObjects500_limitMag26_NbOfSimulation0_Beta2.22.dat'

print "path_of_sample: t = asc.read('"+str(path_of_sample)+"', guess=False)"

t = asc.read(path_of_sample, guess=False)
slist, clist, Mlist, m_app, redshift, redshift_l, Volume, Volume_l, Volumebias, Mbias, m_bias, redshiftbias, lifetime, lookbacktime, redshift_1, redshift_2, lookbacktime_1, lookbacktime_2, Weighting, factorWeighting, redred, z1, m_thlistplt, degree_of_Polyn, Ncoefficients_of_KcorrPolyn = t['slist'], t['clist'], t['Mlist'], t['m_app'], t['redshift'], t['redshift_l'], t['Volume'], t['Volume_l'], t['Volumebias'], t['Mbias'], t['m_bias'], t['redshiftbias'], t['lifetime'], t['lookbacktime'], t['redshift_1'], t['redshift_2'], t['lookbacktime_1'], t['lookbacktime_2'], t['Weighting'], t['factorWeighting'], t['redred'], t['z1'], t['m_thlistplt'], t['degree_of_Polyn'], t['Ncoefficients_of_KcorrPolyn']

alfa, bbeta = 0.14, 3.1
alfa, bbeta = 3.14, 5.1

MM = Mlist - alfa*(slist - 1) + bbeta*clist

zeta_vir = m_app - MM


omegalist = linspace(np.float64(10**-6), np.float64(1.0), 30) 
lamdalist = linspace(np.float64(10**-6), np.float64(2.0), 30) 
alfalist = linspace(np.float64(10**-6), np.float64(0.5), 30) 
bbetalist = linspace(np.float64(2.5), np.float64(3.5), 30) 
M_0list = linspace(np.float64(-21), np.float64(-18), 30) 


All_models = []
for i in range(len(lamdalist)):                    
    for j in range(len(omegalist)):
        All_models.append((lamdalist[i], omegalist[j]))

All_models = np.asarray(All_models)
X = All_models[:,1]
Y = All_models[:,0]


DD = []

for i in range(len(X)):                           
    Bounced_limit = Functions(Y[i], X[i]).NoBigBang(X[i], 'permission')
    if Bounced_limit==0:
        DD.append((Y[i], X[i], nan))
    else:
        D = []     
        #zeta_tild = ApparentMagnitude(redshift_SNLS, Y[i], X[i], 0, 'Realdata_NO', 0, 0).m_theor(redshift_SNLS, Y[i], X[i], 0,0,0,0,0)
        zeta_tild = ApparentMagnitude(redshift, Y[i], X[i], 0, 'Realdata_NO', 0, 0).m_theor(redshift, Y[i], X[i], 0,0,0,0,0)
        for k in range(len(zeta_tild)):
            x = zeta_tild[k]
            Heavisidelist = Functions(0.7,0.3).Heavisidefuction(x - zeta_tild)
            F_x = sum(Heavisidelist)/float(len(zeta_tild))
            Heavisidelist = Functions(0.7,0.3).Heavisidefuction(x - zeta_vir)
            F_H_x = sum(Heavisidelist)/float(len(zeta_tild))
            D.append(abs(F_H_x - F_x))
        DD.append((Y[i], X[i], max(D)))
 
DDD = array(DD)[:,2]
DDD = list(array(DD)[:,2])
print "OmegaOfCorrNull[DDD.index(min(DDD))] = ", X[DDD.index(min(DDD))], "\n"
print "LambdaOfCorrNull[DDD.index(min(DDD))] = ", Y[DDD.index(min(DDD))], "\n"

print "D_max of zeta_vir and zeta_tild distributions = ", min(DDD), "\n"
print "For 1 %  level = ", sqrt(-0.5*np.log(0.01/2.)/len(redshift)), "\n"
print "For 2 %  level = ", sqrt(-0.5*np.log(0.02/2.)/len(redshift)), "\n"
print "For 5 %  level = ", sqrt(-0.5*np.log(0.05/2.)/len(redshift)), "\n"
print "For 8 %  level = ", sqrt(-0.5*np.log(0.08/2.)/len(redshift)), "\n"
print "For 10 %  level = ", sqrt(-0.5*np.log(0.10/2.)/len(redshift)), "\n"
print "For 15 %  level = ", sqrt(-0.5*np.log(0.15/2.)/len(redshift)), "\n"
print "For 20 %  level = ", sqrt(-0.5*np.log(0.20/2.)/len(redshift)), "\n"
print "For 25 %  level = ", sqrt(-0.5*np.log(0.25/2.)/len(redshift)), "\n"
print "For 30 %  level = ", sqrt(-0.5*np.log(0.30/2.)/len(redshift)), "\n"
print "For 55 %  level = ", sqrt(-0.5*np.log(0.55/2.)/len(redshift)), "\n"
print "For 85 %  level = ", sqrt(-0.5*np.log(0.85/2.)/len(redshift)), "\n"
print "For 95 %  level = ", sqrt(-0.5*np.log(0.95/2.)/len(redshift)), "\n"

# Plot 3D

DDD = array(DD)[:,2]
DDD_r = DDD*1.0
DDD_r = DDD_r.reshape(sqrt(len(DDD_r)), sqrt(len(DDD_r)))
X_r = X*1.0
X_r = X_r.reshape(sqrt(len(X_r)), sqrt(len(X_r)))
Y_r = Y*1.0
Y_r = Y_r.reshape(sqrt(len(Y_r)), sqrt(len(Y_r)))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X_r, Y_r, DDD_r, rstride=4, cstride=4, color='g')


fig = plt.figure()
ax = fig.add_subplot(111)
#col = plt.scatter(All_models[:,1], All_models[:,0], marker='.', s=150, c=np.exp(-DDD), linewidths=0.3, cmap=plt.cm.Spectral_r) 	# Pour ploter le maximum
col = plt.scatter(All_models[:,1], All_models[:,0], marker='.', s=150, c=DDD, linewidths=0.3, cmap=plt.cm.Spectral_r)
cax = fig.colorbar(col, orientation='vertical',format='%.6f')
plt.xlabel('$\\Omega_{0}$', fontsize=14)
plt.ylabel('$\\lambda_{0}$', fontsize=14)
Lamdalimit=[]
for i in range(len(omegalist)):                                  
	Lamdalimit.append(Functions(0.7, 0.3).NoBigBang(omegalist[i], 'limit'))

plt.fill_between(omegalist, Lamdalimit, max(Lamdalimit), color=(0.5,0.5,0.5))
plt.xlim(0.01, max(omegalist))                                                                                                                                                                                                                                                                                                                                                                                                               
plt.ylim(0.01, max(lamdalist))
#plt.plot([0, 0.3],[0.71, 0.71], '-', color='r')
#plt.plot([0.3, 0.3],[0, 0.71], '-', color='r')
#plt.plot([0.3],[0.71], marker='o', color='r')

#============**** To make the contours of confidence levels ****=============  
All_modelslambda = All_models[:,0].reshape(sqrt(len(All_models[:,0])), sqrt(len(All_models[:,0]))) # For matrix of lambda
All_modelsomega = All_models[:,1].reshape(sqrt(len(All_models[:,1])), sqrt(len(All_models[:,0]))) # For matrix of omega

CS = pl.contour(All_modelsomega, All_modelslambda, DDD_r, 3, linewidths=np.arange(.5, 4, .5), colors=('r', 'green', 'blue', (1,1,0), '#afeeee', '0.5'))
plt.clabel(CS, inline=1, fontsize=10)
#============***************************************************=============  


plt.show()
