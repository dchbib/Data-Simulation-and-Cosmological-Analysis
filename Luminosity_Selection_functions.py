#!/usr/bin/python2.7.3
# -*- coding: latin-1 -*-

#==============================================================================
#title           :main.py
#description     :Bulding of the luminosity function  
#author          :Dyaa Chbib
#date            :2014_11_28
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
from scipy.integrate import *
from scipy import stats
from scipy.stats import norm
#from decimal import *
from astropy.table import Table, Column
import socket

import Methods
from Methods import Functions
from Cosmological_Model import Model
from HubbleDiagram import ApparentMagnitude
from GradientMethod2 import Gradient
from GradientMethod3 import Gradient
from Savedata import Save
from Plot import Plotting
from Mock_Sample import SIMULATED
from Mock_Sample_SN_2D_CDF import SIMULATED


Betalist = linspace(np.float64(0.01), np.float64(5), 500)
Beta = array([array([a]) for a in Betalist])

#==============K-corrections=======================================
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

tt = asc.read(Name_of_machine+'Zip_file/Python/Simulations/POO/OO/NullCorrelation/NullCorrelation_Realdata/SDSS-DR3/K-correction.dat', guess=False)
		
degree_of_Polyn, Ncoefficients_of_KcorrPolyn = tt['redsh'], tt['kcorre']

dkcorr = [] 

for i in range(len(tt['redsh'])):
	if i==0 or i==(len(tt['redsh'])-1):
		dkcorr.append(tt['kcorre'][i])
	else:
		dkcor = (tt['kcorre'][i+1] - tt['kcorre'][i-1])/(tt['redsh'][i+1] - tt['redsh'][i-1])
		dkcorr.append(dkcor)
				
dkcorr = array(dkcorr)
#====================================================================

lambda0, omega0, M_0, sigmaM, m_lim, zform, z_max, z_min, Sizeofsample, NumbOfSamples = 0.7, 0.3, -19.01, 0.3, 25.04, 16, 16, 10**-15.95, 300, 3
lambda0, omega0 = 1.2, 0.1
lambda0, omega0 = 1.275, 0.132
lambda0, omega0 = 0.7, 0.3
#lambda0, omega0 = 0.8, 0.2
#lambda0, omega0 = 1.3, 0.1
#lambda0, omega0 = 1.03, 0.4
#lambda0, omega0 = 0.85, 0.25
#lambda0, omega0 = 1.636, 0.2571
alfa, bbeta = 1.141, 3.14
alfa, bbeta = 0.141, 3.14
Realdata = 'Realdata_NO'
degree_of_Polyn, Ncoefficients_of_KcorrPolyn = 0, 0
M0, s0, c0 = -19.00, 0, 0

#Direct
alfa, bbeta = 0.141, 3.14
#alfa, bbeta = 1.2, 2.1
sigmaM, sigmas, sigmac, sigma_eps = 0.25, 0.99, 0.085, 0.09

#Indirect
#alfa, bbeta = 5, 3.14
#sigmaM, sigmas, sigmac, sigma_eps = 0.25, 0.99, 0.085, 0.6

degree_of_Polyn, Ncoefficients_of_KcorrPolyn = [], []

paths_of_samples = SIMULATED(lambda0, omega0).Sample(alfa, bbeta, s0, c0, sigmas, sigmac, M0, sigmaM, sigma_eps, m_lim, zform, z_max, z_min, Sizeofsample, NumbOfSamples, 'Realdata_NO', 'Direct', degree_of_Polyn, Ncoefficients_of_KcorrPolyn)

#plt.show()

for bb in range(0, NumbOfSamples):
	
	Number_of_simulation = bb
	
	path_of_sample = paths_of_samples[bb]
	print "path_of_sample: t = asc.read('"+str(path_of_sample)+"', guess=False)"
	
	t = asc.read(path_of_sample, guess=False)
	"""
	m_app, redshift, = t['m_app'], t['redshift']
	
	try:
		beta = float(path_of_sample[len(path_of_sample)-8:len(path_of_sample)-4])
	
	except ValueError:
		beta = float(path_of_sample[len(path_of_sample)-7:len(path_of_sample)-4])
	"""
	
	slist, clist, Mlist, m_app, redshift, redshift_l, Volume, Volume_l, Volumebias, Mbias, m_bias, redshiftbias, lifetime, lookbacktime, redshift_1, redshift_2, lookbacktime_1, lookbacktime_2, factorWeighting_s, factorWeighting_c, redred, z1, m_thlistplt, degree_of_Polyn, Ncoefficients_of_KcorrPolyn = t['slist'], t['clist'], t['Mlist'], t['m_app'], t['redshift'], t['redshift_l'], t['Volume'], t['Volume_l'], t['Volumebias'], t['Mbias'], t['m_bias'], t['redshiftbias'], t['lifetime'], t['lookbacktime'], t['redshift_1'], t['redshift_2'], t['lookbacktime_1'], t['lookbacktime_2'], t['factorWeighting_s'], t['factorWeighting_c'], t['redred'], t['z1'], t['m_thlistplt'], t['degree_of_Polyn'], t['Ncoefficients_of_KcorrPolyn']
	factorWeighting = t['factorWeighting']
	m_SNLS = m_app*1.0
	m_app = m_app - bbeta*clist
        Mlist = Mlist + bbeta*clist
        

        print "len(m_app)", len(m_app), "\n"
        print max(m_app)
        print max(redshift)


	#===============================SUITE_of_K-corrections===========================================
	degree_of_Polyn, Ncoefficients_of_KcorrPolyn = 0, 0
	#=================================================================================================
		
	#============determine of beta (constant used in the weighting factor term)=======================
        Betalist = linspace(np.float64(0.01), np.float64(3.9), 350)
        Beta = array([array([a]) for a in Betalist])
        WeightingManifold = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, Beta)
        #ten_power_m_M = (10**((m_app - bbeta*clist - Mlist)/5.))**Beta
        #Jacob = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Jacobian(redshift)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
        #WeightingManifold = Jacob * ten_power_m_M
        WeightingSum = [a/sum(a) for a in WeightingManifold]
        #DetermineBeta = [abs(max(b) - min(b)) for b in WeightingSum]
        DetermineBeta = [abs((max(b) - min(b))/(max(b) + min(b))) for b in WeightingSum]
        beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
        Weighting = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, beta)
        #ten_power_m_M_beta = (10**((m_app - bbeta*clist - Mlist)/5.))**beta
        #Weighting = (Jacob*ten_power_m_M_beta)
        SumWeighting = sum(Weighting)
        factorWeighting = array(Weighting)/SumWeighting
        #===================================================================================================
	print "beta = ", beta, "\n"
	m_app_c = m_app - Ncoefficients_of_KcorrPolyn
	#Mlist = m_app_c - ApparentMagnitude(redshift, lambda0, omega0, 0.0, 'Realdata_NO', degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(redshift, lambda0, omega0, 0.0)
	Mlist = Mlist*1.0
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Mlist = m_app - ApparentMagnitude(redshift, lambda0, omega0, 0.0, 'Realdata_OK', degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(redshift, lambda0, omega0, 0.0)
	#Volume = Model(lambda0, omega0).curvature(redshift, lambda0, omega0, 0.0)[1]
	#plt.figure()
	#for i in range(len(redshift)):               
	#	plt.plot(Mlist[i], Volume[i], marker='+', color='r')
	#plt.xlabel('$Absolute$'+' '+'$magnitude$'+' '+'$M$', fontsize=16)
	#plt.ylabel('$Volume$', fontsize=16)
	#plt.show()
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	plt.figure(bb)
	
	Mlist1, redshift1, factorWeighting1, m_app1,m_app_c1, m_SNLS1 = list(Mlist*1.0), list(redshift*1.0), list(factorWeighting*1.0), list(m_app*1.0), list(m_app_c*1.0), list(m_SNLS)
	Mlist11, m_app11, m_app_c11, factorWeighting11, redshift11, m_SNLS11 = [], [], [], [], [], []
	for i in range(len(m_app)):
	    m_app11.append(min(m_app1))
	    m_SNLS11.append(m_SNLS1[m_app1.index(min(m_app1))])
	    m_SNLS1.remove(m_SNLS1[m_app1.index(min(m_app1))])
	    m_app1.remove(min(m_app1))
	    m_app_c11.append(min(m_app_c1))
	    m_app_c1.remove(min(m_app_c1))
	
	
	m_app11 = array(m_app11) 
	m_SNLS11 = array(m_SNLS11)
	
	Fst =[]
	for i in range(len(m_app11)):
	    Heavisidelist = Functions(lambda0, omega0).Heavisidefuction(m_app11[i] - m_app)
	    F_st2 = sum(factorWeighting*Heavisidelist)
	    Fst.append(F_st2)

	Fst_of_phistar = Fst
	plt.subplot(234)
	plt.step(m_app11, Fst, 'g', label="Empirical CDF")
	plt.legend(loc='upper left',prop={'size':12})
	plt.xlabel('$Apparent$'+' '+'$magnitude$'+' '+'$m$', fontsize=13)
	plt.ylabel('$CDF$'+' '+'$of$'+' '+'$\\varphi_{w}(m)$', fontsize=13)
	plt.xlim(14, max(m_app11)+1)
	plt.ylim(-0.0001, 1.1)
	
	f_gauss, m_app3, m_app_c3, m_SNLS3 = [], [], [], []
	for i in range(int(sqrt(Sizeofsample)/2), len(m_app11)-int(sqrt(Sizeofsample)/2)):
	    h1, h2 = m_app11[i+int(sqrt(Sizeofsample)/2)] - m_app11[i], m_app11[i] - m_app11[i-int(sqrt(Sizeofsample)/2)]
	    ghg=1
	    while h1==0 or h2==0:
		dfg = i+int(sqrt(Sizeofsample)/2)+ghg
		if dfg < len(m_app11):
			h1, h2 = m_app11[i+int(sqrt(Sizeofsample)/2)+ghg] - m_app11[i], m_app11[i] - m_app11[i-int(sqrt(Sizeofsample)/2)-ghg]
			ghg = ghg + 20
		else:
			break
		      
	    m1 = m_app11[i] + h1
	    m2 = m_app11[i] - h2
	    Heavisidelist = Functions(lambda0, omega0).Heavisidefuction(m1 - m_app)
	    F_st1 = sum(factorWeighting*Heavisidelist)
	    Heavisidelist = Functions(lambda0, omega0).Heavisidefuction(m2 - m_app)
	    F_st2 = sum(factorWeighting*Heavisidelist)
	    fm = (F_st1 - F_st2)/(h1+h2)
	    f_gauss.append(fm)
	    m_app3.append(m_app11[i])
	    m_app_c3.append(m_app_c11[i])
	    m_SNLS3.append(m_SNLS11[i])
	
	phistar = f_gauss
	f_gauss = array(f_gauss)
	m_app3 = array(m_app3)
	m_app_c3 = array(m_app_c3)
	m_SNLS3 = array(m_SNLS3)
	
	plt.subplot(235)
	
	for i in range(len(m_app3)):
		plt.plot(m_app3[i], f_gauss[i]/mean(f_gauss), marker='+', color='g')

	plt.plot(m_app3, f_gauss/mean(f_gauss), 'g', label="$\\varphi_{w}(m)$")
	plt.legend(loc='upper left',prop={'size':12})
	plt.xlabel('$Apparent$'+' '+'$magnitude$'+' '+'$m$', fontsize=13)
	plt.ylabel('$\\varphi_{w}(m)$', fontsize=13)
	plt.xlim(16.5, max(m_app11)+1)
	plt.ylim(-0.2,2.5)
	
	plt.subplot(236)
	plt.plot(linspace(0.1, max(m_app11)+10, 10000), Functions(lambda0, omega0).Heavisidefuction(m_lim - linspace(0.1, max(m_app11)+10, 10000)), 'r')
	plt.text(m_lim, -0.16, '$m_{lim}$', color='r')
	#f_gauss = f_gauss * 10**(-beta*m_app3/5.)
	f_gauss = f_gauss * 10**(-beta*m_SNLS3/5.)
	phistar_time_TenPower = f_gauss
	
	for i in range(len(m_app3)):
		plt.plot(m_app3[i], f_gauss[i]/mean(f_gauss), marker='+', color='m')

	plt.plot(m_app3, f_gauss/mean(f_gauss), 'm', label="$\\varphi(m) \\propto \\varphi_{w}(m).10^\\frac{-\\beta.m}{5}$")
	plt.legend(loc='upper right',prop={'size':12})
	plt.xlabel('$Apparent$'+' '+'$magnitude$'+' '+'$m$', fontsize=13)
	plt.ylabel('$\\varphi(m)$', fontsize=13)
	plt.xlim(16.5, max(m_app11)+1)
	plt.ylim(-0.2,2.5)
	
	
	plt.subplot(231)   
	
	Mlist1, redshift1, factorWeighting1, m_app1 = list(Mlist*1.0), list(redshift*1.0), list(factorWeighting*1.0), list(m_app*1.0)
	Mlist11, factorWeighting11, redshift11 = [], [], []

	for i in range(len(Mlist)):                             
	    Mlist11.append(min(Mlist1))
	    Mlist1.remove(min(Mlist1))
 
	Fst =[]
	for i in range(len(Mlist11)):
	    Heavisidelist = Functions(lambda0, omega0).Heavisidefuction(Mlist11[i] - Mlist)
	    F_st2 = sum(factorWeighting*Heavisidelist)
	    Fst.append(F_st2)
	    
	Fst_of_fMstar = Fst
	plt.step(Mlist11, Fst, 'g', label="Empirical CDF")
	plt.legend(loc='upper left',prop={'size':12})
	plt.xlabel('$Absolute$'+' '+'$magnitude$'+' '+'$M$', fontsize=13)
	plt.ylabel('$CDF$'+' '+'$of$'+' '+'$f_{w}(M)$', fontsize=13)
	
	plt.subplot(232)

	f_gauss, Mlist3 = [], []
	for i in range(int(sqrt(Sizeofsample)/2), len(Mlist11)-int(sqrt(Sizeofsample)/2)):
	    h1, h2 = Mlist11[i+int(sqrt(Sizeofsample)/2)] - Mlist11[i], Mlist11[i] - Mlist11[i-int(sqrt(Sizeofsample)/2)]
	    M1 = Mlist11[i] + h1
	    M2 = Mlist11[i] - h2
	    Heavisidelist = Functions(lambda0, omega0).Heavisidefuction(M1 - Mlist)
	    F_st1 = sum(factorWeighting*Heavisidelist)
	    Heavisidelist = Functions(lambda0, omega0).Heavisidefuction(M2 - Mlist)
	    F_st2 = sum(factorWeighting*Heavisidelist)
	    fM = (F_st1 - F_st2)/(h1+h2)
	    f_gauss.append(fM)
	    Mlist3.append(Mlist11[i])

	for i in range(len(Mlist3)):
		plt.plot(Mlist3[i], f_gauss[i]/max(f_gauss), marker='+', color='g')
		
	plt.plot(Mlist3, f_gauss/max(f_gauss), 'g', label="$f_{w}(M)$")
	plt.legend(loc='upper left',prop={'size':12})
	plt.xlabel('$Absolute$'+' '+'$magnitude$'+' '+'$M$', fontsize=13)
	plt.ylabel('$Weighted$'+' '+'$luminosity$'+' '+'$function$'+' '+'$f_{w}(M)$', fontsize=13)

	plt.subplot(233)
	fMstar = f_gauss
	
	f_gauss = array(f_gauss)
	Mlist3 = array(Mlist3)
	f_gauss = f_gauss * 10**(beta*Mlist3/5.)
		
	fMstar_time_TenPower = f_gauss
	
	for i in range(len(Mlist3)):
		plt.plot(Mlist3[i], f_gauss[i]/max(f_gauss), marker='+', color='m')
		
	plt.plot(Mlist3, f_gauss/max(f_gauss), 'm', label="$f(M) \\propto f_{w}(M).10^\\frac{\\beta.M}{5}$")
	plt.legend(loc='upper left',prop={'size':12})
	plt.xlabel('$Absolute$'+' '+'$magnitude$'+' '+'$M$', fontsize=13)
	plt.ylabel('$Luminosity$'+' '+'$function$'+' '+'$f(M)$', fontsize=13)
	
	filename = 'fM_phim_'+str(lambda0)+'-'+str(omega0)+'_NbOfObjects'+str(Sizeofsample)+'_limitMag'+str(m_lim)+'Beta'+str(beta)+'.dat'

	listofarrays = [Mlist11, m_app11, m_app_c11, factorWeighting, Fst_of_fMstar, Fst_of_phistar, Mlist3, m_app3, m_app_c3, fMstar, fMstar_time_TenPower, phistar, phistar_time_TenPower]
	stringlistofarrays = ['Mlist11', 'm_app11', 'm_app_c11', 'factorWeighting', 'Fst_of_fMstar', 'Fst_of_phistar', 'Mlist3', 'm_app3', 'm_app_c3', 'fMstar', 'fMstar_time_TenPower', 'phistar', 'phistar_time_TenPower']
	

	path_tosave = Save(listofarrays, stringlistofarrays, filename).ReadWrite(True, 'work')

	#subprocess.call(['cp', path_tosave, Name_of_machine+'Zip_file/Python/Simulations/POO/OO/SNIa/NullCorrelation/NullCorrelation_Realdata/JLA/SN_SNLS/'+filename], shell = False)

	print path_tosave
	print Name_of_machine+'Zip_file/Python/Simulations/POO/OO/SNIa/NullCorrelation/NullCorrelation_Simulation/'+filename

plt.show()
 




"""



