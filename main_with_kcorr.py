#!/usr/bin/python2.7.3
# -*- coding: latin-1 -*-

#==============================================================================
#title           :main.py
#description     :It is the main file
#author          :Dyaa Chbib
#date            :2014_11_28
#version         :0.1
#python_version  :2.7.3
#==============================================================================

import time
print time.ctime()


import Methods
from Methods import Functions
from Cosmological_Model import Model
from HubbleDiagram import ApparentMagnitude
from GradientMethod3 import Gradient
#from GradientMethod4 import Gradient
#from Grid_Class import Gradient
from Savedata import Save
from Plot import Plotting
from Mock_Sample import SIMULATED
from Mock_Sample_SN_2D_CDF import SIMULATED

lambda0, omega0, M_0, sigma, m_lim, zform, z_max, z_min, Sizeofsample, NumbOfSamples = 0.8, 0.4, -19.7, 0.3, 26, 16, 16, 10**-15.95, 300, 1


omegalist = linspace(np.float64(10**-6), np.float64(1.0), 30) 
lamdalist = linspace(np.float64(10**-6), np.float64(2.0), 30) 

M_0list = linspace(np.float64(-21), np.float64(-18), 30) 

alfalist = linspace(np.float64(10**-6), np.float64(5), 10) 
bbetalist = linspace(np.float64(10**-6), np.float64(5), 10) 

All_models = []
for i in range(len(alfalist)):
	for j in range(len(bbetalist)):
		All_models.append((alfalist[i], bbetalist[j]))
		
All_models = array(All_models)		
#=======================================================================================================================
zz = np.linspace(z_min, z_max, 10**3)	#10**5 #10**-8; For some models: if zz=0 --> z=1 --> tau=0 --> 1/tau = infini ou nan (informatic problem)

ZETA = ApparentMagnitude(zz, lambda0, omega0, 0.0, 'Realdata_NO', 0, 0).m_theor(zz, lambda0, omega0, 0,0,0,0,0)

kappa0 = Functions(lambda0, omega0).kappa_0

#=================================================================================================================


lambda0, omega0, M_0, sigmaM, m_lim, zform, z_max, z_min, Sizeofsample, NumbOfSamples = 0.7, 0.3, -19.01, 0.3, 25.04, 16, 16, 10**-15.95, 240, 1

alfa, bbeta = 0.141, 3.14
Realdata = 'Realdata_NO'
degree_of_Polyn, Ncoefficients_of_KcorrPolyn = 0, 0
M0, s0, c0 = -19.00, 0, 0

sigmaM, sigmas, sigmac, sigma_eps = 0.25, 0.99, 0.085, 0.09

paths_of_samples = SIMULATED(lambda0, omega0).Sample(alfa, bbeta, s0, c0, sigmas, sigmac, M0, sigmaM, sigma_eps, m_lim, zform, z_max, z_min, Sizeofsample, NumbOfSamples, 'Realdata_NO', 'Direct', degree_of_Polyn, Ncoefficients_of_KcorrPolyn)
plt.show()

COR = []
for bb in range(0, NumbOfSamples):

	Number_of_simulation = bb
	path_of_sample = paths_of_samples[bb]

	t = asc.read(path_of_sample, guess=False)

	slist, clist, Mlist, m_app, redshift, redshift_l, Volume, Volume_l, Volumebias, Mbias, m_bias, redshiftbias, lifetime, lookbacktime, redshift_1, redshift_2, lookbacktime_1, lookbacktime_2, factorWeighting_s, factorWeighting_c, redred, z1, m_thlistplt, degree_of_Polyn, Ncoefficients_of_KcorrPolyn = t['slist'], t['clist'], t['Mlist'], t['m_app'], t['redshift'], t['redshift_l'], t['Volume'], t['Volume_l'], t['Volumebias'], t['Mbias'], t['m_bias'], t['redshiftbias'], t['lifetime'], t['lookbacktime'], t['redshift_1'], t['redshift_2'], t['lookbacktime_1'], t['lookbacktime_2'], t['factorWeighting_s'], t['factorWeighting_c'], t['redred'], t['z1'], t['m_thlistplt'], t['degree_of_Polyn'], t['Ncoefficients_of_KcorrPolyn']
	factorWeighting = t['factorWeighting']
	MM = Mlist + bbeta*clist

        #=====Weighting=====================================        
        Betalist = linspace(np.float64(0.01), np.float64(3.9), 150)
        Beta = array([array([a]) for a in Betalist])
	WeightingManifold = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, Beta)
	

	WeightingSum = [a/sum(a) for a in WeightingManifold]
	DetermineBeta = [abs((max(b) - min(b))/(max(b) + min(b))) for b in WeightingSum]
	beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]

	Weighting = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, beta)

	SumWeighting = sum(Weighting)
	factorWeighting = array(Weighting)/SumWeighting
        #=====Weighting=====================================
        
	#=================================================================
	#@m_app = m_app - bbeta*clist
	Mmean = sum(factorWeighting*MM)
	Mmean_m = sum(factorWeighting*m_app)

	Cov_M = sum(factorWeighting*((MM - Mmean)*(m_app - Mmean_m)))
	Cor_M = Cov_M

	#============================================================================================================
	
	#=================================================================
	Mmean = sum(factorWeighting_c*clist)
	Mmean_m = sum(factorWeighting_c*m_app)

	Cov_c = sum(factorWeighting_c*((clist - Mmean)*(m_app - Mmean_m)))
	Cor_c = Cov_c
	#============================================================================================================
	 
	#=================================================================
	smean = sum(factorWeighting*slist)
	Mmean_m = sum(factorWeighting*m_app)

	Cov_s = sum(factorWeighting*((slist - smean)*(m_app - Mmean_m)))
	print Cov_s
	#============================================================================================================
	 
	
