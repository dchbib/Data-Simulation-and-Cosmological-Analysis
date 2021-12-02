#!/usr/bin/python2.7.3
# -*- coding: latin-1 -*-

# Nom de fichier: Mock_Sample.py
#==============================================================================
#title           :Mock_Sample.py
#description     :This file is to create a mock sample.
#author          :Dyaa Chbib
#date            :2014_11_28
#version         :0.1
#python_version  :2.7.3
#==============================================================================

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
import socket

import Methods
from Methods import Functions
from Cosmological_Model import Model
from HubbleDiagram import ApparentMagnitude
from Savedata import Save

H_0 = np.float64(70.0) #Km.s-1.Mpc-1
LightVelocity = np.float64(2.99792458*10**5.0) # Km.s-1

class SIMULATED:
	"""Creating of the mock sample"""
	def __init__(self, lambda_0, omega_0):
	  
		self.param = Functions(lambda_0, omega_0)
                 
	def Sample(self, alfa, bbeta, s_0, c_0, sigmas, sigmac, M_0, sigmaM, sigma_eps, m_lim, zform, z_max, z_min, Sizeofsample, NumbOfSamples, Realdata, Direct_Indirect, degree_of_Polyn, Ncoefficients_of_KcorrPolyn):
	
		lambda0, omega0 = self.param.lambda_0, self.param.omega_0
		kappa0 = self.param.kappa_0
		print "kappa0 = ", kappa0, "\n"
		print Functions(lambda0, omega0).NoBigBang(omega0, 'permission')
		
		#Betalist = linspace(np.float64(0.01), np.float64(6), 600)
		Betalist = linspace(np.float64(0.01), np.float64(3.9), 150)
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

	        t = asc.read(Name_of_machine+'Zip_file/Python/Simulations/POO/OO/NullCorrelation/NullCorrelation_Realdata/SDSS-DR3/K-correction.dat', guess=False)

		degree_of_Polyn, Ncoefficients_of_KcorrPolyn = t['redsh'], t['kcorre']

		dkcorr = []

		for i in range(len(t['redsh'])):
			if i==0 or i==(len(t['redsh'])-1):
				dkcorr.append(t['kcorre'][i])
			else:
				dkcor = (t['kcorre'][i+1] - t['kcorre'][i-1])/(t['redsh'][i+1] - t['redsh'][i-1])
				dkcorr.append(dkcor)
				
		dkcorr = array(dkcorr)
		#====================================================================
		
		zz = np.linspace(z_min, zform, 10**4)	#10**5 #10**-8; For some models: if zz=0 --> z=1 --> tau=0 --> 1/tau = infini ou nan (informatic problem)
		#==============SUITE_of_K-corrections=======================================
		if Realdata == 'Realdata_NO':
			degree_of_Polyn_l, Ncoefficients_of_KcorrPolyn_l = [], []
			pass
		elif Realdata == 'Realdata_OK':
			List_of_Poly = []
			for i in range(int(len(t['kcorre'])/5)):
			    j = i*5
			    fit4 = np.polyfit(t['redsh'][j:j+5+1], t['kcorre'][j:j+5+1], int(len(t['kcorre'][j:j+5+1])))
			    fit_fn = np.poly1d(fit4)
			    List_of_Poly.append((min(t['redsh'][j:j+5+1]), max(t['redsh'][j:j+5+1]), fit_fn))

			    
			def Fit_of_Kcorrection(List_of_Poly, z):
				i = 0
				Condition = z>List_of_Poly[i][1]
				while Condition == True and i<=len(List_of_Poly)-2:
					i = i+1
					Condition = z>List_of_Poly[i][1]
				if Condition == True:
					value = -0.2135964
				else:
					value = List_of_Poly[i][2](z)
				return value
				
				
			degree_of_Polyn_l = []

			for i in range(len(zz)):                                                                                                      
			    if i==0 or i==(len(zz)-1):
				degree_of_Polyn_l.append(Fit_of_Kcorrection(List_of_Poly, zz[i]))
			    else:
				degree_of_Polyn_l.append((Fit_of_Kcorrection(List_of_Poly, zz[i+1])-Fit_of_Kcorrection(List_of_Poly, zz[i-1]))/(zz[i+1]-zz[i-1]))

			degree_of_Polyn_l = array(degree_of_Polyn_l)
			      
			Ncoefficients_of_KcorrPolyn_l = []
				
			for i in range(len(zz)):                                                                                                      
				Ncoefficients_of_KcorrPolyn_l.append(Fit_of_Kcorrection(List_of_Poly, zz[i]))
			    
			Ncoefficients_of_KcorrPolyn_l = array(Ncoefficients_of_KcorrPolyn_l)
		else:
			print "You must to determine the condition of the Realdata as 'Realdata_OK' or 'Realdata_NO'. ", "\n"
		#====================================================================
		ZETA = ApparentMagnitude(zz, lambda0, omega0, 0.0, Realdata, 0, 0).m_theor(zz, lambda0, omega0, 0.0, 0, 0, 0, 0)	
		
		paths_of_samples = []

		for bb in range(0, NumbOfSamples):
			
			# NOTE: if you make the test of hundred of simulation, so note the number of the simulation

			Number_of_simulation = bb

			#lists to save the complete sample, biased sample
			Mbias = []
			m_app, m_bias = [], []
			redshift, redshiftbias = [], []
			Volumebias = []
			redred  = []

			#================================Gaussian distribution of the absolute magnitude=================================

			print "Avant M", "\n"
			Mlist = np.random.normal(M_0, sigmaM, 2*Sizeofsample)
			print "Après M", "\n"
			
			#shape, scale = 2., 2.
			#Gamma_distrib = np.random.gamma(shape, scale, Sizeofsample)
			#Gammalist = []

			#========================================================================================================================

			max_of_F_M_V = 1
			print "max_of_F_M_V = ", max_of_F_M_V, "\n"
			max_of_F_M_V = [max_of_F_M_V] #*Sizeofsample

			#============================================================================================
			
			Mlist, Volume, Volume_l, redshift_l = [], [], [], []
			Volume_l1, Volume_l2 = [], []

			F_cap = np.random.uniform(low=0.0, high=1, size=(Sizeofsample))
			
			H = []
			H.append(10**-4.)
			for k in range(1, Sizeofsample):
				H.append(k*(1./Sizeofsample))
				
			HH = H
			#H = np.random.uniform(low=0.0, high=1, size=(Sizeofsample))
			#Mlist = np.random.normal(M_0, sigma, Sizeofsample)
			Mlistlist, slistlist, clistlist,  = [], [], []
			#Mlist = np.random.gumbel(M_0, sigma, Sizeofsample)
			#Mlist = np.random.triangular(-22, M_0, -16, Sizeofsample)
			shape, scale = 2., 2.
			#Mlist = -1*np.random.gamma(shape, scale, Sizeofsample)
			#=====generate the real distribution of Luminosity function=====
			
			s_lim, c_lim = 0.85, 0.14
			
			Mlist1 = np.linspace(-24, -16, 10)

			slist1 = np.linspace(-3, 3, 10)

			clist1 = np.linspace(-2, 2, 10)
			
			#slist1 = np.linspace(30, 60, 100)
			#clist1 = np.linspace(-2, 2, 100)
			#M_0, s_0, c_0 = -19.01, 40, 0
			#sigmas, sigmac = 1.5, 0.5
			
			"""
			p_th = dblquad(Functions(0.7, 0.3).Pro_DensityFunction_SN, -np.inf,  np.inf, lambda sy: -np.inf, lambda sy: np.inf, args=(M_0, c_0, 0, sigmaM, sigmac, bbeta, z_l, m_lim))[0]

			
			cdf_G_M_s_c = []
			for k in range(Sizeofsample):
			      G_M_s_c = dblquad(Functions(0.7, 0.3).Pro_DensityFunction_SN, -np.inf,  Mlist[k], lambda sy: -np.inf, lambda sy: clist[k], args=(M_0, c_0, 0, sigmaM, sigmac, bbeta, z_l, m_lim))[0]
			      cdf_G_M_s_c.append(G_M_s_c)
			
			
			
			F_M = array(F_M)/p_th
			for i in range(Sizeofsample):
				F_cap = np.random.uniform(low=0.0, high=1, size=(1))
				M_k = Functions(lambda0, omega0).Nonlinear_Interpolation(Mlist1, F_M, F_cap[i])
				Mlist.append(M_k)
			
				      
				      
			M_k = M_c - beta*col
			Functions(0.7, 0.3).Heavisidefuction(m_lim - M_k - zeta)
			
			
			
			"""
			
			"""
			SSlist, CClist = np.mgrid[min(slist1):max(slist1):10j, min(clist1):max(clist1):10j]
			
			CDF_G_M_s_c = []

			for j in range(len(slist1)):
				cdf_G_M_s_c = []
				for k in range(len(clist1)):
					#G_M_s_c = dblquad(Functions(0.7, 0.3).G_M_s_c_2D, min(slist1),  slist1[j], lambda sy: min(clist1), lambda sy: clist1[k], args=(M_0, s_0, c_0, sigmaM, sigmas, sigmac,alfa, bbeta))[0]
					M = M_0 - alfa*(slist1[j] - s_0) + bbeta*clist1[k]
					z_l = Functions(0.7, 0.3).Nonlinear_Interpolation(zz, ZETA, m_lim - M)
					#pd_of_V = quad(ApparentMagnitude(z_l, lambda0, omega0, 0.0, 'Realdata_NO', 0, 0).VOLUME, -30, M, args=(zz, ZETA, m_lim))[0]
					pd_of_V = 1.0
					G_M_s_c = dblquad(Functions(0.7, 0.3).Pro_DensityFunction_SN, min(slist1),  slist1[j], lambda sy: min(clist1), lambda sy: clist1[k], args=(M_0, s_0, c_0, 0, sigmas, sigmac,alfa, bbeta, z_l, m_lim))[0]
					cdf_G_M_s_c.append(G_M_s_c*pd_of_V)
				CDF_G_M_s_c.append(cdf_G_M_s_c)
				
				
			print "len(CDF_G_M_s_c) = ", len(CDF_G_M_s_c)
			for i in range(len(slist1)):
				print "len(CDF_G_M_s_c[",i,"]) = ", len(CDF_G_M_s_c[i]), "\n"

			CDF_G_M_s_c = array(CDF_G_M_s_c)
			CDF_of_G_M_s_c = CDF_G_M_s_c*1.0

			print "CDF_of_G_M_s_c = ", CDF_of_G_M_s_c
		      # # # Marginalization on s and c # # #
			Gs_tilda, s_s, Gc_tilda, c_c = [],[],[],[]
			for i in range(len(SSlist)):
				Gc_tilda.append(sum(CDF_of_G_M_s_c[i]))
				c_c.append(SSlist[i][0])
				Gs_tilda.append(sum(CDF_of_G_M_s_c[:,i]))                                                                                                               
				s_s.append(CClist[:,i][0])


			Gs_tilda, s_s, Gc_tilda, c_c = array(Gs_tilda)/max(Gs_tilda), array(s_s), array(Gc_tilda)/max(Gc_tilda), array(c_c)
			
			Mlist, slist, clist, Volume, Volume_l, redshift_l = [], [], [], [], [], []
			Volume_l1, Volume_l2 = [], []
			
			G_cap = np.random.uniform(low=0.0, high=1, size=(Sizeofsample))
			
			for i in range(Sizeofsample):
				Mlist.append(M_0)
				s_k = Functions(lambda0, omega0).Nonlinear_Interpolation(s_s, Gs_tilda, G_cap[i])
				slist.append(s_k)
				c_k = Functions(lambda0, omega0).Nonlinear_Interpolation(c_c, Gc_tilda, G_cap[i])
				clist.append(c_k)
			"""
			
			if Direct_Indirect == 'Direct':
				slist = np.random.normal(s_0, sigmas, Sizeofsample)
				epsilon = np.random.normal(0, sigma_eps, Sizeofsample) 
				Mlist = -alfa*slist + M_0 + epsilon
				clist = np.random.normal(c_0, sigmac, Sizeofsample)
				
				Mlist, slist, clist= array(Mlist), array(slist), array(clist) #*-1.0
				
				#computing of calibration's parameters
				print "=====Direct method====="
				COVA_riance = sum((1./(len(slist) - 1))*(slist - sum(slist/len(slist) ))*(Mlist - sum(Mlist/len(slist))))
				stand_M = sqrt(sum((1./(len(slist) - 1))*(Mlist - sum(Mlist/len(slist) ))*(Mlist - sum(Mlist/len(slist)))))
				stand_s = sqrt(sum((1./(len(slist) - 1))*(slist - sum(slist/len(slist) ))*(slist - sum(slist/len(slist)))))
				print "stand_s = ", stand_s, "\n"
				alfa_est = COVA_riance/stand_s**2.
				print "alfa_est = ", alfa_est, "\n"
				COR_elation = COVA_riance/(stand_s*stand_M)
				Mmen = sum(Mlist/len(slist))
				smen = sum(slist/len(slist))
				stand_xi = stand_M*sqrt(1-COR_elation**2.)
				print "stand_xi = ", stand_xi, "\n"
				M0_est = Mmen - alfa_est*smen + 3*stand_xi**2.
				print "M0_est = ", M0_est, "\n"
				print "=====Direct method====="

			elif Direct_Indirect == 'Indirect':
				Mlist = np.random.normal(M_0, sigmaM, Sizeofsample)
				epsilon = np.random.normal(0, sigma_eps, Sizeofsample)
				#slist = -alfa*Mlist + M_0 + epsilon  # Normalement il faut cette equation, c.a.d : slist = -alfa_IBS*Mlist + M_0_IBS + epsilon
				slist = -alfa*Mlist + (alfa*M_0) + epsilon
				clist = np.random.normal(c_0, sigmac, Sizeofsample)
				
				Mlist, slist, clist= array(Mlist), array(slist), array(clist) #*-1.0
				
				#computing of calibration's parameters
				print "=====Indirect method====="
				COVA_riance = sum((1./(len(slist) - 1))*(slist - sum(slist/len(slist) ))*(Mlist - sum(Mlist/len(slist))))
				stand_M = sqrt(sum((1./(len(slist) - 1))*(Mlist - sum(Mlist/len(slist) ))*(Mlist - sum(Mlist/len(slist)))))
				print "stand_M = ", stand_M, "\n"
				#alfa_est = stand_M**2./COVA_riance
				alfa_est = COVA_riance/stand_M**2. 
				print "alfa_est = ", alfa_est, "\n"
				stand_s = sqrt(sum((1./(len(slist) - 1))*(slist - sum(slist/len(slist) ))*(slist - sum(slist/len(slist)))))
				COR_elation = COVA_riance/(stand_s*stand_M)
				Mmen = sum(Mlist/len(slist))
				smen = sum(slist/len(slist))
				#stand_xi = stand_M*sqrt((1./COR_elation**2.) - 1)
				stand_xi = stand_s*sqrt(1-COR_elation**2.)
				print "stand_xi = ", stand_xi, "\n"
				M0_est = smen - alfa_est*Mmen
				print "M0_est = ", M0_est, "\n"

				print "=====Indirect method====="

			else: 
				print "You must to determine the method of simulation as 'Direct' or 'Indirect'. ", "\n"
			
			
			
			#slist = np.random.normal(s_0, sigmas, Sizeofsample)
			#clist = np.random.normal(c_0, sigmac, Sizeofsample)

			#@plt.figure(1)
			#@plt.hist(Mlist, bins=50, histtype='step', normed=False, color='b')
			#@plt.hist(slist, bins=50, histtype='step', normed=False, color='r')
			#@plt.hist(clist, bins=50, histtype='step', normed=False, color='g')
			
			#plt.figure(2)
			#plt.plot(s_s, Gs_tilda, color='r')
			#plt.plot(c_c, Gc_tilda, color='g')
			
			#fig = plt.figure(3)
			#ax = fig.add_subplot(111, projection='3d')
			#ax.plot_surface(CClist, SSlist, CDF_of_G_M_s_c, rstride=4, cstride=4, color='r')

			#@plt.figure(4)
			correlation_coefficient_s_c = np.corrcoef(slist, clist)
			print "correlation coefficient(s,c) = ", correlation_coefficient_s_c
			#@for i in range(Sizeofsample):
			#@	plt.plot(slist[i], clist[i], marker='+', color='r')
			
			#@plt.show()
			"""
			Mlistselection, slistselection, clistselection = [], [], []
			for i in range(len(slist)):
				if slist[i] < s_lim or clist[i] > c_lim:
					pass
				else:
					Mlistselection.append(Mlist[i])
					slistselection.append(slist[i])
					clistselection.append(clist[i])
				
			Mlist, slist, clist = array(Mlistselection), array(slistselection), array(clistselection)
			Sizeofsample = len(slist)
			"""
			#==================================================================

			if Realdata == 'Realdata_NO':
				
				if kappa0 > 0:
				  
					Angular_distance = Model(lambda0, omega0).curvature(np.linspace(z_min, 500, 10**4), lambda0, omega0, 0.0)[0]*sqrt(kappa0)
					z_pi = Functions(lambda0, omega0).Nonlinear_Interpolation(np.linspace(z_min, 500, 10**4), Angular_distance, pi-0.05) # redshift_of_antipode; En fait je cherche pi-0.05 parce que je ne veux pas que zform soit exactement le redshift d'antipode. Sinon, toute les objets vont etre compris entre M_star et Mform (ZETA_star et ZETA_form).
					if zform > z_pi:
						zform = z_pi*1.0
						zz = np.linspace(z_min, zform, 10**4)
						ZETA = ApparentMagnitude(zz, lambda0, omega0, 0.0, Realdata, 0, 0).m_theor(zz, lambda0, omega0, 0.0, 0, 0, 0, 0)
						
					else:
						pass
					
					M_star = m_lim - max(ZETA)
					ZETA_form = ApparentMagnitude(zform, lambda0, omega0, 0.0, Realdata, 0, 0).m_theor(zform, lambda0, omega0, 0.0, 0, 0, 0, 0)
					M_form = m_lim - ZETA_form
					
					for ii in range(Sizeofsample):
						M_k = Mlist[ii] + bbeta*clist[ii]
						s_k = slist[ii]
						c_k = clist[ii]
						#mu_lim = m_lim - M_k + alfa*(s_k - 1) - bbeta*c_k
						mu_lim = m_lim - M_k 
						
						if M_k>M_star and M_k<M_form:
						  
							ZETA_jj = ZETA*1.0
							zz_jj = zz*1.0
							z_l_jj = []
							for jj in range(10):
								z_ljj = Functions(lambda0, omega0).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim)
								z_l_jj.append(z_ljj)
								index_of_mu_lim = list(abs(ZETA_jj - mu_lim)).index(min(abs(ZETA_jj - mu_lim)))
								ZETA_jj = list(ZETA_jj*1)
								zz_jj = list(zz_jj*1)
								ZETA_jj.remove(ZETA_jj[index_of_mu_lim])
								zz_jj.remove(zz_jj[index_of_mu_lim])
								ZETA_jj = array(ZETA_jj*1)
								zz_jj = array(zz_jj*1)
								
							z_l_jj = array(z_l_jj)
							z_l_jj.sort()
							Distanguish_of_solutions = [str(gg)[0:5] for gg in z_l_jj]
							z_l_l = []
							z_l_l.append(z_l_jj[0])
							for i in range(1, len(z_l_jj)-1):
								condition = (Distanguish_of_solutions[i]==Distanguish_of_solutions[i-1])
								if condition == False:
									z_l_l.append(z_l_jj[i])
								else:
									pass
								
							Distanguish_of_solutions = [str(gg)[0:4] for gg in z_l_l]
							z_l = []
							z_l.append(z_l_l[0])
							for i in range(1, len(z_l_l)):
								condition1 = (Distanguish_of_solutions[i]==Distanguish_of_solutions[i-1])
								condition2 = (abs(float(Distanguish_of_solutions[i]) - float(Distanguish_of_solutions[i-1]))>= 0.06)
								if condition1 == False and condition2 == True:
									z_l.append(z_l_l[i])
								else:
									pass
							
							if (len(z_l) % 2 == 0): #even 
								print "je suis paire!!"
								print "len(z_l)", len(z_l)
								z_l.append(zform)
								z_l.sort()
								z_l = array(z_l)
								Volume_list = list(Model(lambda0, omega0).curvature(z_l, lambda0, omega0, 0.0)[1])
								V_l = Volume_list[0]
								RandomChoiceList_of_Volume = []
								for i in range(1, int((len(Volume_list)-1)/2)+1):
									V_l += Volume_list[2*i] - Volume_list[2*i-1]
									
								Volume_l.append(V_l)
								
								#V_K = HH[ii]*min(Volume_list) # V_z = V_K
								#V_K = HH[ii]*random.choice(Volume_list) # V_z = V_K
								
								V_K = HH[ii]*Volume_list[0] # V_z = V_K
								Volume.append(V_K)
								Mlistlist.append(Mlist[ii])
								slistlist.append(slist[ii])
								clistlist.append(clist[ii])
								
								V_K = Volume_list[1] + HH[ii]*(Volume_list[2]-Volume_list[1]) # V_z = V_K
								Volume.append(V_K)
								Mlistlist.append(Mlist[ii])
								slistlist.append(slist[ii])
								clistlist.append(clist[ii])
								
								"""
								ZoneProhibited_of_Volume = []
								ZoneProhibited_of_redshift = []
								for i in range(int(len(Volume_list)/2)):
									ZoneProhibited_of_Volume.append((Volume_list[2*i], Volume_list[2*i+1]))
									ZoneProhibited_of_redshift.append((z_l[2*i], z_l[2*i+1]))
									
								HH_list = []
								condition_on_HH = (len(HH_list) == 0)
								ij = -1
								V_z_list = []
								condition_on_V_z = (len(V_z_list) == 0)
								while condition_on_V_z == True and condition_on_HH == True:
									ij = ij + 1
									if ij > len(HH)-1:
										break
									else:
										V_K = HH[ij]*V_l
										V_z = V_K
										RandomChoiceList_of_Volume.append(V_z)
										for i in range(len(Volume_list)-1):
											if (i % 2 == 0): #even 
												V_z += - Volume_list[i]
											else: #odd
												V_z += Volume_list[i]
												RandomChoiceList_of_Volume.append(V_z)

										for dd in range(len(RandomChoiceList_of_Volume)):
											V_z = RandomChoiceList_of_Volume[dd]
											
											for i in range(len(ZoneProhibited_of_Volume)):
												condition = (V_z > min(ZoneProhibited_of_Volume[i]) and V_z < max(ZoneProhibited_of_Volume[i]))
												if condition == True:
													continue
												else:
													if V_z > max(ZoneProhibited_of_Volume[i]):
														continue 
													else: # ya3nei V_z as8ar mn Volume_list[2*i]
														V_z_list.append(V_z)
														HH_list.append(HH[ij])
														HH.remove(HH[ij])
														break
											condition_on_HH = (len(HH_list) == 0)
											condition_on_V_z = (len(V_z_list) == 0)
											if condition_on_V_z == True and condition_on_HH == True:
												continue
											else:
												pass
											
											V_z = V_z_list[0]
											Volume.append(V_z)
											Mlistlist.append(Mlist[ii])
											slistlist.append(slist[ii])
											clistlist.append(clist[ii])
											
								"""
									
							else: #odd
								print "je suis immpaire!!"
								z_l.sort()
								z_l = z_l[0]
								redshift_l.append(z_l)
								V_l = Model(lambda0, omega0).curvature(z_l, lambda0, omega0, 0.0)[1]
								Volume_l.append(V_l)
								V_min = Model(lambda0, omega0).curvature(z_min, lambda0, omega0, 0.0)[1]
								#V_K = np.random.uniform(low=V_min, high=V_l, size=(1))[0]
								#V_K = np.random.uniform(low=0, high=V_l, size=(1))[0]
								#V_K = H[i]*V_l
								#V_K = HH[0]*V_l  # V_K = V_z ds cet cas
								V_K = HH[ii]*V_l
								#HH.remove(HH[0])
								Volume.append(V_K)
								Mlistlist.append(Mlist[ii])
								slistlist.append(slist[ii])
								clistlist.append(clist[ii])
							
						else:
							print "je suis ni paire ni impaire!!"
							z_l = Functions(lambda0, omega0).Nonlinear_Interpolation(zz, ZETA, mu_lim)
							redshift_l.append(z_l)
							V_l = Model(lambda0, omega0).curvature(z_l, lambda0, omega0, 0.0)[1]
							Volume_l.append(V_l)
							V_min = Model(lambda0, omega0).curvature(z_min, lambda0, omega0, 0.0)[1]
							#V_K = np.random.uniform(low=V_min, high=V_l, size=(1))[0]
							#V_K = np.random.uniform(low=0, high=V_l, size=(1))[0]
							#V_K = H[i]*V_l
							#V_K = HH[0]*V_l  # V_K = V_z ds cet cas
							V_K = HH[ii]*V_l
							#HH.remove(HH[0])
							Volume.append(V_K)
							Mlistlist.append(Mlist[ii])
							slistlist.append(slist[ii])
							clistlist.append(clist[ii])
				  
				else:
					for ii in range(Sizeofsample):
						M_k = Mlist[ii] + bbeta*clist[ii]
						s_k = slist[ii]
						c_k = clist[ii]
						#mu_lim = m_lim - M_k + alfa*(s_k - 1) - bbeta*c_k
						mu_lim = m_lim - M_k 
						
						z_l = Functions(lambda0, omega0).Nonlinear_Interpolation(zz, ZETA, mu_lim)
						redshift_l.append(z_l)
						V_l = Model(lambda0, omega0).curvature(z_l, lambda0, omega0, 0.0)[1]
						Volume_l.append(V_l)
						V_min = Model(lambda0, omega0).curvature(z_min, lambda0, omega0, 0.0)[1]
						#V_K = np.random.uniform(low=V_min, high=V_l, size=(1))[0]
						#V_K = np.random.uniform(low=0, high=V_l, size=(1))[0]
						#V_K = H[i]*V_l
						V_K = HH[ii]*V_l
						#HH.remove(HH[0])
						Volume.append(V_K)
						Mlistlist.append(Mlist[ii])
						slistlist.append(slist[ii])
						clistlist.append(clist[ii])
			
			elif Realdata == 'Realdata_OK':
				
				if kappa0 > 0:
				  
					Angular_distance = Model(lambda0, omega0).curvature(np.linspace(z_min, 500, 10**4), lambda0, omega0, 0.0)[0]*sqrt(kappa0)
					z_pi = Functions(lambda0, omega0).Nonlinear_Interpolation(np.linspace(z_min, 500, 10**4), Angular_distance, pi-0.05) # redshift_of_antipode; En fait je cherche pi-0.05 parce que je ne veux pas que zform soit exactement le redshift d'antipode. Sinon, toute les objets vont etre compris entre M_star et Mform (ZETA_star et ZETA_form).
					if zform > z_pi:
						zform = z_pi*1.0
						zz = np.linspace(z_min, zform, 10**4)
						Ncoefficients_of_KcorrPolyn_l = array([Fit_of_Kcorrection(List_of_Poly, redsh) for redsh in zz])
						degree_of_Polyn_l = []
						for i in range(len(zz)):                                                                                                      
						    if i==0 or i==(len(zz)-1):
							degree_of_Polyn_l.append(Fit_of_Kcorrection(List_of_Poly, zz[i]))
						    else:
							degree_of_Polyn_l.append((Fit_of_Kcorrection(List_of_Poly, zz[i+1])-Fit_of_Kcorrection(List_of_Poly, zz[i-1]))/(zz[i+1]-zz[i-1]))

						degree_of_Polyn_l = array(degree_of_Polyn_l)
						ZETA = ApparentMagnitude(zz, lambda0, omega0, 0.0, Realdata, degree_of_Polyn_l, Ncoefficients_of_KcorrPolyn_l).m_theor(zz, lambda0, omega0, 0.0, 0, 0, 0, 0)
						
					else:
						pass
					
					M_star = m_lim - max(ZETA)
					Ncoefficients_of_KcorrPolyn_zform = Fit_of_Kcorrection(List_of_Poly,zform)
					degree_of_Polyn_zform = 0 # pas importante tant que je ne calcul pas le W_k
					ZETA_form = ApparentMagnitude(zform, lambda0, omega0, 0.0, Realdata, 0, Ncoefficients_of_KcorrPolyn_zform).m_theor(zform, lambda0, omega0, 0.0, 0, 0, 0, 0)
					M_form = m_lim - ZETA_form
					
					for ii in range(Sizeofsample):
						M_k = Mlist[ii] + bbeta*clist[ii]
						s_k = slist[ii]
						c_k = clist[ii]
						#mu_lim = m_lim - M_k + alfa*(s_k - 1) - bbeta*c_k
						mu_lim = m_lim - M_k 
						
						if M_k>M_star and M_k<M_form:
						  
							ZETA_jj = ZETA*1.0
							zz_jj = zz*1.0
							z_l_jj = []
							for jj in range(20):
								z_ljj = Functions(lambda0, omega0).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim)
								z_l_jj.append(z_ljj)
								index_of_mu_lim = list(abs(ZETA_jj - mu_lim)).index(min(abs(ZETA_jj - mu_lim)))
								ZETA_jj = list(ZETA_jj*1)
								zz_jj = list(zz_jj*1)
								ZETA_jj.remove(ZETA_jj[index_of_mu_lim])
								zz_jj.remove(zz_jj[index_of_mu_lim])
								ZETA_jj = array(ZETA_jj*1)
								zz_jj = array(zz_jj*1)
								
							z_l_jj = array(z_l_jj)
							z_l_jj.sort()
							Distanguish_of_solutions = [str(gg)[0:5] for gg in z_l_jj]
							z_l_l = []
							z_l_l.append(z_l_jj[0])
							for i in range(1, len(z_l_jj)-1):
								condition = (Distanguish_of_solutions[i]==Distanguish_of_solutions[i-1])
								if condition == False:
									z_l_l.append(z_l_jj[i])
								else:
									pass
								
							Distanguish_of_solutions = [str(gg)[0:4] for gg in z_l_l]
							z_l = []
							z_l.append(z_l_l[0])
							for i in range(1, len(z_l_l)):
								condition1 = (Distanguish_of_solutions[i]==Distanguish_of_solutions[i-1])
								condition2 = (abs(float(Distanguish_of_solutions[i]) - float(Distanguish_of_solutions[i-1]))>= 0.06)
								if condition1 == False and condition2 == True:
									z_l.append(z_l_l[i])
								else:
									pass
							
							if (len(z_l) % 2 == 0): #even 
								z_l.append(zform)
								z_l.sort()
								z_l = array(z_l)
								Volume_list = list(Model(lambda0, omega0).curvature(z_l, lambda0, omega0, 0.0)[1])
								V_l = Volume_list[0]
								RandomChoiceList_of_Volume = []
								for i in range(1, int((len(Volume_list)-1)/2)+1):
									V_l += Volume_list[2*i] - Volume_list[2*i-1]
									
								Volume_l.append(V_l)
								
								V_K = HH[ii]*min(Volume_list) # V_z = V_K
								Volume.append(V_K)
								Mlistlist.append(Mlist[ii])
								slistlist.append(slist[ii])
								clistlist.append(clist[ii])
								"""
								ZoneProhibited_of_Volume = []
								ZoneProhibited_of_redshift = []
								for i in range(int(len(Volume_list)/2)):
									ZoneProhibited_of_Volume.append((Volume_list[2*i], Volume_list[2*i+1]))
									ZoneProhibited_of_redshift.append((z_l[2*i], z_l[2*i+1]))
								
								HH_list = []
								condition_on_HH = (len(HH_list) == 0)
								ij = -1
								V_z_list = []
								condition_on_V_z = (len(V_z_list) == 0)
								while condition_on_V_z == True and condition_on_HH == True:
									ij = ij + 1
									if ij > len(HH)-1:
										break
									else:
										V_K = HH[ij]*V_l
										V_z = V_K
										RandomChoiceList_of_Volume.append(V_z)
										for i in range(len(Volume_list)-1):
											if (i % 2 == 0): #even 
												V_z += - Volume_list[i]
											else: #odd
												V_z += Volume_list[i]
												RandomChoiceList_of_Volume.append(V_z)

										for dd in range(len(RandomChoiceList_of_Volume)):
											V_z = RandomChoiceList_of_Volume[dd]
											
											for i in range(len(ZoneProhibited_of_Volume)):
												condition = (V_z > min(ZoneProhibited_of_Volume[i]) and V_z < max(ZoneProhibited_of_Volume[i]))
												if condition == True:
													continue
												else:
													if V_z > max(ZoneProhibited_of_Volume[i]):
														continue 
													else: # ya3nei V_z as8ar mn Volume_list[2*i]
														V_z_list.append(V_z)
														HH_list.append(HH[ij])
														HH.remove(HH[ij])
														break
											condition_on_HH = (len(HH_list) == 0)
											condition_on_V_z = (len(V_z_list) == 0)
											if condition_on_V_z == True and condition_on_HH == True:
												continue
											else:
												pass
											
											V_z = V_z_list[0]
											Volume.append(V_z)
											Mlistlist.append(Mlist[ii])
											slistlist.append(slist[ii])
											clistlist.append(clist[ii])
								 """
							else: #odd
								print "je suis immpaire!!"
								z_l.sort()
								z_l = z_l[0]
								redshift_l.append(z_l)
								V_l = Model(lambda0, omega0).curvature(z_l, lambda0, omega0, 0.0)[1]
								Volume_l.append(V_l)
								V_min = Model(lambda0, omega0).curvature(z_min, lambda0, omega0, 0.0)[1]
								#V_K = np.random.uniform(low=V_min, high=V_l, size=(1))[0]
								#V_K = np.random.uniform(low=0, high=V_l, size=(1))[0]
								#V_K = H[i]*V_l
								#V_K = HH[0]*V_l  # V_K = V_z ds cet cas
								V_K = HH[ii]*V_l
								#HH.remove(HH[0])
								Volume.append(V_K)
								Mlistlist.append(Mlist[ii])	
								slistlist.append(slist[ii])
								clistlist.append(clist[ii])
						else:
						  
							ZETA_jj = ZETA*1.0
							zz_jj = zz*1.0
							z_l_jj = []
							for jj in range(20):
								z_ljj = Functions(lambda0, omega0).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim)
								z_l_jj.append(z_ljj)
								index_of_mu_lim = list(abs(ZETA_jj - mu_lim)).index(min(abs(ZETA_jj - mu_lim)))
								ZETA_jj = list(ZETA_jj*1)
								zz_jj = list(zz_jj*1)
								ZETA_jj.remove(ZETA_jj[index_of_mu_lim])
								zz_jj.remove(zz_jj[index_of_mu_lim])
								ZETA_jj = array(ZETA_jj*1)
								zz_jj = array(zz_jj*1)
								
							z_l_jj = array(z_l_jj)
							z_l_jj.sort()
							Distanguish_of_solutions = [str(gg)[0:5] for gg in z_l_jj]
							z_l_l = []
							z_l_l.append(z_l_jj[0])
							for i in range(1, len(z_l_jj)-1):
								condition = (Distanguish_of_solutions[i]==Distanguish_of_solutions[i-1])
								if condition == False:
									z_l_l.append(z_l_jj[i])
								else:
									pass
								
							Distanguish_of_solutions = [str(gg)[0:4] for gg in z_l_l]
							z_l = []
							z_l.append(z_l_l[0])
							for i in range(1, len(z_l_l)):
								condition1 = (Distanguish_of_solutions[i]==Distanguish_of_solutions[i-1])
								condition2 = (abs(float(Distanguish_of_solutions[i]) - float(Distanguish_of_solutions[i-1]))>= 0.06)
								if condition1 == False and condition2 == True:
									z_l.append(z_l_l[i])
								else:
									pass
										      
							if (len(z_l) % 2 == 0): #even 
								z_l = z_l[0]
								redshift_l.append(z_l)
								V_l = Model(lambda0, omega0).curvature(z_l, lambda0, omega0, 0.0)[1]
								Volume_l.append(V_l)
								V_min = Model(lambda0, omega0).curvature(z_min, lambda0, omega0, 0.0)[1]
								#V_K = np.random.uniform(low=V_min, high=V_l, size=(1))[0]
								#V_K = np.random.uniform(low=0, high=V_l, size=(1))[0]
								#V_K = H[i]*V_l
								V_K = HH[ii]*V_l
								#HH.remove(HH[0])
								Volume.append(V_K)
								Mlistlist.append(Mlist[ii])
								slistlist.append(slist[ii])
								clistlist.append(clist[ii])
							else: #odd
							
								z_l.sort()
								z_l = array(z_l)
								Volume_list = list(Model(lambda0, omega0).curvature(z_l, lambda0, omega0, 0.0)[1])
								V_l = Volume_list[0]
								RandomChoiceList_of_Volume = []
								for i in range(1, int((len(Volume_list)-1)/2)+1):
									V_l += Volume_list[2*i] - Volume_list[2*i-1]
									
								Volume_l.append(V_l)
								
								V_K = HH[ii]*min(Volume_list) # V_z = V_K
								Volume.append(V_K)
								Mlistlist.append(Mlist[ii])
								slistlist.append(slist[ii])
								clistlist.append(clist[ii])
								"""
								ZoneProhibited_of_Volume = []
								ZoneProhibited_of_redshift = []
								for i in range(int(len(Volume_list)/2)):
									ZoneProhibited_of_Volume.append((Volume_list[2*i], Volume_list[2*i+1]))
									ZoneProhibited_of_redshift.append((z_l[2*i], z_l[2*i+1]))
								
								HH_list = []
								condition_on_HH = (len(HH_list) == 0)
								ij = -1
								V_z_list = []
								condition_on_V_z = (len(V_z_list) == 0)
								while condition_on_V_z == True and condition_on_HH == True:
									ij = ij + 1
									if ij > len(HH)-1:
										break
									else:
										V_K = HH[ij]*V_l
										V_z = V_K
										RandomChoiceList_of_Volume.append(V_z)
										for i in range(len(Volume_list)-1):
											if (i % 2 == 0): #even 
												V_z += - Volume_list[i]
											else: #odd
												V_z += Volume_list[i]
												RandomChoiceList_of_Volume.append(V_z)

										for dd in range(len(RandomChoiceList_of_Volume)):
											V_z = RandomChoiceList_of_Volume[dd]
											
											for i in range(len(ZoneProhibited_of_Volume)):
												condition = (V_z > min(ZoneProhibited_of_Volume[i]) and V_z < max(ZoneProhibited_of_Volume[i]))
												if condition == True:
													continue
												else:
													if V_z > max(ZoneProhibited_of_Volume[i]):
														continue 
													else: # ya3nei V_z as8ar mn Volume_list[2*i]
														V_z_list.append(V_z)
														HH_list.append(HH[ij])
														HH.remove(HH[ij])
														break
											condition_on_HH = (len(HH_list) == 0)
											condition_on_V_z = (len(V_z_list) == 0)
											if condition_on_V_z == True and condition_on_HH == True:
												continue
											else:
												pass
											
											V_z = V_z_list[0]
											Volume.append(V_z)
											Mlistlist.append(Mlist[ii])
											slistlist.append(slist[ii])
											clistlist.append(clist[ii])
								"""
	
				else:
					for ii in range(Sizeofsample):
						M_k = Mlist[ii] + bbeta*clist[ii]
						s_k = slist[ii]
						c_k = clist[ii]
						#mu_lim = m_lim - M_k + alfa*(s_k - 1) - bbeta*c_k
						mu_lim = m_lim - M_k 
						
						ZETA_jj = ZETA*1.0
						zz_jj = zz*1.0
						z_l_jj = []
						for jj in range(20):
							z_ljj = Functions(lambda0, omega0).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim)
							z_l_jj.append(z_ljj)
							index_of_mu_lim = list(abs(ZETA_jj - mu_lim)).index(min(abs(ZETA_jj - mu_lim)))
							ZETA_jj = list(ZETA_jj*1)
							zz_jj = list(zz_jj*1)
							ZETA_jj.remove(ZETA_jj[index_of_mu_lim])
							zz_jj.remove(zz_jj[index_of_mu_lim])
							ZETA_jj = array(ZETA_jj*1)
							zz_jj = array(zz_jj*1)
							
						z_l_jj = array(z_l_jj)
						z_l_jj.sort()
						Distanguish_of_solutions = [str(gg)[0:5] for gg in z_l_jj]
						z_l_l = []
						z_l_l.append(z_l_jj[0])
						for i in range(1, len(z_l_jj)-1):
							condition = (Distanguish_of_solutions[i]==Distanguish_of_solutions[i-1])
							if condition == False:
								z_l_l.append(z_l_jj[i])
							else:
								pass
							
						Distanguish_of_solutions = [str(gg)[0:4] for gg in z_l_l]
						z_l = []
						z_l.append(z_l_l[0])
						for i in range(1, len(z_l_l)):
							condition1 = (Distanguish_of_solutions[i]==Distanguish_of_solutions[i-1])
							condition2 = (abs(float(Distanguish_of_solutions[i]) - float(Distanguish_of_solutions[i-1]))>= 0.06)
							if condition1 == False and condition2 == True:
								z_l.append(z_l_l[i])
							else:
								pass
						
						if (len(z_l) % 2 == 0): #even 
							z_l = z_l[0]
							redshift_l.append(z_l)
							V_l = Model(lambda0, omega0).curvature(z_l, lambda0, omega0, 0.0)[1]
							Volume_l.append(V_l)
							V_min = Model(lambda0, omega0).curvature(z_min, lambda0, omega0, 0.0)[1]
							#V_K = np.random.uniform(low=V_min, high=V_l, size=(1))[0]
							#V_K = np.random.uniform(low=0, high=V_l, size=(1))[0]
							#V_K = H[i]*V_l
							V_K = HH[ii]*V_l
							#HH.remove(HH[0])
							Volume.append(V_K)
							Mlistlist.append(Mlist[ii])
							slistlist.append(slist[ii])
							clistlist.append(clist[ii])
						else: #odd
							z_l.sort()
							z_l = array(z_l)
							Volume_list = list(Model(lambda0, omega0).curvature(z_l, lambda0, omega0, 0.0)[1])
							V_l = Volume_list[0]
							RandomChoiceList_of_Volume = []
							for i in range(1, int((len(Volume_list)-1)/2)+1):
								V_l += Volume_list[2*i] - Volume_list[2*i-1]
								
							Volume_l.append(V_l)

							V_K = HH[ii]*min(Volume_list) # V_z = V_K
							Volume.append(V_K)
							Mlistlist.append(Mlist[ii])
							slistlist.append(slist[ii])
							clistlist.append(clist[ii])
							"""
							ZoneProhibited_of_Volume = []
							ZoneProhibited_of_redshift = []
							for i in range(int(len(Volume_list)/2)):
								ZoneProhibited_of_Volume.append((Volume_list[2*i], Volume_list[2*i+1]))
								ZoneProhibited_of_redshift.append((z_l[2*i], z_l[2*i+1]))
							
							HH_list = []
							condition_on_HH = (len(HH_list) == 0)
							ij = -1
							V_z_list = []
							condition_on_V_z = (len(V_z_list) == 0)
							while condition_on_V_z == True and condition_on_HH == True:
								ij = ij + 1
								if ij > len(HH)-1:
									break
								else:
									V_K = HH[ij]*V_l
									V_z = V_K
									RandomChoiceList_of_Volume.append(V_z)
									for i in range(len(Volume_list)-1):
										if (i % 2 == 0): #even 
											V_z += - Volume_list[i]
										else: #odd
											V_z += Volume_list[i]
											RandomChoiceList_of_Volume.append(V_z)
											
										for dd in range(len(RandomChoiceList_of_Volume)):
											V_z = RandomChoiceList_of_Volume[dd]
											
											for i in range(len(ZoneProhibited_of_Volume)):
												condition = (V_z > min(ZoneProhibited_of_Volume[i]) and V_z < max(ZoneProhibited_of_Volume[i]))
												if condition == True:
													continue
												else:
													if V_z > max(ZoneProhibited_of_Volume[i]):
														continue 
													else: # ya3nei V_z as8ar mn Volume_list[2*i]
														V_z_list.append(V_z)
														HH_list.append(HH[ij])
														HH.remove(HH[ij])
														break
											condition_on_HH = (len(HH_list) == 0)
											condition_on_V_z = (len(V_z_list) == 0)
											if condition_on_V_z == True and condition_on_HH == True:
												continue
											else:
												pass
											
											V_z = V_z_list[0]
											Volume.append(V_z)
											Mlistlist.append(Mlist[ii])
											slistlist.append(slist[ii])
											clistlist.append(clist[ii])
							"""
			else:
				print "You must to determine the condition of the Realdata as 'Realdata_OK' or 'Realdata_NO'. ", "\n"
			
			Mlistlist = np.asarray(Mlistlist)
			slistlist = np.asarray(slistlist)
			clistlist = np.asarray(clistlist)
			
			Mlist = Mlistlist*1.0
			slist = slistlist*1.0
			clist = clistlist*1.0
			#==================
			
			#slist = np.random.normal(1, sigmas, Sizeofsample)
			#clist = np.random.normal(0.01, sigmac, Sizeofsample)
			
			#==================
			Volume = np.asarray(Volume)
			Sizeofsample = len(Volume)
			#=====================Seek of ComovingDistance, redshift and computation of apparent magnitude========
			Mlist_second, slist_second, clist_second = [], [], []
			Volume_second = []
			Volume_l_second = []
			redshift_l_second = []
			V = Model(lambda0, omega0).curvature(zz, lambda0, omega0, 0.0)[1]
			print "len(Volume), len(Mlist) = ", len(Volume), len(Mlist), "\n"
			for ii in range(Sizeofsample):
				
				red = Functions(lambda0, omega0).Nonlinear_Interpolation(zz, V, Volume[ii])
				redshift.append(red)
				redred.append(red)
				Mlist_second.append(Mlist[ii])
				slist_second.append(slist[ii])
				clist_second.append(clist[ii])
				      
			Mlist = array(Mlist_second)
			slist = array(slist_second)
			clist = array(clist_second)
			
			Sizeofsample = len(redshift)
			redshift = array(redshift)
			lookbacktime = (978./H_0)*Model(lambda0, omega0).curvature(redshift, lambda0, omega0, 0.0)[3]
			lifetime, lookbacktime_1, lookbacktime_2 = [], [], []
			lookbacktime_of_redshift_atformation_epoch = (978./H_0)*Model(lambda0, omega0).curvature(16, lambda0, omega0, 0.0)[3]
			"""
			for iik in range(Sizeofsample):
				shape, scale = 5, 14
				#shape, scale = 100, 200
				etta = np.random.gamma(shape, scale, 6000)
				Date = np.random.uniform(low=0.0, high=lookbacktime_of_redshift_atformation_epoch, size=(6000))
				print 'min(Date), mean(Date), max(Date) = ', min(Date), mean(Date), max(Date), '\n'
				print 'min(lookbacktime), mean(lookbacktime), max(lookbacktime) = ', min(lookbacktime), mean(lookbacktime), max(lookbacktime), '\n'
				half_dt = (etta/(2.*365))*10**-9
				t_1 = Date - half_dt
				t_2 = Date + half_dt
				for ijk in range(len(Date)):
					if lookbacktime[iik]>t_1[ijk] and lookbacktime[iik]<t_2[ijk]:
						lifetime.append(etta[ijk])
						lookbacktime_1.append(t_1[ijk])
						lookbacktime_2.append(t_2[ijk])
						break
					else:
						pass
			"""
			for iik in range(Sizeofsample):
				shape, scale = 5, 14
				etta = np.random.gamma(shape, scale, 1)[0]	# En jours
				date_of_observation_factor = np.random.uniform(low=0.0, high=1, size=(1))[0]	# je prend un valeurs entre 0 et 1 pour que ne sera pas toutes les SNs detéctés dans la moitié de leurs vies (si je me limite sur la division de etta sur 2).
				part_1_dt = ((etta*date_of_observation_factor)/365.)*10**-9		# En Gyrs
				part_2_dt = ((etta*(1-date_of_observation_factor))/365.)*10**-9		# En Gyrs, (1-date_of_observation_factor) = (duré de vie de SN - date_of_observation_factor)
				t_1 = lookbacktime[iik] - part_1_dt		# En Gyrs
				t_2 = lookbacktime[iik] + part_2_dt		# En Gyrs
				lifetime.append(etta)
				lookbacktime_1.append(t_1)
				lookbacktime_2.append(t_2)
			
			zzz = np.linspace(z_min, zform, 10**4)
			Theoretical_lookbacktime = (978./H_0)*Model(lambda0, omega0).curvature(zzz, lambda0, omega0, 0.0)[3]
			redshift_1, redshift_2 = [], []
			for igk in range(len(lookbacktime_1)):
				red_1 = Functions(lambda0, omega0).Nonlinear_Interpolation(zzz, Theoretical_lookbacktime, lookbacktime_1[igk])
				redshift_1.append(red_1)
				red_2 = Functions(lambda0, omega0).Nonlinear_Interpolation(zzz, Theoretical_lookbacktime, lookbacktime_2[igk])
				redshift_2.append(red_2)
			
			lifetime, redshift_1, redshift_2, lookbacktime_1, lookbacktime_2 = array(lifetime), array(redshift_1), array(redshift_2), array(lookbacktime_1), array(lookbacktime_2)
			print "len(Volume), len(Mlist), len(lifetime) = ", len(Volume), len(Mlist), len(lifetime), "\n"
			#==============SUITE_of_K-corrections=======================================
			if Realdata == 'Realdata_NO':
				degree_of_Polyn, Ncoefficients_of_KcorrPolyn = [0]*Sizeofsample, [0]*Sizeofsample
				pass
			elif Realdata == 'Realdata_OK':
				kcorr = []
				d_kcorr = []
				for i in range(len(redshift)):
					#print i
					k_cor = Functions(lambda0, omega0).Nonlinear_Interpolation(t['kcorre'], t['redsh'], redshift[i])
					kcorr.append(k_cor)
					
				for i in range(len(redshift)):
					d_k_cor = Functions(lambda0, omega0).Nonlinear_Interpolation(dkcorr, t['redsh'], redshift[i])
					d_kcorr.append(d_k_cor)
					
					
				kcorr = array(kcorr)
				d_kcorr = array(d_kcorr)
				degree_of_Polyn, Ncoefficients_of_KcorrPolyn = d_kcorr, kcorr
			else:
				print "You must to determine the condition of the Realdata as 'Realdata_OK' or 'Realdata_NO'. ", "\n"
			#====================================================================

			observed_lifetime = (1 + redshift)*lifetime
			m_app = ApparentMagnitude(redshift, lambda0, omega0, Mlist, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(redshift, lambda0, omega0, Mlist, 0, bbeta, slist, clist)

			#=================================================
			redred.sort()
			Mlist = np.asarray(Mlist)
			Mbias = np.asarray(Mbias)
			m_app = np.asarray(m_app)
			m_bias = np.asarray(m_bias)

			redred = np.asarray(redred)
			redshift = np.asarray(redshift)
			redshiftbias = np.asarray(redshiftbias)

			print "min(redshift) = ", min(redshift)
			print "max(redshift) = ", max(redshift)
			print "max(m_app) = ", max(m_app)
			print "len(m_app) = ", len(m_app),"\n"
			print "len(m_bias) = ", len(m_bias),"\n"

			Mlist1 = Mlist*1.0
			m_app1 = m_app*1.0
			Mlist1.sort()
			m_app1.sort()

			#c'est pour faire un plot
			if z_min==0.0:
			    z1 = linspace(z_min+10**-6, z_max, Sizeofsample)
			else:
			    z1 = linspace(z_min, z_max, Sizeofsample)
			
			m_thlistplt = ApparentMagnitude(z1, lambda0, omega0, M_0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z1, lambda0, omega0, M_0, alfa, bbeta, s_0, c_0)
			print "len(m_thlistplt) = ", len(m_thlistplt),"\n"
			print "type(m_thlistplt) = ", type(m_thlistplt),"\n"
			print "dimension of m_thlistplt = ", m_thlistplt.shape,"\n"

			m_thlistplt2 = ApparentMagnitude(redshift, lambda0, omega0, M_0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(redshift, lambda0, omega0, M_0, alfa, bbeta, s_0, c_0)

			#==================================================================================
			
			#============determine of beta (constant used in the weighting factor term)=======================
			#@ten_power_m_M = (10**((m_app - bbeta*clist - Mlist)/5.))**Beta
			#@Jacob = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Jacobian(redshift)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
			#@WeightingManifold = Jacob * ten_power_m_M
			WeightingManifold = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, Beta)
			
			ten_power_s_m = (10**((bbeta*clist)/5.))**Beta
			ten_power_c_m = (10**((-alfa*slist)/5.))**Beta
			
			WeightingManifold_s = WeightingManifold*ten_power_s_m
			WeightingManifold_c = WeightingManifold*ten_power_c_m
						
			WeightingSum = [a/sum(a) for a in WeightingManifold]
			#DetermineBeta = [abs(max(b) - min(b)) for b in WeightingSum]
			DetermineBeta = [abs((max(b) - min(b))/(max(b) + min(b))) for b in WeightingSum]
			beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]

			#@ten_power_m_M_beta = (10**((m_app - bbeta*clist - Mlist)/5.))**beta
			#@Weighting = (Jacob*ten_power_m_M_beta)
			Weighting = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, beta)
			
			SumWeighting = sum(Weighting)
			factorWeighting = array(Weighting)/SumWeighting
			
			#----------------------pour slist -----------------------
			WeightingSum_s = [a/sum(a) for a in WeightingManifold_s]
			DetermineBeta = [abs(max(b) - min(b)) for b in WeightingSum_s]
			beta_s = Beta[DetermineBeta.index(min(DetermineBeta))][0]

			ten_power_s_m = (10**((bbeta*clist)/5.))**beta_s
			Weighting_s = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, beta_s)*ten_power_s_m
			SumWeighting = sum(Weighting_s)
			factorWeighting_s = array(Weighting_s)/SumWeighting
			#---------------------------------------------------------
			
			#----------------------pour clist -----------------------
			WeightingSum_c = [a/sum(a) for a in WeightingManifold_c]
			DetermineBeta = [abs(max(b) - min(b)) for b in WeightingSum_c]
			beta_c = Beta[DetermineBeta.index(min(DetermineBeta))][0]
			
			ten_power_c_m = (10**((-alfa*slist)/5.))**beta_c
			Weighting_c = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, beta_c)*ten_power_c_m
			SumWeighting = sum(Weighting_c)
			factorWeighting_c = array(Weighting_c)/SumWeighting
			#---------------------------------------------------------

			#===================================================================================================

			#=======save all data to (plot, ...) in /media/Elements/Zip_file/Python/Simulations/POO/OO/NullCorrelation/NullCorrelation_Simulation/Savedata.dat========

			#filename = 'SecondHundredSimulation/SimulationData_'+str(lambda0)+'-'+str(omega0)+'_NbOfObjects'+str(Sizeofsample)+'_limitMag'+str(m_lim)+'_NbOfSimulation'+str(Number_of_simulation)+'.dat'
			#filename = 'OneHundredSimulation/OneHundredSimulation_2000_Objects/SimulationData_'+str(lambda0)+'-'+str(omega0)+'_NbOfObjects'+str(Sizeofsample)+'_limitMag'+str(m_lim)+'_NbOfSimulation'+str(Number_of_simulation)+'.dat'
			#filename = 'SimulationData_'+str(lambda0)+'-'+str(omega0)+'_NbOfObjects'+str(Sizeofsample)+'_limitMag'+str(m_lim)+'_NbOfSimulation'+str(Number_of_simulation)+'.dat'
			#filename = 'Grid_around_model/SimulationData_'+str(lambda0)+'-'+str(omega0)+'_NbOfObjects'+str(Sizeofsample)+'_limitMag'+str(m_lim)+'_NbOfSimulation'+str(Number_of_simulation)+'.dat'
			#filename = 'OneHundredSimulation/OneHundredSimulation_350_Objects/SimulationData_'+str(lambda0)+'-'+str(omega0)+'_NbOfObjects'+str(Sizeofsample)+'_limitMag'+str(m_lim)+'_NbOfSimulation'+str(Number_of_simulation)+'Beta'+str(beta)+'.dat'
			#filename = 'OneHundredSimulation/Different_Number_of_Simulations_350_Objects/SimulationData_'+str(lambda0)+'-'+str(omega0)+'_NbOfObjects'+str(Sizeofsample)+'_limitMag'+str(m_lim)+'_NbOfSimulation'+str(Number_of_simulation)+'Beta'+str(beta)+'.dat'
			filename = 'SimulatedData_'+str(lambda0)+'-'+str(omega0)+'_alfa'+str(alfa)+'_bbeta'+str(bbeta)+'_s0'+str(s_0)+'_c0'+str(c_0)+'_sigmas'+str(sigmas)+'_sigmac'+str(sigmac)+'_M0'+str(M_0)+'_NbOfObjects'+str(Sizeofsample)+'_limitMag'+str(m_lim)+'_NbOfSimulation'+str(Number_of_simulation)+'_Beta'+str(beta)+'.dat'
			
			listofarrays = [slist, clist, Mlist, m_app, redshift, redshift_l, Volume, Volume_l, Volumebias, Mbias, m_bias, redshiftbias, lifetime, lookbacktime, redshift_1, redshift_2, lookbacktime_1, lookbacktime_2, Weighting, Weighting_s, Weighting_c, factorWeighting_s, factorWeighting_c, factorWeighting, redred, z1, m_thlistplt, degree_of_Polyn, Ncoefficients_of_KcorrPolyn, max_of_F_M_V]
			stringlistofarrays = ['slist', 'clist', 'Mlist', 'm_app', 'redshift', 'redshift_l', 'Volume', 'Volume_l', 'Volumebias', 'Mbias', 'm_bias', 'redshiftbias', 'lifetime', 'lookbacktime', 'redshift_1', 'redshift_2', 'lookbacktime_1', 'lookbacktime_2', 'Weighting', 'Weighting_s', 'Weighting_c', 'factorWeighting_s', 'factorWeighting_c', 'factorWeighting', 'redred', 'z1', 'm_thlistplt', 'degree_of_Polyn', 'Ncoefficients_of_KcorrPolyn', 'max_of_F_M_V']

			Save(listofarrays, stringlistofarrays, filename).ReadWrite(True, 'work')
			path = Save(listofarrays, stringlistofarrays, filename).ReadWrite(True, 'work')
			paths_of_samples.append(path)
			#============================================================================================================
			print "Number_of_simulation = ", Number_of_simulation, "\n"
			print "Data set in this path: ", "\n"
			print "path = ", path, "\n"
			
		return paths_of_samples
			


"""
### ### ### ### ### ### ### ### ###Pour plusieurs test de calibration ### ### ### ### ### ### ### ### ###

Direct_Indirect = 'Direct'
Direct_Indirect = 'Indirect'
for i in range(10):
	if Direct_Indirect == 'Direct':
		slist = np.random.normal(s0, sigmas, Sizeofsample)
		epsilon = np.random.normal(0, sigma_eps, Sizeofsample) 
		Mlist = -alfa*slist + M0 + epsilon
		clist = np.random.normal(c0, sigmac, Sizeofsample)
		
		Mlist, slist, clist= array(Mlist), array(slist), array(clist) #*-1.0
		
		#computing of calibration's parameters
		print "=====Direct method====="
		COVA_riance = sum((1./(len(slist) - 1))*(slist - sum(slist/len(slist) ))*(Mlist - sum(Mlist/len(slist))))
		stand_M = sqrt(sum((1./(len(slist) - 1))*(Mlist - sum(Mlist/len(slist) ))*(Mlist - sum(Mlist/len(slist)))))
		stand_s = sqrt(sum((1./(len(slist) - 1))*(slist - sum(slist/len(slist) ))*(slist - sum(slist/len(slist)))))
		print "stand_s = ", stand_s, "\n"
		alfa_est = COVA_riance/stand_s**2.
		print "alfa_est = ", alfa_est, "\n"
		COR_elation = COVA_riance/(stand_s*stand_M)
		Mmen = sum(Mlist/len(slist))
		smen = sum(slist/len(slist))
		stand_xi = stand_M*sqrt(1-COR_elation**2.)
		print "stand_xi = ", stand_xi, "\n"
		M0_est = Mmen - alfa_est*smen + 3*stand_xi**2.
		print "M0_est = ", M0_est, "\n"
		print "=====Direct method====="

	elif Direct_Indirect == 'Indirect':
		Mlist = np.random.normal(M0, sigmaM, Sizeofsample)
		epsilon = np.random.normal(0, sigma_eps, Sizeofsample)
		#slist = -alfa*Mlist + M0 + epsilon  # Normalement il faut cette equation, c.a.d : slist = -alfa_IBS*Mlist + M_0_IBS + epsilon
		slist = -alfa*Mlist + 2*M0 + epsilon
		clist = np.random.normal(c0, sigmac, Sizeofsample)
		
		Mlist, slist, clist= array(Mlist), array(slist), array(clist) #*-1.0
		
		#computing of calibration's parameters
		print "=====Indirect method====="
		COVA_riance = sum((1./(len(slist) - 1))*(slist - sum(slist/len(slist) ))*(Mlist - sum(Mlist/len(slist))))
		stand_M = sqrt(sum((1./(len(slist) - 1))*(Mlist - sum(Mlist/len(slist) ))*(Mlist - sum(Mlist/len(slist)))))
		print "stand_M = ", stand_M, "\n"
		alfa_est = COVA_riance/stand_M**2.
		print "alfa_est = ", alfa_est, "\n"
		stand_s = sqrt(sum((1./(len(slist) - 1))*(slist - sum(slist/len(slist) ))*(slist - sum(slist/len(slist)))))
		COR_elation = COVA_riance/(stand_s*stand_M)
		Mmen = sum(Mlist/len(slist))
		smen = sum(slist/len(slist))
		#stand_xi = stand_M*sqrt((1./COR_elation**2.) - 1)
		stand_xi = stand_s*sqrt(1-COR_elation**2.)
		print "stand_xi = ", stand_xi, "\n"
		M0_est = smen - alfa_est*Mmen
		print "M0_est = ", M0_est, "\n"
		print "=====Indirect method====="
	else: 
		print "You must to determine the method of simulation as 'Direct' or 'Indirect'. ", "\n"

		
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


#computing of calibration's parameters
print "=====Direct method====="
COVA_riance = sum((1./(len(slist) - 1))*(slist - sum(slist/len(slist) ))*(Mlist - sum(Mlist/len(slist))))
stand_M = sqrt(sum((1./(len(slist) - 1))*(Mlist - sum(Mlist/len(slist) ))*(Mlist - sum(Mlist/len(slist)))))
stand_s = sqrt(sum((1./(len(slist) - 1))*(slist - sum(slist/len(slist) ))*(slist - sum(slist/len(slist)))))
print "stand_s = ", stand_s, "\n"
alfa_est = COVA_riance/stand_s**2.
print "alfa_est = ", alfa_est, "\n"
COR_elation = COVA_riance/(stand_s*stand_M)
Mmen = sum(Mlist/len(slist))
smen = sum(slist/len(slist))
stand_xi = stand_M*sqrt(1-COR_elation**2.)
print "stand_xi = ", stand_xi, "\n"
M0_est = Mmen - alfa_est*smen + 3*stand_xi**2.
print "M0_est = ", M0_est, "\n"
print "=====Direct method====="

#computing of calibration's parameters
print "=====Indirect method====="
COVA_riance = sum((1./(len(slist) - 1))*(slist - sum(slist/len(slist) ))*(Mlist - sum(Mlist/len(slist))))
stand_M = sqrt(sum((1./(len(slist) - 1))*(Mlist - sum(Mlist/len(slist) ))*(Mlist - sum(Mlist/len(slist)))))
print "stand_M = ", stand_M, "\n"
alfa_est = COVA_riance/stand_M**2.
print "alfa_est = ", alfa_est, "\n"
stand_s = sqrt(sum((1./(len(slist) - 1))*(slist - sum(slist/len(slist) ))*(slist - sum(slist/len(slist)))))
COR_elation = COVA_riance/(stand_s*stand_M)
Mmen = sum(Mlist/len(slist))
smen = sum(slist/len(slist))
#stand_xi = stand_M*sqrt((1./COR_elation**2.) - 1)
stand_xi = stand_s*sqrt(1-COR_elation**2.)
print "stand_xi = ", stand_xi, "\n"
M0_est = smen - alfa_est*Mmen
print "M0_est = ", M0_est, "\n"
print "=====Indirect method====="



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#============determine of beta (constant used in the weighting factor term)=======================
Betalist = linspace(np.float64(0.01), np.float64(3.9), 350)
Beta = array([array([a]) for a in Betalist])
#WeightingManifold = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, Beta)
ten_power_m_M = (10**((m_app - bbeta*clist - Mlist)/5.))**Beta
Jacob = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Jacobian(redshift)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
WeightingManifold = Jacob * ten_power_m_M
WeightingSum = [a/sum(a) for a in WeightingManifold]
#DetermineBeta = [abs(max(b) - min(b)) for b in WeightingSum]
DetermineBeta = [abs((max(b) - min(b))/(max(b) + min(b))) for b in WeightingSum]
beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
#Weighting = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, beta)
ten_power_m_M_beta = (10**((m_app - bbeta*clist - Mlist)/5.))**beta
Weighting = (Jacob*ten_power_m_M_beta)
SumWeighting = sum(Weighting)
factorWeighting = array(Weighting)/SumWeighting
#===================================================================================================
plt.figure()
plt.subplot(223)
for i in range(len(redshift)):
         plt.plot(redshift[i], 239*factorWeighting[i], marker='+', color='r')

plt.xlabel('$redshift$'+' '+'$z$', fontsize=13)
plt.ylabel("$\\omega_{k}$", fontsize=16)
Weighting = (Jacob*1.0)
SumWeighting = sum(Weighting)
factorWeighting = array(Weighting)/SumWeighting
#===================================================================================================
plt.subplot(221)
for i in range(len(redshift)):
         plt.plot(redshift[i], 239*factorWeighting[i], marker='+', color='b')

plt.xlabel('$redshift$'+' '+'$z$', fontsize=13)
plt.ylabel("$\\rho^{-1}$", fontsize=16)
WeightingManifold = ten_power_m_M*1.0
WeightingSum = [a/sum(a) for a in WeightingManifold]
#DetermineBeta = [abs(max(b) - min(b)) for b in WeightingSum]
DetermineBeta = [abs((max(b) - min(b))/(max(b) + min(b))) for b in WeightingSum]
beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
#Weighting = ApparentMagnitude(redshift, lambda0, omega0, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, beta)
ten_power_m_M_beta = (10**((m_app - bbeta*clist - Mlist)/5.))**beta
Weighting = (ten_power_m_M_beta*1.0)
SumWeighting = sum(Weighting)
factorWeighting = array(Weighting)/SumWeighting
#===================================================================================================
plt.subplot(222)
for i in range(len(redshift)):
         plt.plot(redshift[i], 239*factorWeighting[i], marker='+', color='g')

plt.xlabel('$redshift$'+' '+'$z$', fontsize=13)
plt.ylabel("$10^\\frac{\\mu \\zeta_{est}}{5}$", fontsize=16)



"""