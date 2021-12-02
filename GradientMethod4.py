#!/usr/bin/python2.7.3
# -*- coding: latin-1 -*-

# Nom de fichier: GradientMethod.py
#==============================================================================
#title           :GradientMethod.py
#description     :This file is to find the null correlation curve on the grid of the all models.
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
from scipy.integrate import quad
from scipy import stats
from scipy.stats import norm
#from decimal import *


import Methods
from Methods import Functions
from Cosmological_Model import Model
from HubbleDiagram import ApparentMagnitude


class Gradient:
	"""Gradient Method to seekin of cosmological parameters"""
	def __init__(self):
		#self, Booleen
		#self.Booleen = Booleen
		pass

		
	def SeekOfParameters(self, zz, ZETA, zform, m_lim, m_app, redshift, M_0, alfa, bbeta, slist, clist, MM, Cov_s, Cov_c, lambda0_Model, omega0_Model, couleur, KS, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn, Scanning_on_parameters, Variables): # fobs, ftheo are array type. If you use a number, we must change chi2 = sum(.....

		#m_app = m_thlistplt2*1.0
		ZETA_virt = m_app - MM
		#******************************compute of all weighting factors*************************

		omegalist = linspace(np.float64(10**-6), np.float64(1.0), 30) 
		lamdalist = linspace(np.float64(10**-6), np.float64(2.0), 30) 
		alfalist = linspace(np.float64(2), np.float64(2.5), 15) 
		bbetalist = linspace(np.float64(2.8), np.float64(3.5), 30) 
		M_0list = linspace(np.float64(-21), np.float64(-18), 20) 
		
		#Betalist = linspace(np.float64(2.1), np.float64(3.1), 100) #100
		Betalist = linspace(np.float64(0.01), np.float64(3.9), 150)
		#Betalist = linspace(np.float64(0.01), np.float64(3.9), 100)
		
		Beta = array([array([a]) for a in Betalist])
		
		if Scanning_on_parameters == 'lambda0_omega0':
			if Variables == 'stretch_and_color':
				pass
				
			elif Variables == 'stretch_only':
				pass

			elif Variables == 'color_only':
				All_models_lambda_omega, All_models_omega_bbeta, All_models_bbeta_lambda = [], [], []
				
				
				for i in range(len(lamdalist)):
					for j in range(len(omegalist)):
						All_models_lambda_omega.append((lamdalist[i], omegalist[j]))
						
				for i in range(len(lamdalist)):
					for j in range(len(omegalist)):
						All_models_omega_bbeta.append((bbetalist[i], omegalist[j]))
						
				for i in range(len(lamdalist)):
					for j in range(len(omegalist)):
						All_models_bbeta_lambda.append((lamdalist[i], bbetalist[j]))
						
				
				m_app_no_correct = m_app*1.0
				X_omega, Y_lambda, Z_bbeta = [], [], []

				def func(a, z, bet):     # a = All_models;  z = redhsift;  bet = Beta
					aa=[]
					#beta = bet
					for i in range(len(a)):
						
						Bounced_limit = Functions(a[i][0],a[i][1]).NoBigBang(a[i][1], 'permission')       # si n'est pas valable, alors sera remplacer par le Modele standard (à voir après)
						#if lamdalist[i] >= Bounced_limit:
						if Bounced_limit==0:
							aa.append(Cov_c)
						else:
							#Mlist_tild = m_app - ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,bbeta,0,0)
							Mlist_tild = MM*1.0
							#ZETA_tild = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,0,0,0)
							
							ten_power_m_M = (10**((m_app - Mlist_tild)/5.))**bet
							Jacob = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Jacobian(z)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
							WeightingManifold = Jacob * ten_power_m_M
							#WeightingManifold = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, bet)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
							
							WeightingSum = [b/sum(b) for b in WeightingManifold]
							DetermineBeta = [abs((max(b) - min(b))/(max(b) + min(b))) for b in WeightingSum]
							beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
							#L_1_Manifold = [1.0 + (1./np.log(len(b))*sum(b*np.log(b))) for b in WeightingSum]
							#beta = Beta[L_1_Manifold.index(min(L_1_Manifold))][0]
							beta_pliot = beta*1.0

							#w_k = (ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta))/sum((ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta)))
							ten_power_m_M_beta = (10**((m_app - Mlist_tild)/5.))**beta
							w_k = (Jacob*ten_power_m_M_beta)/sum(Jacob*ten_power_m_M_beta)

							COVA_riance = sum(w_k*(Mlist_tild - sum(w_k*Mlist_tild ))*(m_app - sum(w_k*m_app)))

							COR_elation = COVA_riance
							#COR_elation = (1./((len(Mlist_tild)/(len(Mlist_tild) - 1.))*sqrt(sum(w_k*(Mlist_tild - sum(w_k*Mlist_tild ))**2.0))*sqrt(sum(w_k*(m_app - sum(w_k*m_app))**2.0))))*COVA_riance

							aa.append(COR_elation)

					return aa

				def func1(lambdaAxe, a, z, bet):     # a = All_models;  z = redhsift;  bet = Beta
					aa=[]
					#beta = bet
					for i in range(len(a)):
						m_app = m_app_no_correct - a[i][0]*clist
						Bounced_limit = Functions(lambdaAxe, a[i][1]).NoBigBang(a[i][1], 'permission')       # si n'est pas valable, alors sera remplacer par le Modele standard (à voir après)
						#if lamdalist[i] >= Bounced_limit:
						if Bounced_limit==0:
							aa.append(Cov_c)
						else:
							#Mlist_tild = m_app - ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,bbeta,0,0)
							Mlist_tild = MM*1.0
							#ZETA_tild = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,0,0,0)
							
							ten_power_m_M = (10**((m_app - Mlist_tild)/5.))**bet
							Jacob = ApparentMagnitude(z, lambdaAxe, a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Jacobian(z)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
							WeightingManifold = Jacob * ten_power_m_M
							#WeightingManifold = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, bet)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
							
							WeightingSum = [b/sum(b) for b in WeightingManifold]
							DetermineBeta = [abs((max(b) - min(b))/(max(b) + min(b))) for b in WeightingSum]
							beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
							#L_1_Manifold = [1.0 + (1./np.log(len(b))*sum(b*np.log(b))) for b in WeightingSum]
							#beta = Beta[L_1_Manifold.index(min(L_1_Manifold))][0]
							beta_pliot = beta*1.0

							#w_k = (ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta))/sum((ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta)))
							ten_power_m_M_beta = (10**((m_app - Mlist_tild)/5.))**beta
							w_k = (Jacob*ten_power_m_M_beta)/sum(Jacob*ten_power_m_M_beta)

							COVA_riance = sum(w_k*(Mlist_tild - sum(w_k*Mlist_tild ))*(m_app - sum(w_k*m_app)))

							COR_elation = COVA_riance
							#COR_elation = (1./((len(Mlist_tild)/(len(Mlist_tild) - 1.))*sqrt(sum(w_k*(Mlist_tild - sum(w_k*Mlist_tild ))**2.0))*sqrt(sum(w_k*(m_app - sum(w_k*m_app))**2.0))))*COVA_riance

							aa.append(COR_elation)

					return aa

				def func2(omegaAxe, a, z, bet):     # a = All_models;  z = redhsift;  bet = Beta
					aa=[]
					#beta = bet
					for i in range(len(a)):
						
						m_app = m_app_no_correct - a[i][1]*clist
						Bounced_limit = Functions(a[i][0],omegaAxe).NoBigBang(omegaAxe, 'permission')       # si n'est pas valable, alors sera remplacer par le Modele standard (à voir après)
						#if lamdalist[i] >= Bounced_limit:
						if Bounced_limit==0:
							aa.append(Cov_c)
						else:
							#Mlist_tild = m_app - ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,bbeta,0,0)
							Mlist_tild = MM*1.0
							#ZETA_tild = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,0,0,0)
							
							ten_power_m_M = (10**((m_app - Mlist_tild)/5.))**bet
							Jacob = ApparentMagnitude(z, a[i][0], omegaAxe, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Jacobian(z)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
							WeightingManifold = Jacob * ten_power_m_M
							#WeightingManifold = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, bet)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
							
							WeightingSum = [b/sum(b) for b in WeightingManifold]
							DetermineBeta = [abs((max(b) - min(b))/(max(b) + min(b))) for b in WeightingSum]
							beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
							#L_1_Manifold = [1.0 + (1./np.log(len(b))*sum(b*np.log(b))) for b in WeightingSum]
							#beta = Beta[L_1_Manifold.index(min(L_1_Manifold))][0]
							beta_pliot = beta*1.0

							#w_k = (ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta))/sum((ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta)))
							ten_power_m_M_beta = (10**((m_app - Mlist_tild)/5.))**beta
							w_k = (Jacob*ten_power_m_M_beta)/sum(Jacob*ten_power_m_M_beta)

							COVA_riance = sum(w_k*(Mlist_tild - sum(w_k*Mlist_tild ))*(m_app - sum(w_k*m_app)))

							COR_elation = COVA_riance
							#COR_elation = (1./((len(Mlist_tild)/(len(Mlist_tild) - 1.))*sqrt(sum(w_k*(Mlist_tild - sum(w_k*Mlist_tild ))**2.0))*sqrt(sum(w_k*(m_app - sum(w_k*m_app))**2.0))))*COVA_riance

							aa.append(COR_elation)

					return aa

				
				
				
				#Heaviside = Functions(0.7, 0.3).Heavisidefuction(mass_hst-10)
				#mass_hst = t['3rdvar']	#log10(Mass_of_host_galaxy) it must be > 10
				
##################################################################################################################################################################################################

				for ij in range(len(bbetalist)):
				#for ij in range(1):
				#	bbetalist[ij] = bbeta*1.0
				
					All_models = array(All_models_lambda_omega)*1.0
					
					m_app = m_app_no_correct - bbetalist[ij]*clist
					Cov_All_Models = func(All_models,redshift, Beta)

					#plt.figure(152)
					#plt.plot(bbetalist[ij], Cov_All_Models[0], marker='.', color='r')

					#WeightingNormalisation_All_Models = func(All_models,redshift, Beta)
					#Cov_All_Models = [(1./(sqrt(sum(b[0]*(b[1] - sum(b[0]*b[1]))**2.0))*sqrt(sum(b[0]*(m_app - sum(b[0]*m_app))**2.0))))*sum(b[0]*(b[1] - sum(b[0]*b[1]))*(m_app - sum(b[0]*m_app))) for b in WeightingNormalisation_All_Models]	# (b[0], b[1]) = (w_k, Mlist_tild)
				
					# compute of mean of all true models
					#MeanCov_All_Models = np.mean(Cov_All_Models)
					MeanCov_All_Models = Cov_c
					
					#print "Cov_All_Models = ", Cov_All_Models, "\n"
					print "len(Cov_All_Models) = ", len(Cov_All_Models), "\n"
					print "max(Cov_All_Models) = ", max(Cov_All_Models), "\n"
					print "min(Cov_All_Models) = ", min(Cov_All_Models), "\n"
					print "MeanCov_All_Models = ", MeanCov_All_Models, "\n"
					
					print "len(All_models) = ", len(All_models), "\n"
					print "bbetalist[ij] = ", bbetalist[ij], "\n"
					
					All_models = np.asarray(All_models)
					Cov_All_Models = np.asarray(Cov_All_Models)
					
					Cov_All_Models = Cov_All_Models - MeanCov_All_Models
					
					#***********grid of covariance**********************
					
					
					#===========For the save of the data====================
					OmegaOfGrid = All_models[:,1]*1.0
					LambdaOfGrid = All_models[:,0]*1.0
					CorrCoefficient = Cov_All_Models*1.0
					#=======================================================

					
					#***************************************************************************************

					Covariance_Vert = []
					Covariance_Hor = []
					
					omega0Covariance_Vert = []
					lambda0Covariance_Vert = []
					q0Covariance_Vert = []

					omega0Covariance_Hor = []
					lambda0Covariance_Hor = []
					q0Covariance_Hor = []
					
					# make the array as matrix to search the null correlation curve on the grid!
					Cov_All_Models = Cov_All_Models.reshape(sqrt(len(All_models[:,1])), sqrt(len(All_models[:,0])))
					
					All_modelslambda = All_models[:,0].reshape(sqrt(len(All_models[:,0])), sqrt(len(All_models[:,0]))) # For matrix of lambda
					All_modelsomega = All_models[:,1].reshape(sqrt(len(All_models[:,1])), sqrt(len(All_models[:,0]))) # For matrix of omega

					#============**** To make the contours of confidence levels ****=============  
					#CS = pl.contour(All_modelsomega, All_modelslambda, abs(Cov_All_Models), 3, linewidths=np.arange(.5, 4, .5), colors=('r', 'green', 'blue', (1,1,0), '#afeeee', '0.5'))
					#plt.clabel(CS, inline=1, fontsize=10)
					#============***************************************************=============  
					#plt.show()
					#search of covariance zero by HORIZONTAL interpolation
					omegalist = linspace(np.float64(10**-6), np.float64(1.0), 30) 
					lamdalist = linspace(np.float64(10**-6), np.float64(2.0), 30)
					for i in range(len(omegalist)):
						#----------------------------------------------------
						
						Condition_1 = any(Cov_All_Models[i] == 0.0)
						if Condition_1 == True:
							rr = Cov_All_Models[i]
							bb = All_modelsomega[i] # j'insiste sur omegalist parce que les lignes dans la matrice lambdalist sont toujours divisé avec des valeurs égaux. (interpolation horizontale de la grille (c.a.d de la matrice))# pareil pour omegalist si je fait aussi l'interpolation vertical.
							
							while Condition_1 == True:
								rr = list(rr)
								bb = list(bb)
								
								bb.remove(bb[rr.index(0.0)])
								rr.remove(0.0)
								
								rr = array(rr)
								bb = array(bb)
								
								Condition_1 = any(rr == 0.0)
							
						else:
							rr = Cov_All_Models[i]
							bb = All_modelsomega[i]
							
							
						rr = array(rr)
						bb = array(bb)
						
						Condition_2 = all(sign(rr)> 0)
						Condition_3 = all(sign(rr)< 0)
						
						if Condition_1 == False and Condition_2 == False:
							Cov_0 = 0.0
							lambda0 = All_modelslambda[i][0]# lambda est constant ici / En réalité, toujours lamda va etre le meme dans All_modelslambda[i] parce que depuis le début, j'ai fait les deux boucles en fixant lambda0 et faire tourner la deuxième boucle sur omega0.
							#omega0 = Functions(0.7, 0.3).Nonlinear_Interpolation(bb, rr, Cov_0)
							#===================================
							for bi in range(1, len(rr)-2):
								if sign(rr[bi])!= sign(rr[bi-1]):
									omega0 = Functions(0.7, 0.3).Nonlinear_Interpolation_one_by_one(bb[bi-1:bi+2], rr[bi-1:bi+2], Cov_0)
									q0 = omega0/2. - lambda0
									Covariance_Hor.append(Cov_0)
									omega0Covariance_Hor.append(omega0)
									lambda0Covariance_Hor.append(lambda0)
									q0Covariance_Hor.append(q0)
								else:
									pass
							#===================================	      
							"""
							q0 = omega0/2. - lambda0
							
							Covariance_Hor.append(Cov_0)
							omega0Covariance_Hor.append(omega0)
							lambda0Covariance_Hor.append(lambda0)
							q0Covariance_Hor.append(q0)
							"""
						else:
							pass
							
						#------------------------------------------------------

					#search of covariance zero by VERTICAL interpolation
					for i in range(len(lamdalist)):
						
						#----------------------------------------------------
						
						Condition_1 = any(Cov_All_Models[:,i] == 0.0)
						if Condition_1 == True:
							rr = Cov_All_Models[:,i]
							bb = All_modelslambda[:,i] # j'insiste sur omegalist parce que les lignes dans la matrice lambdalist sont toujours divisé avec des valeurs égaux. (interpolation horizontale de la grille (c.a.d de la matrice))# pareil pour omegalist si je fait aussi l'interpolation vertical.
							
							while Condition_1 == True:
								rr = list(rr)
								bb = list(bb)
								
								bb.remove(bb[rr.index(0.0)])
								rr.remove(0.0)
								
								rr = array(rr)
								bb = array(bb)
								
								Condition_1 = any(rr == 0.0)
							
						else:
							rr = Cov_All_Models[:,i]
							bb = All_modelslambda[:,i]
							
							
						rr = array(rr)
						bb = array(bb)
						
						Condition_2 = all(sign(rr)> 0)
						Condition_3 = all(sign(rr)< 0)
						
						if Condition_1 == False and Condition_2 == False:
							Cov_0 = 0.0
							#lambda0 = Functions(0.7, 0.3).Nonlinear_Interpolation(bb, rr, Cov_0) # lambda est constant ici / En réalité, toujours lamda va etre le meme dans All_modelslambda[i] parce que depuis le début, j'ai fait les deux boucles en fixant lambda0 et faire tourner la deuxième boucle sur omega0.
							omega0 = All_modelsomega[:,i][0]
							#===================================
							for bi in range(1, len(rr)-2):
								if sign(rr[bi])!= sign(rr[bi-1]):
									lambda0 = Functions(0.7, 0.3).Nonlinear_Interpolation_one_by_one(bb[bi-1:bi+2], rr[bi-1:bi+2], Cov_0)
									q0 = omega0/2. - lambda0
									Covariance_Vert.append(Cov_0)
									omega0Covariance_Vert.append(omega0)
									lambda0Covariance_Vert.append(lambda0)
									q0Covariance_Vert.append(q0)
								
								else:
									pass
							#===================================	 
							"""
							q0 = omega0/2. - lambda0
							
							Covariance_Vert.append(Cov_0)
							omega0Covariance_Vert.append(omega0)
							lambda0Covariance_Vert.append(lambda0)
							q0Covariance_Vert.append(q0)
							"""
						else:
							pass
							
						#------------------------------------------------------
				
					
					Y_Fit, X_Fit = [], []
					lambda0Covariance = lambda0Covariance_Hor + lambda0Covariance_Vert
					omega0Covariance = omega0Covariance_Hor + omega0Covariance_Vert
					#"""
					lambda0Covariance1 = list(array(lambda0Covariance)*1.0)
					omega0Covariance1 = list(array(omega0Covariance)*1.0)
					
					if len(lambda0Covariance1) > 3 and len(omega0Covariance1) > 3:
						print "lambda0Covariance1 = ", lambda0Covariance1
						print "omega0Covariance1 = ", omega0Covariance1
						
						for i in range(len(lambda0Covariance)):	
							X_Fit.append(min(omega0Covariance1))
							Y_Fit.append(lambda0Covariance1[omega0Covariance1.index(min(omega0Covariance1))])
							
							lambda0Covariance1.remove(lambda0Covariance1[omega0Covariance1.index(min(omega0Covariance1))])
							omega0Covariance1.remove(min(omega0Covariance1))
						
						Y_Fit2, X_Fit2 = list(array(lambda0Covariance)*1.0), list(array(omega0Covariance)*1.0)
						X_Fit3, Y_Fit3 = [], []
						
						for i in range(len(X_Fit2)):
							if X_Fit2[i] == min(X_Fit2) or Y_Fit2[i] == min(Y_Fit2):
								pass
							elif X_Fit2[i] == max(X_Fit2) or Y_Fit2[i] == max(Y_Fit2):
								pass
							      
							else:
								X_Fit3.append(X_Fit2[i])
								Y_Fit3.append(Y_Fit2[i])
								
						X, Y = [], []
						for i in range(len(X_Fit3)):
							Y.append(min(Y_Fit3))
							X.append(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
							X_Fit3.remove(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
							Y_Fit3.remove(min(Y_Fit3))

						print "len(X) = ", len(X)
						print "len(Y) = ", len(Y)
						
						print "X = ", X
						print "Y = ", Y
						
						X_omega = X_omega + X
						Y_lambda = Y_lambda + Y
						Z_bbeta = Z_bbeta + [bbetalist[ij]]*len(X)
						
						plt.figure(125)
						plt.plot(X, Y)
						
						Lamdalimit=[]
						omegalistlimit = linspace(min(omegalist), max(omegalist), 200)
						
						for i in range(len(omegalistlimit)):
							Lamdalimit.append(Functions(0.7, 0.3).NoBigBang(omegalistlimit[i], 'limit'))

						plt.fill_between(omegalistlimit, Lamdalimit, max(Lamdalimit), color=(0.5,0.5,0.5))
						plt.text(0.001, 1.468, 'No Big Bang', color='w')
						plt.text(0.7, 1.32, '($\\Omega_{\\circ}$, $\\lambda_{\\circ}$) = ('+str(omega0_Model)+', '+str(lambda0_Model)+')', color='r')
						
						plt.plot([omega0_Model],[lambda0_Model], marker='o', color='r')
						plt.xlabel('$\\Omega_{\\circ}$', fontsize=16)
						plt.ylabel('$\\lambda_{\\circ}$', fontsize=16)
						grid(True)


						
					else:
						pass
      
##################################################################################################################################################################################################

##################################################################################################################################################################################################

				for ij in range(len(lamdalist)):

					All_models = array(All_models_omega_bbeta)*1.0
					
					Cov_All_Models = func1(lamdalist[ij], All_models,redshift, Beta)

					MeanCov_All_Models = Cov_c
					
					#print "Cov_All_Models = ", Cov_All_Models, "\n"
					print "len(Cov_All_Models) = ", len(Cov_All_Models), "\n"
					print "max(Cov_All_Models) = ", max(Cov_All_Models), "\n"
					print "min(Cov_All_Models) = ", min(Cov_All_Models), "\n"
					print "MeanCov_All_Models = ", MeanCov_All_Models, "\n"
					
					print "len(All_models) = ", len(All_models), "\n"
					print "bbetalist[ij] = ", lamdalist[ij], "\n"
					
					All_models = np.asarray(All_models)
					Cov_All_Models = np.asarray(Cov_All_Models)
					
					Cov_All_Models = Cov_All_Models - MeanCov_All_Models
					
					#***********grid of covariance**********************
					
					
					#===========For the save of the data====================
					OmegaOfGrid = All_models[:,1]*1.0
					LambdaOfGrid = All_models[:,0]*1.0
					CorrCoefficient = Cov_All_Models*1.0
					#=======================================================

					
					#***************************************************************************************

					Covariance_Vert = []
					Covariance_Hor = []
					
					omega0Covariance_Vert = []
					lambda0Covariance_Vert = []
					q0Covariance_Vert = []

					omega0Covariance_Hor = []
					lambda0Covariance_Hor = []
					q0Covariance_Hor = []
					
					# make the array as matrix to search the null correlation curve on the grid!
					Cov_All_Models = Cov_All_Models.reshape(sqrt(len(All_models[:,1])), sqrt(len(All_models[:,0])))
					
					All_modelslambda = All_models[:,0].reshape(sqrt(len(All_models[:,0])), sqrt(len(All_models[:,0]))) # For matrix of lambda
					All_modelsomega = All_models[:,1].reshape(sqrt(len(All_models[:,1])), sqrt(len(All_models[:,0]))) # For matrix of omega

					#============**** To make the contours of confidence levels ****=============  
					#CS = pl.contour(All_modelsomega, All_modelslambda, abs(Cov_All_Models), 3, linewidths=np.arange(.5, 4, .5), colors=('r', 'green', 'blue', (1,1,0), '#afeeee', '0.5'))
					#plt.clabel(CS, inline=1, fontsize=10)
					#============***************************************************=============  
					#plt.show()
					#search of covariance zero by HORIZONTAL interpolation
					omegalist = linspace(np.float64(10**-6), np.float64(1.0), 30) 
					lamdalist = linspace(np.float64(10**-6), np.float64(2.0), 30)
					for i in range(len(omegalist)):
						#----------------------------------------------------
						
						Condition_1 = any(Cov_All_Models[i] == 0.0)
						if Condition_1 == True:
							rr = Cov_All_Models[i]
							bb = All_modelsomega[i] # j'insiste sur omegalist parce que les lignes dans la matrice lambdalist sont toujours divisé avec des valeurs égaux. (interpolation horizontale de la grille (c.a.d de la matrice))# pareil pour omegalist si je fait aussi l'interpolation vertical.
							
							while Condition_1 == True:
								rr = list(rr)
								bb = list(bb)
								
								bb.remove(bb[rr.index(0.0)])
								rr.remove(0.0)
								
								rr = array(rr)
								bb = array(bb)
								
								Condition_1 = any(rr == 0.0)
							
						else:
							rr = Cov_All_Models[i]
							bb = All_modelsomega[i]
							
							
						rr = array(rr)
						bb = array(bb)
						
						Condition_2 = all(sign(rr)> 0)
						Condition_3 = all(sign(rr)< 0)
						
						if Condition_1 == False and Condition_2 == False:
							Cov_0 = 0.0
							lambda0 = All_modelslambda[i][0]# lambda est constant ici / En réalité, toujours lamda va etre le meme dans All_modelslambda[i] parce que depuis le début, j'ai fait les deux boucles en fixant lambda0 et faire tourner la deuxième boucle sur omega0.
							#omega0 = Functions(0.7, 0.3).Nonlinear_Interpolation(bb, rr, Cov_0)
							#===================================
							for bi in range(1, len(rr)-2):
								if sign(rr[bi])!= sign(rr[bi-1]):
									omega0 = Functions(0.7, 0.3).Nonlinear_Interpolation_one_by_one(bb[bi-1:bi+2], rr[bi-1:bi+2], Cov_0)
									q0 = omega0/2. - lambda0
									Covariance_Hor.append(Cov_0)
									omega0Covariance_Hor.append(omega0)
									lambda0Covariance_Hor.append(lambda0)
									q0Covariance_Hor.append(q0)
								else:
									pass
							#===================================	      
							"""
							q0 = omega0/2. - lambda0
							
							Covariance_Hor.append(Cov_0)
							omega0Covariance_Hor.append(omega0)
							lambda0Covariance_Hor.append(lambda0)
							q0Covariance_Hor.append(q0)
							"""
						else:
							pass
							
						#------------------------------------------------------

					#search of covariance zero by VERTICAL interpolation
					for i in range(len(lamdalist)):
						
						#----------------------------------------------------
						
						Condition_1 = any(Cov_All_Models[:,i] == 0.0)
						if Condition_1 == True:
							rr = Cov_All_Models[:,i]
							bb = All_modelslambda[:,i] # j'insiste sur omegalist parce que les lignes dans la matrice lambdalist sont toujours divisé avec des valeurs égaux. (interpolation horizontale de la grille (c.a.d de la matrice))# pareil pour omegalist si je fait aussi l'interpolation vertical.
							
							while Condition_1 == True:
								rr = list(rr)
								bb = list(bb)
								
								bb.remove(bb[rr.index(0.0)])
								rr.remove(0.0)
								
								rr = array(rr)
								bb = array(bb)
								
								Condition_1 = any(rr == 0.0)
							
						else:
							rr = Cov_All_Models[:,i]
							bb = All_modelslambda[:,i]
							
							
						rr = array(rr)
						bb = array(bb)
						
						Condition_2 = all(sign(rr)> 0)
						Condition_3 = all(sign(rr)< 0)
						
						if Condition_1 == False and Condition_2 == False:
							Cov_0 = 0.0
							#lambda0 = Functions(0.7, 0.3).Nonlinear_Interpolation(bb, rr, Cov_0) # lambda est constant ici / En réalité, toujours lamda va etre le meme dans All_modelslambda[i] parce que depuis le début, j'ai fait les deux boucles en fixant lambda0 et faire tourner la deuxième boucle sur omega0.
							omega0 = All_modelsomega[:,i][0]
							#===================================
							for bi in range(1, len(rr)-2):
								if sign(rr[bi])!= sign(rr[bi-1]):
									lambda0 = Functions(0.7, 0.3).Nonlinear_Interpolation_one_by_one(bb[bi-1:bi+2], rr[bi-1:bi+2], Cov_0)
									q0 = omega0/2. - lambda0
									Covariance_Vert.append(Cov_0)
									omega0Covariance_Vert.append(omega0)
									lambda0Covariance_Vert.append(lambda0)
									q0Covariance_Vert.append(q0)
								
								else:
									pass
							#===================================	 
							"""
							q0 = omega0/2. - lambda0
							
							Covariance_Vert.append(Cov_0)
							omega0Covariance_Vert.append(omega0)
							lambda0Covariance_Vert.append(lambda0)
							q0Covariance_Vert.append(q0)
							"""
						else:
							pass
							
						#------------------------------------------------------
				
					
					Y_Fit, X_Fit = [], []
					lambda0Covariance = lambda0Covariance_Hor + lambda0Covariance_Vert
					omega0Covariance = omega0Covariance_Hor + omega0Covariance_Vert
					#"""
					lambda0Covariance1 = list(array(lambda0Covariance)*1.0)
					omega0Covariance1 = list(array(omega0Covariance)*1.0)
					
					if len(lambda0Covariance1) > 3 and len(omega0Covariance1) > 3:
						print "lambda0Covariance1 = ", lambda0Covariance1
						print "omega0Covariance1 = ", omega0Covariance1
						
						for i in range(len(lambda0Covariance)):	
							X_Fit.append(min(omega0Covariance1))
							Y_Fit.append(lambda0Covariance1[omega0Covariance1.index(min(omega0Covariance1))])
							
							lambda0Covariance1.remove(lambda0Covariance1[omega0Covariance1.index(min(omega0Covariance1))])
							omega0Covariance1.remove(min(omega0Covariance1))
						
						Y_Fit2, X_Fit2 = list(array(lambda0Covariance)*1.0), list(array(omega0Covariance)*1.0)
						X_Fit3, Y_Fit3 = [], []
						
						for i in range(len(X_Fit2)):
							if X_Fit2[i] == min(X_Fit2) or Y_Fit2[i] == min(Y_Fit2):
								pass
							elif X_Fit2[i] == max(X_Fit2) or Y_Fit2[i] == max(Y_Fit2):
								pass
							      
							else:
								X_Fit3.append(X_Fit2[i])
								Y_Fit3.append(Y_Fit2[i])
								
						X, Y = [], []
						for i in range(len(X_Fit3)):
							Y.append(min(Y_Fit3))
							X.append(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
							X_Fit3.remove(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
							Y_Fit3.remove(min(Y_Fit3))

						print "len(X) = ", len(X)
						print "len(Y) = ", len(Y)
						
						print "X = ", X
						print "Y = ", Y
						X_omega = X_omega + X
						Y_lambda = Y_lambda + [lamdalist[ij]]*len(X)
						Z_bbeta = Z_bbeta + Y
						
						plt.figure(126)
						plt.plot(X, Y)

						plt.plot([omega0_Model],[bbeta], marker='o', color='r')
						plt.xlabel('$\\Omega_{\\circ}$', fontsize=16)
						plt.ylabel('$\\beta$', fontsize=16)
						grid(True)


						
					else:
						pass
      
##################################################################################################################################################################################################

##################################################################################################################################################################################################

				for ij in range(len(omegalist)):

					All_models = array(All_models_bbeta_lambda)*1.0
					
					Cov_All_Models = func2(omegalist[ij], All_models,redshift, Beta)

					MeanCov_All_Models = Cov_c
					
					#print "Cov_All_Models = ", Cov_All_Models, "\n"
					print "len(Cov_All_Models) = ", len(Cov_All_Models), "\n"
					print "max(Cov_All_Models) = ", max(Cov_All_Models), "\n"
					print "min(Cov_All_Models) = ", min(Cov_All_Models), "\n"
					print "MeanCov_All_Models = ", MeanCov_All_Models, "\n"
					
					print "len(All_models) = ", len(All_models), "\n"
					print "bbetalist[ij] = ", omegalist[ij], "\n"
					
					All_models = np.asarray(All_models)
					Cov_All_Models = np.asarray(Cov_All_Models)
					
					Cov_All_Models = Cov_All_Models - MeanCov_All_Models
					
					#***********grid of covariance**********************
					
					
					#===========For the save of the data====================
					OmegaOfGrid = All_models[:,1]*1.0
					LambdaOfGrid = All_models[:,0]*1.0
					CorrCoefficient = Cov_All_Models*1.0
					#=======================================================

					
					#***************************************************************************************

					Covariance_Vert = []
					Covariance_Hor = []
					
					omega0Covariance_Vert = []
					lambda0Covariance_Vert = []
					q0Covariance_Vert = []

					omega0Covariance_Hor = []
					lambda0Covariance_Hor = []
					q0Covariance_Hor = []
					
					# make the array as matrix to search the null correlation curve on the grid!
					Cov_All_Models = Cov_All_Models.reshape(sqrt(len(All_models[:,1])), sqrt(len(All_models[:,0])))
					
					All_modelslambda = All_models[:,0].reshape(sqrt(len(All_models[:,0])), sqrt(len(All_models[:,0]))) # For matrix of lambda
					All_modelsomega = All_models[:,1].reshape(sqrt(len(All_models[:,1])), sqrt(len(All_models[:,0]))) # For matrix of omega

					#============**** To make the contours of confidence levels ****=============  
					#CS = pl.contour(All_modelsomega, All_modelslambda, abs(Cov_All_Models), 3, linewidths=np.arange(.5, 4, .5), colors=('r', 'green', 'blue', (1,1,0), '#afeeee', '0.5'))
					#plt.clabel(CS, inline=1, fontsize=10)
					#============***************************************************=============  
					#plt.show()
					#search of covariance zero by HORIZONTAL interpolation
					omegalist = linspace(np.float64(10**-6), np.float64(1.0), 30) 
					lamdalist = linspace(np.float64(10**-6), np.float64(2.0), 30)
					for i in range(len(omegalist)):
						#----------------------------------------------------
						
						Condition_1 = any(Cov_All_Models[i] == 0.0)
						if Condition_1 == True:
							rr = Cov_All_Models[i]
							bb = All_modelsomega[i] # j'insiste sur omegalist parce que les lignes dans la matrice lambdalist sont toujours divisé avec des valeurs égaux. (interpolation horizontale de la grille (c.a.d de la matrice))# pareil pour omegalist si je fait aussi l'interpolation vertical.
							
							while Condition_1 == True:
								rr = list(rr)
								bb = list(bb)
								
								bb.remove(bb[rr.index(0.0)])
								rr.remove(0.0)
								
								rr = array(rr)
								bb = array(bb)
								
								Condition_1 = any(rr == 0.0)
							
						else:
							rr = Cov_All_Models[i]
							bb = All_modelsomega[i]
							
							
						rr = array(rr)
						bb = array(bb)
						
						Condition_2 = all(sign(rr)> 0)
						Condition_3 = all(sign(rr)< 0)
						
						if Condition_1 == False and Condition_2 == False:
							Cov_0 = 0.0
							lambda0 = All_modelslambda[i][0]# lambda est constant ici / En réalité, toujours lamda va etre le meme dans All_modelslambda[i] parce que depuis le début, j'ai fait les deux boucles en fixant lambda0 et faire tourner la deuxième boucle sur omega0.
							#omega0 = Functions(0.7, 0.3).Nonlinear_Interpolation(bb, rr, Cov_0)
							#===================================
							for bi in range(1, len(rr)-2):
								if sign(rr[bi])!= sign(rr[bi-1]):
									omega0 = Functions(0.7, 0.3).Nonlinear_Interpolation_one_by_one(bb[bi-1:bi+2], rr[bi-1:bi+2], Cov_0)
									q0 = omega0/2. - lambda0
									Covariance_Hor.append(Cov_0)
									omega0Covariance_Hor.append(omega0)
									lambda0Covariance_Hor.append(lambda0)
									q0Covariance_Hor.append(q0)
								else:
									pass
							#===================================	      
							"""
							q0 = omega0/2. - lambda0
							
							Covariance_Hor.append(Cov_0)
							omega0Covariance_Hor.append(omega0)
							lambda0Covariance_Hor.append(lambda0)
							q0Covariance_Hor.append(q0)
							"""
						else:
							pass
							
						#------------------------------------------------------

					#search of covariance zero by VERTICAL interpolation
					for i in range(len(lamdalist)):
						
						#----------------------------------------------------
						
						Condition_1 = any(Cov_All_Models[:,i] == 0.0)
						if Condition_1 == True:
							rr = Cov_All_Models[:,i]
							bb = All_modelslambda[:,i] # j'insiste sur omegalist parce que les lignes dans la matrice lambdalist sont toujours divisé avec des valeurs égaux. (interpolation horizontale de la grille (c.a.d de la matrice))# pareil pour omegalist si je fait aussi l'interpolation vertical.
							
							while Condition_1 == True:
								rr = list(rr)
								bb = list(bb)
								
								bb.remove(bb[rr.index(0.0)])
								rr.remove(0.0)
								
								rr = array(rr)
								bb = array(bb)
								
								Condition_1 = any(rr == 0.0)
							
						else:
							rr = Cov_All_Models[:,i]
							bb = All_modelslambda[:,i]
							
							
						rr = array(rr)
						bb = array(bb)
						
						Condition_2 = all(sign(rr)> 0)
						Condition_3 = all(sign(rr)< 0)
						
						if Condition_1 == False and Condition_2 == False:
							Cov_0 = 0.0
							#lambda0 = Functions(0.7, 0.3).Nonlinear_Interpolation(bb, rr, Cov_0) # lambda est constant ici / En réalité, toujours lamda va etre le meme dans All_modelslambda[i] parce que depuis le début, j'ai fait les deux boucles en fixant lambda0 et faire tourner la deuxième boucle sur omega0.
							omega0 = All_modelsomega[:,i][0]
							#===================================
							for bi in range(1, len(rr)-2):
								if sign(rr[bi])!= sign(rr[bi-1]):
									lambda0 = Functions(0.7, 0.3).Nonlinear_Interpolation_one_by_one(bb[bi-1:bi+2], rr[bi-1:bi+2], Cov_0)
									q0 = omega0/2. - lambda0
									Covariance_Vert.append(Cov_0)
									omega0Covariance_Vert.append(omega0)
									lambda0Covariance_Vert.append(lambda0)
									q0Covariance_Vert.append(q0)
								
								else:
									pass
							#===================================	 
							"""
							q0 = omega0/2. - lambda0
							
							Covariance_Vert.append(Cov_0)
							omega0Covariance_Vert.append(omega0)
							lambda0Covariance_Vert.append(lambda0)
							q0Covariance_Vert.append(q0)
							"""
						else:
							pass
							
						#------------------------------------------------------
				
					
					Y_Fit, X_Fit = [], []
					lambda0Covariance = lambda0Covariance_Hor + lambda0Covariance_Vert
					omega0Covariance = omega0Covariance_Hor + omega0Covariance_Vert
					#"""
					lambda0Covariance1 = list(array(lambda0Covariance)*1.0)
					omega0Covariance1 = list(array(omega0Covariance)*1.0)
					
					if len(lambda0Covariance1) > 3 and len(omega0Covariance1) > 3:
						print "lambda0Covariance1 = ", lambda0Covariance1
						print "omega0Covariance1 = ", omega0Covariance1
						
						for i in range(len(lambda0Covariance)):	
							X_Fit.append(min(omega0Covariance1))
							Y_Fit.append(lambda0Covariance1[omega0Covariance1.index(min(omega0Covariance1))])
							
							lambda0Covariance1.remove(lambda0Covariance1[omega0Covariance1.index(min(omega0Covariance1))])
							omega0Covariance1.remove(min(omega0Covariance1))
						
						Y_Fit2, X_Fit2 = list(array(lambda0Covariance)*1.0), list(array(omega0Covariance)*1.0)
						X_Fit3, Y_Fit3 = [], []
						
						for i in range(len(X_Fit2)):
							if X_Fit2[i] == min(X_Fit2) or Y_Fit2[i] == min(Y_Fit2):
								pass
							elif X_Fit2[i] == max(X_Fit2) or Y_Fit2[i] == max(Y_Fit2):
								pass
							      
							else:
								X_Fit3.append(X_Fit2[i])
								Y_Fit3.append(Y_Fit2[i])
								
						X, Y = [], []
						for i in range(len(X_Fit3)):
							Y.append(min(Y_Fit3))
							X.append(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
							X_Fit3.remove(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
							Y_Fit3.remove(min(Y_Fit3))

						print "len(X) = ", len(X)
						print "len(Y) = ", len(Y)
						
						print "X = ", X
						print "Y = ", Y
						X_omega = X_omega + [omegalist[ij]]*len(X)
						Y_lambda = Y_lambda + Y
						Z_bbeta = Z_bbeta + X
						
						plt.figure(127)
						plt.plot(X, Y)
						
						plt.plot([bbeta],[lambda0_Model], marker='o', color='r')
						plt.xlabel('$\\beta$', fontsize=16)
						plt.ylabel('$\\lambda_{\\circ}$', fontsize=16)
						grid(True)


						
					else:
						pass
      
##################################################################################################################################################################################################


				fig = plt.figure()
				ax = fig.gca(projection='3d')
				ax.plot_trisurf(X_omega, Y_lambda, Z_bbeta, linewidth=0.2, antialiased=True)
				
				ax.scatter(omega0_Model, lambda0_Model, bbeta, color='m', marker='o')
				ax.scatter(omega0_Model, lambda0_Model, bbeta, color='m', marker='o')
				ax.scatter(omega0_Model, lambda0_Model, bbeta, color='m', marker='o')
				ax.scatter(omega0_Model, lambda0_Model, bbeta, color='m', marker='o')
				
				ax.set_xlabel('$\\Omega_{\\circ}$', fontsize=14)
				ax.set_ylabel('$\\lambda_{\\circ}$', fontsize=14)
				ax.set_zlabel('$\\beta$', fontsize=14)
				plt.show()
				
				if KS == 'K-S_OK':
				  
					DD = []
					
				elif KS == 'K-S_NO':
				    
					DD = []
				      
				else:
					print "You must to determine the condition of the Kolmogrov-Smirnov test as 'K-S_OK' or 'K-S_NO'. ", "\n"

				return OmegaOfGrid, LambdaOfGrid, CorrCoefficient, X, Y, DD


			else:
			      print "You must to determine the condition of the ### Variables ### as 'stretch_and_color' or 'stretch_only' or 'color_only'. ", "\n"
			
		elif Scanning_on_parameters == 'alfa_bbeta':
			if Variables == 'stretch_and_color':
				pass
			elif Variables == 'stretch_only':
				pass
				#plt.show()
				
				if KS == 'K-S_OK':
				  
					DD = []
					
				elif KS == 'K-S_NO':
				    
					DD = []
				      
				else:
					print "You must to determine the condition of the Kolmogrov-Smirnov test as 'K-S_OK' or 'K-S_NO'. ", "\n"

				return BetaOfGrid, AlfaOfGrid, CorrCoefficient, X, Y, DD


			elif Variables == 'color_only':
				pass
			      
				if KS == 'K-S_OK':
				  
					DD = []
					
				elif KS == 'K-S_NO':
				    
					DD = []
				      
				else:
					print "You must to determine the condition of the Kolmogrov-Smirnov test as 'K-S_OK' or 'K-S_NO'. ", "\n"

				return BetaOfGrid, AlfaOfGrid, CorrCoefficient, X, Y, DD


			else:
			      print "You must to determine the condition of the ### Variables ### as 'stretch_and_color' or 'stretch_only' or 'color_only'. ", "\n"	
			
		else:
			print "You must to determine the condition of the ### Scanning_on_parameters ### as 'lambda0_omega0' or 'alfa_bbeta'. ", "\n"
	  
		
