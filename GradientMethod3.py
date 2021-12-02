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
		bbetalist = linspace(np.float64(1.5), np.float64(4.5), 30) 
		M_0list = linspace(np.float64(-21), np.float64(-18), 20) 
		
		#Betalist = linspace(np.float64(2.1), np.float64(3.1), 100) #100
		Betalist = linspace(np.float64(0.01), np.float64(3.9), 150)
		#Betalist = linspace(np.float64(0.01), np.float64(3.9), 100)
		
		Beta = array([array([a]) for a in Betalist])
		
		if Scanning_on_parameters == 'lambda0_omega0':
			if Variables == 'stretch_and_color':
				#======================================### Firstly for color ###======================================
				All_models = []
				for i in range(len(lamdalist)):
				    for j in range(len(omegalist)):
					    All_models.append((lamdalist[i], omegalist[j]))

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
							
							#cc = ( m_app - M_0 + alfa*(slist -1) - ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,0,0,0) )/bbeta
							#cc = clist*1.0

							#ten_power_ZETA_tild_ZETA_virt = 10**(bet*(ZETA_virt + ZETA_tild)/2.)
							#ten_power_ZETA_virt = (10**(ZETA_virt/5.))**bet

							#ten_power_c_m = (10**((m_app - bbeta*clist)/5.))**bet
							#ten_power_c_m = (10**((m_app - clist)/5.))**bet
							#ten_power_c_m = (10**((-alfa*slist)/5.))**bet
							
							
							WeightingManifold = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, bet)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
							
							#WeightingManifold = WeightingManifold*ten_power_c_m
							
							WeightingSum = [b/sum(b) for b in WeightingManifold]
							#DetermineBeta = [abs(max(b) - min(b)) for b in WeightingSum]
							DetermineBeta = [abs((max(b) - min(b))/(max(b) + min(b))) for b in WeightingSum]
							beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
							#L_1_Manifold = [1.0 + (1./np.log(len(b))*sum(b*np.log(b))) for b in WeightingSum]
							#beta = Beta[L_1_Manifold.index(min(L_1_Manifold))][0]

							#ten_power_c_m = (10**((m_app - clist)/5.))**beta
							#ten_power_c_m = (10**((-alfa*slist)/5.))**beta
							w_k = (ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta))/sum((ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta)))

							#COVA_riance = sum(w_k*(cc - sum(w_k*cc ))*(m_app - sum(w_k*m_app)))
							#COR_elation = (1./((len(bbeta*clist)/(len(bbeta*clist) - 1.))*sqrt(sum(w_k*(bbeta*clist - sum(w_k*bbeta*clist))**2.0))*sqrt(sum(w_k*(m_app - sum(w_k*m_app))**2.0))))*COVA_riance
							
							COVA_riance = sum(w_k*(Mlist_tild - sum(w_k*Mlist_tild ))*(m_app - sum(w_k*m_app)))

							COR_elation = COVA_riance
							#COR_elation = (1./((len(Mlist_tild)/(len(Mlist_tild) - 1.))*sqrt(sum(w_k*(Mlist_tild - sum(w_k*Mlist_tild ))**2.0))*sqrt(sum(w_k*(m_app - sum(w_k*m_app))**2.0))))*COVA_riance

							aa.append(COR_elation)

							#aa.append((w_k, Mlist_tild))
						
							#~plt.figure(44)
							#z2 = z*1.0
							#z2.sort()
							#ww = ApparentMagnitude(z2, a[i][0],a[i][1], 0.0).Weightingfactor(z2, beta, 0.0)
							#wk = ww/sum(ww)
							#plt.plot(z2, len(z2)*wk)
							#~for hg in range(len(w_k)):
							#~	plt.plot(z[hg], w_k[hg], marker='+')
								
							"""
							plt.figure(5)
							L_1 = 1.0 + (1./np.log(len(wk))*sum(wk*np.log(wk)))
							plt.plot(a[i][0], L_1, marker='+', color='g')
							plt.plot(a[i][1], L_1, marker='+', color='r')
							
							plt.figure(6)
							V = Model(a[i][0], a[i][1]).curvature(z2, a[i][0],a[i][1], 0.0)[1]
							plt.plot(z2, V)
							"""
					return aa


				
				Cov_All_Models = func(All_models,redshift, Beta)
				
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
				
				All_models = np.asarray(All_models)
				Cov_All_Models = np.asarray(Cov_All_Models)
				
				Cov_All_Models = Cov_All_Models - MeanCov_All_Models
				
				#***********grid of covariance**********************
				# figure of 3D(surface)
				
				#All_models[:,0] is a liste of lamda
				#All_models[:,1] is a liste of omega
				
				fig = plt.figure(7)
				ax = fig.add_subplot(111)

				col = plt.scatter(All_models[:,1], All_models[:,0], marker='.', s=150, c=Cov_All_Models, linewidths=0.3, cmap=plt.cm.Spectral_r)
				#for vv, ww, dd in zip(All_models[:,1], All_models[:,0], Cov_All_Models):
				#	pl.text(vv, ww, '%.3f' % dd, ha='center', va='bottom')
				

				# Add a colourbar.
				#cax = fig.colorbar(col, orientation='vertical',format='%.30f')
				cax = fig.colorbar(col, orientation='vertical',format='%.4f')
				cax.set_label('Covariance')
				
				plt.plot([omega0_Model],[lambda0_Model], marker='o', color='m')
				plt.xlabel('$\\Omega_{0}$', fontsize=14)
				plt.ylabel('$\\lambda_{0}$', fontsize=14)
				
				#===========For the save of the data====================
				OmegaOfGrid = All_models[:,1]*1.0
				LambdaOfGrid = All_models[:,0]*1.0
				CorrCoefficient = Cov_All_Models*1.0
				#=======================================================

				
				#plt.figure(8)
				#for i in range(len(Cov_All_Models)):
				#	plt.plot(Cov_All_Models[i], All_models[i][0], marker='o', color='r')
				
				#n, bins, patches = plt.hist(Cov_All_Models,40,normed='True')
				#(mu_hist, sigma_hist) = norm.fit(array(Cov_All_Models))
				#y = mlab.normpdf( bins, mu_hist, sigma_hist)
				#plt.plot(bins, y, 'r--', linewidth=2)
				#plt.title(r'$\rho_{0}=%.7f,\ \sigma_{0}=%.3f$' %(mu_hist, sigma_hist), color='r')
				#plt.xlabel('Correlation coefficient', fontsize=14)
				
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
					
				
				#plt.figure(9)
				
				Y_Fit, X_Fit = [], []
				lambda0Covariance = lambda0Covariance_Hor + lambda0Covariance_Vert
				omega0Covariance = omega0Covariance_Hor + omega0Covariance_Vert
				#"""
				lambda0Covariance1 = list(array(lambda0Covariance)*1.0)
				omega0Covariance1 = list(array(omega0Covariance)*1.0)
				
				print "lambda0Covariance1 = ", lambda0Covariance1
				print "omega0Covariance1 = ", omega0Covariance1
				
				for i in range(len(lambda0Covariance_Vert)):
					plt.plot(omega0Covariance_Vert[i], lambda0Covariance_Vert[i], marker='.', color='r')
					
				for i in range(len(lambda0Covariance_Hor)):
					plt.plot(omega0Covariance_Hor[i], lambda0Covariance_Hor[i], marker='.', color='g')
					
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
				if omega0_Model < min(X):
					omegalist2 = linspace(omega0_Model-0.02, max(X), 200)
					
				elif omega0_Model > max(X):
					omegalist2 = linspace(min(X), omega0_Model+0.02, 200)
						
				else:
					omegalist2 = linspace(min(X), max(X), 200)
						
						
				plt.plot(X, Y, color=couleur)
				
				Lamdalimit=[]
				omegalist = linspace(min(omegalist), max(omegalist), 200)
				
				for i in range(len(omegalist)):
					Lamdalimit.append(Functions(0.7, 0.3).NoBigBang(omegalist[i], 'limit'))

				plt.fill_between(omegalist, Lamdalimit, max(Lamdalimit), color=(0.5,0.5,0.5))
				plt.text(0.001, 1.468, 'No Big Bang', color='w')
				plt.text(0.7, 1.32, '($\\Omega_{0}$, $\\lambda_{0}$) = ('+str(omega0_Model)+', '+str(lambda0_Model)+')', color='r')
				
				plt.plot([omega0_Model],[lambda0_Model], marker='o', color='r')
				plt.xlabel('$\\Omega_{0}$', fontsize=16)
				plt.ylabel('$\\lambda_{0}$', fontsize=16)
				plt.title('The null correlation curve in the ($\\Omega_{0}$, $\\lambda_{0}$) plan', color='b')
				grid(True)
						
				plt.xlim(min(omegalist), max(omegalist))                                                                                                                                                                                                                                                                                                                                                                                                               
				plt.ylim(min(lamdalist), max(lamdalist))      

				
				plt.figure(9)
				
				Y_Fit, X_Fit = [], []
				lambda0Covariance = lambda0Covariance_Hor + lambda0Covariance_Vert
				omega0Covariance = omega0Covariance_Hor + omega0Covariance_Vert
				lambda0Covariance.append(lambda0_Model)
				omega0Covariance.append(omega0_Model)
				#"""
				lambda0Covariance1 = list(array(lambda0Covariance)*1.0)
				omega0Covariance1 = list(array(omega0Covariance)*1.0)
				
				print "lambda0Covariance1 = ", lambda0Covariance1
				print "omega0Covariance1 = ", omega0Covariance1
				
				for i in range(len(lambda0Covariance_Vert)):
					plt.plot(omega0Covariance_Vert[i], lambda0Covariance_Vert[i], marker='.', color='r')
					
				for i in range(len(lambda0Covariance_Hor)):
					plt.plot(omega0Covariance_Hor[i], lambda0Covariance_Hor[i], marker='.', color='g')
					
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
				#---methode de minimum de Y-------
				for i in range(len(X_Fit3)):
					Y.append(min(Y_Fit3))
					X.append(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
					X_Fit3.remove(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
					Y_Fit3.remove(min(Y_Fit3))
					
				#---methode de minimum de X-------
				"""
				for i in range(len(X_Fit3)):
				    X.append(min(X_Fit3))
				    Y.append(Y_Fit3[X_Fit3.index(min(X_Fit3))])
				    Y_Fit3.remove(Y_Fit3[X_Fit3.index(min(X_Fit3))])
				    X_Fit3.remove(min(X_Fit3))
				"""
				#---------------------------------- 
				    
				print "len(X) = ", len(X)
				print "len(Y) = ", len(Y)
					
				print "X = ", X
				print "Y = ", Y
				if omega0_Model < min(X):
					omegalist2 = linspace(omega0_Model-0.02, max(X), 200)
					
				elif omega0_Model > max(X):
					omegalist2 = linspace(min(X), omega0_Model+0.02, 200)
						
				else:
					omegalist2 = linspace(min(X), max(X), 200)
						
						
				plt.plot(X, Y, color=couleur)
				
				Lamdalimit=[]
				omegalist = linspace(min(omegalist), max(omegalist), 200)
				
				for i in range(len(omegalist)):
					Lamdalimit.append(Functions(0.7, 0.3).NoBigBang(omegalist[i], 'limit'))

				plt.fill_between(omegalist, Lamdalimit, max(Lamdalimit), color=(0.5,0.5,0.5))
				plt.text(0.001, 1.468, 'No Big Bang', color='w')
				plt.text(0.7, 1.35, '($\\Omega_{0}$, $\\lambda_{0}$) = ('+str(omega0_Model)+', '+str(lambda0_Model)+')', color='r')
				
				plt.plot([omega0_Model],[lambda0_Model], marker='o', color='r')
				plt.xlabel('$\\Omega_{0}$', fontsize=16)
				plt.ylabel('$\\lambda_{0}$', fontsize=16)
				plt.title('The null correlation curve in the ($\\Omega_{0}$, $\\lambda_{0}$) plan', color='b')
				grid(True)
						
				plt.xlim(min(omegalist), max(omegalist))                                                                                                                                                                                                                                                                                                                                                                                                               
				plt.ylim(min(lamdalist), max(lamdalist)) 
				#plt.show()
				
				
				if KS == 'K-S_OK':
				  
					DD = []
					
				elif KS == 'K-S_NO':
				    
					DD = []
				      
				else:
					print "You must to determine the condition of the Kolmogrov-Smirnov test as 'K-S_OK' or 'K-S_NO'. ", "\n"

				

				
				#=======================================================================================================
				
				#======================================### Secondly for stretch ###=======================================
				omegalist = linspace(np.float64(10**-6), np.float64(1.0), 50) 
				lamdalist = linspace(np.float64(10**-6), np.float64(2.0), 50) 
		
				All_modelss = []
				for ih in range(len(lamdalist)):
					for jh in range(len(omegalist)):
						All_modelss.append((lamdalist[ih], omegalist[jh]))
				
				print "len(omegalist) = ", len(omegalist), "\n"
				print "len(lamdalist) = ", len(lamdalist), "\n"
				print "len(All_modelss) = ", len(All_modelss), "\n"
				
				def func2(a, z, bet):     # a = All_models;  z = redhsift;  bet = Beta
					aa=[]
					#beta = bet
					for i in range(len(a)):
						
						Bounced_limit = Functions(a[i][0],a[i][1]).NoBigBang(a[i][1], 'permission')       # si n'est pas valable, alors sera remplacer par le Modele standard (à voir après)
						#if lamdalist[i] >= Bounced_limit:
						if Bounced_limit==0:
							aa.append(Cov_s)
						else:
							#Mlist_tild = m_app - ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,0,0,0)
							Mlist_tild = MM*1.0
							#ZETA_tild = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,0,0,0)
							
							#ss = (( M_0 + bbeta*clist + ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,0,0,0) - m_app )/alfa )
							ss = slist*1.0

							#ten_power_ZETA_virt = (10**(ZETA_virt/5.))**bet
							
							#ten_power_s_m = (10**((m_app + alfa*(slist-1))/5.))**bet
							#ten_power_s_m = (10**((m_app - slist)/5.))**bet
							#ten_power_s_m = (10**((bbeta*clist)/5.))**bet
							
							
							WeightingManifold = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, bet)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
							
							#WeightingManifold = WeightingManifold*ten_power_s_m
							
							WeightingSum = [b/sum(b) for b in WeightingManifold]
							DetermineBeta = [abs(max(b) - min(b)) for b in WeightingSum]
							beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
							#L_1_Manifold = [1.0 + (1./np.log(len(b))*sum(b*np.log(b))) for b in WeightingSum]
							#beta = Beta[L_1_Manifold.index(min(L_1_Manifold))][0]

							#ten_power_s_m = (10**((m_app - slist)/5.))**beta
							#ten_power_s_m = (10**((bbeta*clist)/5.))**beta
							w_k = (ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta))/sum((ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta)))
							
							COVA_riance = sum(w_k*(ss - sum(w_k*ss))*(m_app - sum(w_k*m_app)))
							#COR_elation = (1./((len(slist)/(len(slist) - 1.))*sqrt(sum(w_k*(slist - sum(w_k*slist))**2.0))*sqrt(sum(w_k*(m_app - sum(w_k*m_app))**2.0))))*COVA_riance
							
							COR_elation = COVA_riance

							aa.append(COR_elation)

							#aa.append((w_k, Mlist_tild))
						
							#~plt.figure(44)
							#z2 = z*1.0
							#z2.sort()
							#ww = ApparentMagnitude(z2, a[i][0],a[i][1], 0.0).Weightingfactor(z2, beta, 0.0)
							#wk = ww/sum(ww)
							#plt.plot(z2, len(z2)*wk)
							#~for hg in range(len(w_k)):
							#~	plt.plot(z[hg], w_k[hg], marker='+')
								
							"""
							plt.figure(5)
							L_1 = 1.0 + (1./np.log(len(wk))*sum(wk*np.log(wk)))
							plt.plot(a[i][0], L_1, marker='+', color='g')
							plt.plot(a[i][1], L_1, marker='+', color='r')
							
							plt.figure(6)
							V = Model(a[i][0], a[i][1]).curvature(z2, a[i][0],a[i][1], 0.0)[1]
							plt.plot(z2, V)
							"""
					return aa


				
				Cov_All_Modelss = func2(All_modelss,redshift, Beta)
				
				#WeightingNormalisation_All_Models = func(All_models,redshift, Beta)
				#Cov_All_Models = [(1./(sqrt(sum(b[0]*(b[1] - sum(b[0]*b[1]))**2.0))*sqrt(sum(b[0]*(m_app - sum(b[0]*m_app))**2.0))))*sum(b[0]*(b[1] - sum(b[0]*b[1]))*(m_app - sum(b[0]*m_app))) for b in WeightingNormalisation_All_Models]	# (b[0], b[1]) = (w_k, Mlist_tild)
			
				# compute of mean of all true models
				#MeanCov_All_Models = np.mean(Cov_All_Models)
				MeanCov_All_Models = Cov_s
				
				#print "Cov_All_Models = ", Cov_All_Models, "\n"
				print "len(Cov_All_Modelss) = ", len(Cov_All_Modelss), "\n"
				print "max(Cov_All_Modelss) = ", max(Cov_All_Modelss), "\n"
				print "min(Cov_All_Modelss) = ", min(Cov_All_Modelss), "\n"
				print "MeanCov_All_Models = ", MeanCov_All_Models, "\n"
				
				print "len(All_modelss) = ", len(All_modelss), "\n"
				
				All_modelss = np.asarray(All_modelss)
				Cov_All_Modelss = np.asarray(Cov_All_Modelss)
				
				Cov_All_Modelss = Cov_All_Modelss - MeanCov_All_Models
				
				#***********grid of covariance**********************
				# figure of 3D(surface)
				
				#All_models[:,0] is a liste of lamda
				#All_models[:,1] is a liste of omega
				
				fig = plt.figure(7)
				ax = fig.add_subplot(111)

				col = plt.scatter(All_modelss[:,1], All_modelss[:,0], marker='.', s=150, c=Cov_All_Modelss, linewidths=0.3, cmap=plt.cm.Spectral_r)
				#for vv, ww, dd in zip(All_modelss[:,1], All_modelss[:,0], Cov_All_Modelss):
				#	pl.text(vv, ww, '%.3f' % dd, ha='center', va='bottom')
				

				# Add a colourbar.
				#cax = fig.colorbar(col, orientation='vertical',format='%.30f')
				cax = fig.colorbar(col, orientation='vertical',format='%.4f')
				cax.set_label('Covariance')
				
				plt.plot([omega0_Model],[lambda0_Model], marker='o', color='m')
				plt.xlabel('$\\Omega_{0}$', fontsize=14)
				plt.ylabel('$\\lambda_{0}$', fontsize=14)
				
				#===========For the save of the data====================
				OmegaOfGrid = All_modelss[:,1]*1.0
				LambdaOfGrid = All_modelss[:,0]*1.0
				CorrCoefficient = Cov_All_Modelss*1.0
				#=======================================================

				
				#plt.figure(8)
				#for i in range(len(Cov_All_Models)):
				#	plt.plot(Cov_All_Models[i], All_models[i][0], marker='o', color='r')
				
				#n, bins, patches = plt.hist(Cov_All_Models,40,normed='True')
				#(mu_hist, sigma_hist) = norm.fit(array(Cov_All_Models))
				#y = mlab.normpdf( bins, mu_hist, sigma_hist)
				#plt.plot(bins, y, 'r--', linewidth=2)
				#plt.title(r'$\rho_{0}=%.7f,\ \sigma_{0}=%.3f$' %(mu_hist, sigma_hist), color='r')
				#plt.xlabel('Correlation coefficient', fontsize=14)
				
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
				Cov_All_Modelss = Cov_All_Modelss.reshape(sqrt(len(All_modelss[:,1])), sqrt(len(All_modelss[:,0])))
				
				All_modelslambda = All_modelss[:,0].reshape(sqrt(len(All_modelss[:,0])), sqrt(len(All_modelss[:,0]))) # For matrix of lambda
				All_modelsomega = All_modelss[:,1].reshape(sqrt(len(All_modelss[:,1])), sqrt(len(All_modelss[:,0]))) # For matrix of omega

				#============**** To make the contours of confidence levels ****=============  
				#CS = pl.contour(All_modelsomega, All_modelslambda, abs(Cov_All_Models), 3, linewidths=np.arange(.5, 4, .5), colors=('r', 'green', 'blue', (1,1,0), '#afeeee', '0.5'))
				#plt.clabel(CS, inline=1, fontsize=10)
				#============***************************************************=============  
				#plt.show()
				#search of covariance zero by HORIZONTAL interpolation
				for i in range(len(omegalist)):

					#----------------------------------------------------
					
					Condition_1 = any(Cov_All_Modelss[i] == 0.0)
					if Condition_1 == True:
						rr = Cov_All_Modelss[i]
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
						rr = Cov_All_Modelss[i]
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
					
					Condition_1 = any(Cov_All_Modelss[:,i] == 0.0)
					if Condition_1 == True:
						rr = Cov_All_Modelss[:,i]
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
						rr = Cov_All_Modelss[:,i]
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
					
				
				#plt.figure(9)
				
				Y_Fit, X_Fit = [], []
				lambda0Covariance = lambda0Covariance_Hor + lambda0Covariance_Vert
				omega0Covariance = omega0Covariance_Hor + omega0Covariance_Vert
				#"""
				lambda0Covariance1 = list(array(lambda0Covariance)*1.0)
				omega0Covariance1 = list(array(omega0Covariance)*1.0)
				
				print "lambda0Covariance1 = ", lambda0Covariance1
				print "omega0Covariance1 = ", omega0Covariance1
				
				for i in range(len(lambda0Covariance_Vert)):
					plt.plot(omega0Covariance_Vert[i], lambda0Covariance_Vert[i], marker='.', color='r')
					
				for i in range(len(lambda0Covariance_Hor)):
					plt.plot(omega0Covariance_Hor[i], lambda0Covariance_Hor[i], marker='.', color='g')
					
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
				if omega0_Model < min(X):
					omegalist2 = linspace(omega0_Model-0.02, max(X), 200)
					
				elif omega0_Model > max(X):
					omegalist2 = linspace(min(X), omega0_Model+0.02, 200)
						
				else:
					omegalist2 = linspace(min(X), max(X), 200)
						
						
				plt.plot(X, Y, color=couleur)
				
				Lamdalimit=[]
				omegalist = linspace(min(omegalist), max(omegalist), 200)
				
				for i in range(len(omegalist)):
					Lamdalimit.append(Functions(0.7, 0.3).NoBigBang(omegalist[i], 'limit'))

				plt.fill_between(omegalist, Lamdalimit, max(Lamdalimit), color=(0.5,0.5,0.5))
				plt.text(0.001, 1.468, 'No Big Bang', color='w')
				plt.text(0.7, 1.32, '($\\Omega_{0}$, $\\lambda_{0}$) = ('+str(omega0_Model)+', '+str(lambda0_Model)+')', color='r')
				
				plt.plot([omega0_Model],[lambda0_Model], marker='o', color='r')
				plt.xlabel('$\\Omega_{0}$', fontsize=16)
				plt.ylabel('$\\lambda_{0}$', fontsize=16)
				plt.title('The null correlation curve in the ($\\Omega_{0}$, $\\lambda_{0}$) plan', color='b')
				grid(True)
						
				plt.xlim(min(omegalist), max(omegalist))                                                                                                                                                                                                                                                                                                                                                                                                               
				plt.ylim(min(lamdalist), max(lamdalist))      

				
				plt.figure(9)
				
				Y_Fit, X_Fit = [], []
				lambda0Covariance = lambda0Covariance_Hor + lambda0Covariance_Vert
				omega0Covariance = omega0Covariance_Hor + omega0Covariance_Vert
				lambda0Covariance.append(lambda0_Model)
				omega0Covariance.append(omega0_Model)
				#"""
				lambda0Covariance1 = list(array(lambda0Covariance)*1.0)
				omega0Covariance1 = list(array(omega0Covariance)*1.0)
				
				print "lambda0Covariance1 = ", lambda0Covariance1
				print "omega0Covariance1 = ", omega0Covariance1
				
				for i in range(len(lambda0Covariance_Vert)):
					plt.plot(omega0Covariance_Vert[i], lambda0Covariance_Vert[i], marker='.', color='r')
					
				for i in range(len(lambda0Covariance_Hor)):
					plt.plot(omega0Covariance_Hor[i], lambda0Covariance_Hor[i], marker='.', color='g')
					
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
				#---methode de minimum de Y-------
				for i in range(len(X_Fit3)):
					Y.append(min(Y_Fit3))
					X.append(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
					X_Fit3.remove(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
					Y_Fit3.remove(min(Y_Fit3))
					
				#---methode de minimum de X-------
				"""
				for i in range(len(X_Fit3)):
				    X.append(min(X_Fit3))
				    Y.append(Y_Fit3[X_Fit3.index(min(X_Fit3))])
				    Y_Fit3.remove(Y_Fit3[X_Fit3.index(min(X_Fit3))])
				    X_Fit3.remove(min(X_Fit3))
				"""
				#---------------------------------- 
				    
				print "len(X) = ", len(X)
				print "len(Y) = ", len(Y)
					
				print "X = ", X
				print "Y = ", Y
				if omega0_Model < min(X):
					omegalist2 = linspace(omega0_Model-0.02, max(X), 200)
					
				elif omega0_Model > max(X):
					omegalist2 = linspace(min(X), omega0_Model+0.02, 200)
						
				else:
					omegalist2 = linspace(min(X), max(X), 200)
						
						
				plt.plot(X, Y, color=couleur)
				
				Lamdalimit=[]
				omegalist = linspace(min(omegalist), max(omegalist), 200)
				
				for i in range(len(omegalist)):
					Lamdalimit.append(Functions(0.7, 0.3).NoBigBang(omegalist[i], 'limit'))

				plt.fill_between(omegalist, Lamdalimit, max(Lamdalimit), color=(0.5,0.5,0.5))
				plt.text(0.001, 1.468, 'No Big Bang', color='w')
				plt.text(0.7, 1.35, '($\\Omega_{0}$, $\\lambda_{0}$) = ('+str(omega0_Model)+', '+str(lambda0_Model)+')', color='r')
				
				plt.plot([omega0_Model],[lambda0_Model], marker='o', color='r')
				plt.xlabel('$\\Omega_{0}$', fontsize=16)
				plt.ylabel('$\\lambda_{0}$', fontsize=16)
				plt.title('The null correlation curve in the ($\\Omega_{0}$, $\\lambda_{0}$) plan', color='b')
				grid(True)
						
				plt.xlim(min(omegalist), max(omegalist))                                                                                                                                                                                                                                                                                                                                                                                                               
				plt.ylim(min(lamdalist), max(lamdalist)) 
				#plt.show()
				
				
				if KS == 'K-S_OK':
				  
					DD = []
					
				elif KS == 'K-S_NO':
				    
					DD = []
				      
				else:
					print "You must to determine the condition of the Kolmogrov-Smirnov test as 'K-S_OK' or 'K-S_NO'. ", "\n"

				return OmegaOfGrid, LambdaOfGrid, CorrCoefficient, X, Y, DD

				
				#=======================================================================================================
				
			elif Variables == 'stretch_only':
				All_models = []
				# regroupement de deux valeurs (lamda, omega) pour représenter un model cosmologique
				for i in range(len(lamdalist)):
				    for j in range(len(omegalist)):
					    All_models.append((lamdalist[i], omegalist[j]))

				def func(a, z, bet):     # a = All_models;  z = redhsift;  bet = Beta
					aa=[]
					#beta = bet
					for i in range(len(a)):
						
						Bounced_limit = Functions(a[i][0],a[i][1]).NoBigBang(a[i][1], 'permission')       # si n'est pas valable, alors sera remplacer par le Modele standard (à voir après)
						#if lamdalist[i] >= Bounced_limit:
						if Bounced_limit==0:
							aa.append(Cov_s)
						else:
							Mlist_tild = m_app - ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,0,0,0)
							#Mlist_tild = MM*1.0
							#ZETA_tild = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,0,0,0)
							
							ss = (( M_0 + bbeta*clist + ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,0,0,0) - m_app )/alfa ) + 1.0
							#ss = slist*1.0

							#ten_power_ZETA_virt = (10**(ZETA_virt/5.))**bet
							
							#ten_power_s_m = (10**((m_app + alfa*(slist-1))/5.))**bet
							#ten_power_s_m = (10**((m_app - slist)/5.))**bet
							ten_power_s_m = (10**((bbeta*clist)/5.))**bet
							
							
							WeightingManifold = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, bet)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
							
							#WeightingManifold = WeightingManifold*ten_power_s_m
							
							WeightingSum = [b/sum(b) for b in WeightingManifold]
							DetermineBeta = [abs(max(b) - min(b)) for b in WeightingSum]
							beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
							#L_1_Manifold = [1.0 + (1./np.log(len(b))*sum(b*np.log(b))) for b in WeightingSum]
							#beta = Beta[L_1_Manifold.index(min(L_1_Manifold))][0]

							#ten_power_s_m = (10**((m_app - slist)/5.))**beta
							#ten_power_s_m = (10**((bbeta*clist)/5.))**beta
							w_k = (ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta))/sum((ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta)))
							
							#COVA_riance = sum(w_k*(ss - sum(w_k*ss))*(m_app - sum(w_k*m_app)))
							#COR_elation = (1./((len(slist)/(len(slist) - 1.))*sqrt(sum(w_k*(slist - sum(w_k*slist))**2.0))*sqrt(sum(w_k*(m_app - sum(w_k*m_app))**2.0))))*COVA_riance
							
							COVA_riance = sum(w_k*(Mlist_tild - sum(w_k*Mlist_tild ))*(m_app - sum(w_k*m_app)))
							COR_elation = COVA_riance

							aa.append(COR_elation)

							#aa.append((w_k, Mlist_tild))
						
							#~plt.figure(44)
							#z2 = z*1.0
							#z2.sort()
							#ww = ApparentMagnitude(z2, a[i][0],a[i][1], 0.0).Weightingfactor(z2, beta, 0.0)
							#wk = ww/sum(ww)
							#plt.plot(z2, len(z2)*wk)
							#~for hg in range(len(w_k)):
							#~	plt.plot(z[hg], w_k[hg], marker='+')
								
							"""
							plt.figure(5)
							L_1 = 1.0 + (1./np.log(len(wk))*sum(wk*np.log(wk)))
							plt.plot(a[i][0], L_1, marker='+', color='g')
							plt.plot(a[i][1], L_1, marker='+', color='r')
							
							plt.figure(6)
							V = Model(a[i][0], a[i][1]).curvature(z2, a[i][0],a[i][1], 0.0)[1]
							plt.plot(z2, V)
							"""
					return aa


				
				Cov_All_Models = func(All_models,redshift, Beta)
				
				#WeightingNormalisation_All_Models = func(All_models,redshift, Beta)
				#Cov_All_Models = [(1./(sqrt(sum(b[0]*(b[1] - sum(b[0]*b[1]))**2.0))*sqrt(sum(b[0]*(m_app - sum(b[0]*m_app))**2.0))))*sum(b[0]*(b[1] - sum(b[0]*b[1]))*(m_app - sum(b[0]*m_app))) for b in WeightingNormalisation_All_Models]	# (b[0], b[1]) = (w_k, Mlist_tild)
			
				# compute of mean of all true models
				#MeanCov_All_Models = np.mean(Cov_All_Models)
				MeanCov_All_Models = Cov_s
				
				#print "Cov_All_Models = ", Cov_All_Models, "\n"
				print "len(Cov_All_Models) = ", len(Cov_All_Models), "\n"
				print "max(Cov_All_Models) = ", max(Cov_All_Models), "\n"
				print "min(Cov_All_Models) = ", min(Cov_All_Models), "\n"
				print "MeanCov_All_Models = ", MeanCov_All_Models, "\n"
				
				print "len(All_models) = ", len(All_models), "\n"
				
				All_models = np.asarray(All_models)
				Cov_All_Models = np.asarray(Cov_All_Models)
				
				Cov_All_Models = Cov_All_Models - MeanCov_All_Models
				
				#***********grid of covariance**********************
				# figure of 3D(surface)
				
				#All_models[:,0] is a liste of lamda
				#All_models[:,1] is a liste of omega
				
				fig = plt.figure(7)
				ax = fig.add_subplot(111)

				col = plt.scatter(All_models[:,1], All_models[:,0], marker='.', s=150, c=Cov_All_Models, linewidths=0.3, cmap=plt.cm.Spectral_r)
				#for vv, ww, dd in zip(All_models[:,1], All_models[:,0], Cov_All_Models):
				#	pl.text(vv, ww, '%.3f' % dd, ha='center', va='bottom')
				

				# Add a colourbar.
				#cax = fig.colorbar(col, orientation='vertical',format='%.30f')
				cax = fig.colorbar(col, orientation='vertical',format='%.4f')
				cax.set_label('Covariance')
				
				plt.plot([omega0_Model],[lambda0_Model], marker='o', color='m')
				plt.xlabel('$\\Omega_{0}$', fontsize=14)
				plt.ylabel('$\\lambda_{0}$', fontsize=14)
				
				#===========For the save of the data====================
				OmegaOfGrid = All_models[:,1]*1.0
				LambdaOfGrid = All_models[:,0]*1.0
				CorrCoefficient = Cov_All_Models*1.0
				#=======================================================

				
				#plt.figure(8)
				#for i in range(len(Cov_All_Models)):
				#	plt.plot(Cov_All_Models[i], All_models[i][0], marker='o', color='r')
				
				#n, bins, patches = plt.hist(Cov_All_Models,40,normed='True')
				#(mu_hist, sigma_hist) = norm.fit(array(Cov_All_Models))
				#y = mlab.normpdf( bins, mu_hist, sigma_hist)
				#plt.plot(bins, y, 'r--', linewidth=2)
				#plt.title(r'$\rho_{0}=%.7f,\ \sigma_{0}=%.3f$' %(mu_hist, sigma_hist), color='r')
				#plt.xlabel('Correlation coefficient', fontsize=14)
				
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
					
				
				#plt.figure(9)
				
				Y_Fit, X_Fit = [], []
				lambda0Covariance = lambda0Covariance_Hor + lambda0Covariance_Vert
				omega0Covariance = omega0Covariance_Hor + omega0Covariance_Vert
				#"""
				lambda0Covariance1 = list(array(lambda0Covariance)*1.0)
				omega0Covariance1 = list(array(omega0Covariance)*1.0)
				
				print "lambda0Covariance1 = ", lambda0Covariance1
				print "omega0Covariance1 = ", omega0Covariance1
				
				for i in range(len(lambda0Covariance_Vert)):
					plt.plot(omega0Covariance_Vert[i], lambda0Covariance_Vert[i], marker='.', color='r')
					
				for i in range(len(lambda0Covariance_Hor)):
					plt.plot(omega0Covariance_Hor[i], lambda0Covariance_Hor[i], marker='.', color='g')
					
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
				if omega0_Model < min(X):
					omegalist2 = linspace(omega0_Model-0.02, max(X), 200)
					
				elif omega0_Model > max(X):
					omegalist2 = linspace(min(X), omega0_Model+0.02, 200)
						
				else:
					omegalist2 = linspace(min(X), max(X), 200)
						
						
				plt.plot(X, Y, color=couleur)
				
				Lamdalimit=[]
				omegalist = linspace(min(omegalist), max(omegalist), 200)
				
				for i in range(len(omegalist)):
					Lamdalimit.append(Functions(0.7, 0.3).NoBigBang(omegalist[i], 'limit'))

				plt.fill_between(omegalist, Lamdalimit, max(Lamdalimit), color=(0.5,0.5,0.5))
				plt.text(0.001, 1.468, 'No Big Bang', color='w')
				plt.text(0.7, 1.32, '($\\Omega_{0}$, $\\lambda_{0}$) = ('+str(omega0_Model)+', '+str(lambda0_Model)+')', color='r')
				
				plt.plot([omega0_Model],[lambda0_Model], marker='o', color='r')
				plt.xlabel('$\\Omega_{0}$', fontsize=16)
				plt.ylabel('$\\lambda_{0}$', fontsize=16)
				plt.title('The null correlation curve in the ($\\Omega_{0}$, $\\lambda_{0}$) plan', color='b')
				grid(True)
						
				plt.xlim(min(omegalist), max(omegalist))                                                                                                                                                                                                                                                                                                                                                                                                               
				plt.ylim(min(lamdalist), max(lamdalist))      

				
				plt.figure(9)
				
				Y_Fit, X_Fit = [], []
				lambda0Covariance = lambda0Covariance_Hor + lambda0Covariance_Vert
				omega0Covariance = omega0Covariance_Hor + omega0Covariance_Vert
				lambda0Covariance.append(lambda0_Model)
				omega0Covariance.append(omega0_Model)
				#"""
				lambda0Covariance1 = list(array(lambda0Covariance)*1.0)
				omega0Covariance1 = list(array(omega0Covariance)*1.0)
				
				print "lambda0Covariance1 = ", lambda0Covariance1
				print "omega0Covariance1 = ", omega0Covariance1
				
				for i in range(len(lambda0Covariance_Vert)):
					plt.plot(omega0Covariance_Vert[i], lambda0Covariance_Vert[i], marker='.', color='r')
					
				for i in range(len(lambda0Covariance_Hor)):
					plt.plot(omega0Covariance_Hor[i], lambda0Covariance_Hor[i], marker='.', color='g')
					
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
				#---methode de minimum de Y-------
				for i in range(len(X_Fit3)):
					Y.append(min(Y_Fit3))
					X.append(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
					X_Fit3.remove(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
					Y_Fit3.remove(min(Y_Fit3))
					
				#---methode de minimum de X-------
				"""
				for i in range(len(X_Fit3)):
				    X.append(min(X_Fit3))
				    Y.append(Y_Fit3[X_Fit3.index(min(X_Fit3))])
				    Y_Fit3.remove(Y_Fit3[X_Fit3.index(min(X_Fit3))])
				    X_Fit3.remove(min(X_Fit3))
				"""
				#---------------------------------- 
				    
				print "len(X) = ", len(X)
				print "len(Y) = ", len(Y)
					
				print "X = ", X
				print "Y = ", Y
				if omega0_Model < min(X):
					omegalist2 = linspace(omega0_Model-0.02, max(X), 200)
					
				elif omega0_Model > max(X):
					omegalist2 = linspace(min(X), omega0_Model+0.02, 200)
						
				else:
					omegalist2 = linspace(min(X), max(X), 200)
						
						
				plt.plot(X, Y, color=couleur)
				
				Lamdalimit=[]
				omegalist = linspace(min(omegalist), max(omegalist), 200)
				
				for i in range(len(omegalist)):
					Lamdalimit.append(Functions(0.7, 0.3).NoBigBang(omegalist[i], 'limit'))

				plt.fill_between(omegalist, Lamdalimit, max(Lamdalimit), color=(0.5,0.5,0.5))
				plt.text(0.001, 1.468, 'No Big Bang', color='w')
				plt.text(0.7, 1.35, '($\\Omega_{0}$, $\\lambda_{0}$) = ('+str(omega0_Model)+', '+str(lambda0_Model)+')', color='r')
				
				plt.plot([omega0_Model],[lambda0_Model], marker='o', color='r')
				plt.xlabel('$\\Omega_{0}$', fontsize=16)
				plt.ylabel('$\\lambda_{0}$', fontsize=16)
				plt.title('The null correlation curve in the ($\\Omega_{0}$, $\\lambda_{0}$) plan', color='b')
				grid(True)
						
				plt.xlim(min(omegalist), max(omegalist))                                                                                                                                                                                                                                                                                                                                                                                                               
				plt.ylim(min(lamdalist), max(lamdalist)) 
				#plt.show()
				
				
				if KS == 'K-S_OK':
				  
					DD = []
					
				elif KS == 'K-S_NO':
				    
					DD = []
				      
				else:
					print "You must to determine the condition of the Kolmogrov-Smirnov test as 'K-S_OK' or 'K-S_NO'. ", "\n"

				return OmegaOfGrid, LambdaOfGrid, CorrCoefficient, X, Y, DD

			elif Variables == 'color_only':
				All_models = []
				for i in range(len(lamdalist)):
					for j in range(len(omegalist)):
						All_models.append((lamdalist[i], omegalist[j]))

				def func(a, z, bet):     # a = All_models;  z = redhsift;  bet = Beta
					aa=[]
					#beta = bet
					for i in range(len(a)):
						
						Bounced_limit = Functions(a[i][0],a[i][1]).NoBigBang(a[i][1], 'permission')       # si n'est pas valable, alors sera remplacer par le Modele standard (à voir après)
						#if lamdalist[i] >= Bounced_limit:
						if Bounced_limit==0:
							aa.append(Cov_c)
						else:
							Mlist_tild = m_app - ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,bbeta,0,0)
							#@Mlist_tild = MM*1.0
							#ZETA_tild = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,0,0,0)
							
							#@ten_power_m_M = (10**((m_app - Mlist_tild)/5.))**bet
							#@Jacob = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Jacobian(z)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
							#@WeightingManifold = Jacob * ten_power_m_M
							WeightingManifold = ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, bet)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
							
							WeightingSum = [b/sum(b) for b in WeightingManifold]
							DetermineBeta = [abs((max(b) - min(b))/(max(b) + min(b))) for b in WeightingSum]
							beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
							#L_1_Manifold = [1.0 + (1./np.log(len(b))*sum(b*np.log(b))) for b in WeightingSum]
							#beta = Beta[L_1_Manifold.index(min(L_1_Manifold))][0]
							beta_pliot = beta*1.0

							w_k = (ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta))/sum((ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta)))
							#@ten_power_m_M_beta = (10**((m_app - Mlist_tild)/5.))**beta
							#@w_k = (Jacob*ten_power_m_M_beta)/sum(Jacob*ten_power_m_M_beta)

							COVA_riance = sum(w_k*(Mlist_tild - sum(w_k*Mlist_tild ))*(m_app - sum(w_k*m_app)))

							COR_elation = COVA_riance
							#COR_elation = (1./((len(Mlist_tild)/(len(Mlist_tild) - 1.))*sqrt(sum(w_k*(Mlist_tild - sum(w_k*Mlist_tild ))**2.0))*sqrt(sum(w_k*(m_app - sum(w_k*m_app))**2.0))))*COVA_riance

							aa.append(COR_elation)

							#aa.append((w_k, Mlist_tild))
						
							#~plt.figure(44)
							#z2 = z*1.0
							#z2.sort()
							#ww = ApparentMagnitude(z2, a[i][0],a[i][1], 0.0).Weightingfactor(z2, beta, 0.0)
							#wk = ww/sum(ww)
							#plt.plot(z2, len(z2)*wk)
							#~for hg in range(len(w_k)):
							#~	plt.plot(z[hg], w_k[hg], marker='+')
								
							"""
							plt.figure(5)
							L_1 = 1.0 + (1./np.log(len(wk))*sum(wk*np.log(wk)))
							plt.plot(a[i][0], L_1, marker='+', color='g')
							plt.plot(a[i][1], L_1, marker='+', color='r')
							
							plt.figure(6)
							V = Model(a[i][0], a[i][1]).curvature(z2, a[i][0],a[i][1], 0.0)[1]
							plt.plot(z2, V)
							"""
					return aa

				m_app_no_correct = m_app*1.0
				
				#Heaviside = Functions(0.7, 0.3).Heavisidefuction(mass_hst-10)
				#mass_hst = t['3rdvar']	#log10(Mass_of_host_galaxy) it must be > 10
				Cov_bbeta = []
				#for ij in range(len(bbetalist)):
				for ij in range(1):
					bbetalist[ij] = bbeta*1.0
				
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
					# figure of 3D(surface)
					
					#All_models[:,0] is a liste of lamda
					#All_models[:,1] is a liste of omega
					
					fig = plt.figure(44)
					ax = fig.add_subplot(111)

					col = plt.scatter(All_models[:,1], All_models[:,0], marker='.', s=150, c=Cov_All_Models, linewidths=0.3, cmap=plt.cm.Spectral_r)
					#for vv, ww, dd in zip(All_models[:,1], All_models[:,0], Cov_All_Models):
					#	pl.text(vv, ww, '%.3f' % dd, ha='center', va='bottom')
					

					# Add a colourbar.
					#cax = fig.colorbar(col, orientation='vertical',format='%.30f')
					cax = fig.colorbar(col, orientation='vertical',format='%.4f')
					cax.set_label('Covariance')
					
					plt.plot([omega0_Model],[lambda0_Model], marker='o', color='m')
					plt.xlabel('$\\Omega_{0}$', fontsize=14)
					plt.ylabel('$\\lambda_{0}$', fontsize=14)
					
					#===========For the save of the data====================
					OmegaOfGrid = All_models[:,1]*1.0
					LambdaOfGrid = All_models[:,0]*1.0
					CorrCoefficient = Cov_All_Models*1.0
					#=======================================================

					
					#plt.figure(8)
					#for i in range(len(Cov_All_Models)):
					#	plt.plot(Cov_All_Models[i], All_models[i][0], marker='o', color='r')
					
					#n, bins, patches = plt.hist(Cov_All_Models,40,normed='True')
					#(mu_hist, sigma_hist) = norm.fit(array(Cov_All_Models))
					#y = mlab.normpdf( bins, mu_hist, sigma_hist)
					#plt.plot(bins, y, 'r--', linewidth=2)
					#plt.title(r'$\rho_{0}=%.7f,\ \sigma_{0}=%.3f$' %(mu_hist, sigma_hist), color='r')
					#plt.xlabel('Correlation coefficient', fontsize=14)
					
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
					Cov_All_Models = Cov_All_Models.reshape(int(sqrt(len(All_models[:,1]))), int(sqrt(len(All_models[:,0]))))
					
					All_modelslambda = All_models[:,0].reshape(int(sqrt(len(All_models[:,0]))), int(sqrt(len(All_models[:,0])))) # For matrix of lambda
					All_modelsomega = All_models[:,1].reshape(int(sqrt(len(All_models[:,1]))), int(sqrt(len(All_models[:,0])))) # For matrix of omega

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
						
					
					#plt.figure(9)
					
					Y_Fit, X_Fit = [], []
					lambda0Covariance = lambda0Covariance_Hor + lambda0Covariance_Vert
					omega0Covariance = omega0Covariance_Hor + omega0Covariance_Vert
					#"""
					lambda0Covariance1 = list(array(lambda0Covariance)*1.0)
					omega0Covariance1 = list(array(omega0Covariance)*1.0)
					
					if len(lambda0Covariance1) > 3 and len(omega0Covariance1) > 3:
						print "lambda0Covariance1 = ", lambda0Covariance1
						print "omega0Covariance1 = ", omega0Covariance1
						
						for i in range(len(lambda0Covariance_Vert)):
							plt.plot(omega0Covariance_Vert[i], lambda0Covariance_Vert[i], marker='.', color='r')
							
						for i in range(len(lambda0Covariance_Hor)):
							plt.plot(omega0Covariance_Hor[i], lambda0Covariance_Hor[i], marker='.', color='g')
							
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
						if omega0_Model < min(X):
							omegalist2 = linspace(omega0_Model-0.02, max(X), 200)
							
						elif omega0_Model > max(X):
							omegalist2 = linspace(min(X), omega0_Model+0.02, 200)
								
						else:
							omegalist2 = linspace(min(X), max(X), 200)
								
								
						plt.plot(X, Y, color=couleur)
						
						Lamdalimit=[]
						omegalist = linspace(min(omegalist), max(omegalist), 200)
						
						for i in range(len(omegalist)):
							Lamdalimit.append(Functions(0.7, 0.3).NoBigBang(omegalist[i], 'limit'))

						plt.fill_between(omegalist, Lamdalimit, max(Lamdalimit), color=(0.5,0.5,0.5))
						plt.text(0.001, 1.468, 'No Big Bang', color='w')
						plt.text(0.7, 1.32, '($\\Omega_{0}$, $\\lambda_{0}$) = ('+str(omega0_Model)+', '+str(lambda0_Model)+')', color='r')
						
						plt.plot([omega0_Model],[lambda0_Model], marker='o', color='r')
						plt.xlabel('$\\Omega_{0}$', fontsize=16)
						plt.ylabel('$\\lambda_{0}$', fontsize=16)
						plt.title('The null correlation curve in the ($\\Omega_{0}$, $\\lambda_{0}$) plan', color='b')
						grid(True)
								
						plt.xlim(min(omegalist), max(omegalist))                                                                                                                                                                                                                                                                                                                                                                                                               
						plt.ylim(min(lamdalist), max(lamdalist))      

						
						plt.figure(77)
						
						Y_Fit, X_Fit = [], []
						lambda0Covariance = lambda0Covariance_Hor + lambda0Covariance_Vert
						omega0Covariance = omega0Covariance_Hor + omega0Covariance_Vert
						lambda0Covariance.append(lambda0_Model)
						omega0Covariance.append(omega0_Model)
						#"""
						lambda0Covariance1 = list(array(lambda0Covariance)*1.0)
						omega0Covariance1 = list(array(omega0Covariance)*1.0)
						
						print "lambda0Covariance1 = ", lambda0Covariance1
						print "omega0Covariance1 = ", omega0Covariance1
						
						for i in range(len(lambda0Covariance_Vert)):
							plt.plot(omega0Covariance_Vert[i], lambda0Covariance_Vert[i], marker='.', color='r')
							
						for i in range(len(lambda0Covariance_Hor)):
							plt.plot(omega0Covariance_Hor[i], lambda0Covariance_Hor[i], marker='.', color='g')
							
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
						#---methode de minimum de Y-------
						for i in range(len(X_Fit3)):
							Y.append(min(Y_Fit3))
							X.append(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
							X_Fit3.remove(X_Fit3[Y_Fit3.index(min(Y_Fit3))])
							Y_Fit3.remove(min(Y_Fit3))
							
						#---methode de minimum de X-------
						"""
						for i in range(len(X_Fit3)):
						    X.append(min(X_Fit3))
						    Y.append(Y_Fit3[X_Fit3.index(min(X_Fit3))])
						    Y_Fit3.remove(Y_Fit3[X_Fit3.index(min(X_Fit3))])
						    X_Fit3.remove(min(X_Fit3))
						"""
						#---------------------------------- 
						    
						print "len(X) = ", len(X)
						print "len(Y) = ", len(Y)
							
						print "X = ", X
						print "Y = ", Y
						if omega0_Model < min(X):
							omegalist2 = linspace(omega0_Model-0.02, max(X), 200)
							
						elif omega0_Model > max(X):
							omegalist2 = linspace(min(X), omega0_Model+0.02, 200)
								
						else:
							omegalist2 = linspace(min(X), max(X), 200)
								
								
						plt.plot(X, Y, color=couleur)
						
						Lamdalimit=[]
						omegalist = linspace(min(omegalist), max(omegalist), 200)
						
						for i in range(len(omegalist)):
							Lamdalimit.append(Functions(0.7, 0.3).NoBigBang(omegalist[i], 'limit'))

						plt.fill_between(omegalist, Lamdalimit, max(Lamdalimit), color=(0.5,0.5,0.5))
						plt.text(0.001, 1.468, 'No Big Bang', color='w')
						plt.text(0.7, 1.35, '($\\Omega_{0}$, $\\lambda_{0}$) = ('+str(omega0_Model)+', '+str(lambda0_Model)+')', color='r')
						plt.text(0.7, 0.9, 'bbetalist[ij] = '+str(bbetalist[ij]), color='g')
						
						plt.plot([omega0_Model],[lambda0_Model], marker='o', color='r')
						plt.xlabel('$\\Omega_{0}$', fontsize=16)
						plt.ylabel('$\\lambda_{0}$', fontsize=16)
						plt.title('The null correlation curve in the ($\\Omega_{0}$, $\\lambda_{0}$) plan', color='b')
						grid(True)
								
						plt.xlim(min(omegalist), max(omegalist))                                                                                                                                                                                                                                                                                                                                                                                                               
						plt.ylim(min(lamdalist), max(lamdalist)) 
						#plt.show()
					else:
						pass

				
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
				#======================================### Firstly for color ###======================================
				All_models = []
				# regroupement de deux valeurs (alfa, bbeta) pour représenter un magnitude absolue M
				for i in range(len(alfalist)):
				    for j in range(len(bbetalist)):
					    All_models.append((alfalist[i], bbetalist[j]))
				   
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

							#aa.append((w_k, Mlist_tild))
						
							plt.figure(44)
							z2 = z*1.0
							z2.sort()
							ww = ApparentMagnitude(z2, a[i][0],a[i][1], 0.0, Realdata, 0, 0).Weightingfactor(z2, beta)
							wk = ww/sum(ww)
							plt.plot(z2, len(z2)*wk)
							plt.text(min(z2), max(len(z2)*wk), '$\\gamma$ = ('+str(beta_pliot)+')')
							#for hg in range(len(w_k)):
							#	plt.plot(z[hg], w_k[hg], marker='+')
								
							"""
							plt.figure(5)
							L_1 = 1.0 + (1./np.log(len(wk))*sum(wk*np.log(wk)))
							plt.plot(a[i][0], L_1, marker='+', color='g')
							plt.plot(a[i][1], L_1, marker='+', color='r')
							
							plt.figure(6)
							V = Model(a[i][0], a[i][1]).curvature(z2, a[i][0],a[i][1], 0.0)[1]
							plt.plot(z2, V)
							"""
					return aa, beta_pliot

				m_app_no_correct = m_app*1.0
				
				#Heaviside = Functions(0.7, 0.3).Heavisidefuction(mass_hst-10)
				#mass_hst = t['3rdvar']	#log10(Mass_of_host_galaxy) it must be > 10
				Cov_bbeta = []
				for ij in range(len(bbetalist)):
				#for ij in range(1):
					#bbetalist[ij] = bbeta*1.0
				
					m_app = m_app_no_correct - bbetalist[ij]*clist
					All_models = [(lambda0_Model, omega0_Model)]	#si je veux faire un test pour determiner bbeta. Si je cherche les parametres cosmologique, il faut reprendre All_models initial (la grille)
					
					Cov_All_Models_operation = func(All_models,redshift, Beta)
					Cov_All_Models = Cov_All_Models_operation[0]
					beta_pliot = Cov_All_Models_operation[1]
					
					#WeightingNormalisation_All_Models = func(All_models,redshift, Beta)
					#Cov_All_Models = [(1./(sqrt(sum(b[0]*(b[1] - sum(b[0]*b[1]))**2.0))*sqrt(sum(b[0]*(m_app - sum(b[0]*m_app))**2.0))))*sum(b[0]*(b[1] - sum(b[0]*b[1]))*(m_app - sum(b[0]*m_app))) for b in WeightingNormalisation_All_Models]	# (b[0], b[1]) = (w_k, Mlist_tild)
				
					# compute of mean of all true models
					#MeanCov_All_Models = np.mean(Cov_All_Models)
					MeanCov_All_Models = Cov_c
					#MeanCov_All_Models = 0
					
					#print "Cov_All_Models = ", Cov_All_Models, "\n"
					print "len(Cov_All_Models) = ", len(Cov_All_Models), "\n"
					print "max(Cov_All_Models) = ", max(Cov_All_Models), "\n"
					print "min(Cov_All_Models) = ", min(Cov_All_Models), "\n"
					print "MeanCov_All_Models = ", MeanCov_All_Models, "\n"
					
					print "len(All_models) = ", len(All_models), "\n"
					print "bbetalist[ij] = ", bbetalist[ij], "\n"
					print "gamma = ", beta_pliot, "\n"
					print "bbeta = ", bbeta, "\n"
					
					All_models = np.asarray(All_models)
					Cov_All_Models = np.asarray(Cov_All_Models)
					
					Cov_All_Models = Cov_All_Models - MeanCov_All_Models
					Cov_bbeta.append(Cov_All_Models[0])
					
					plt.figure(152)
					plt.plot(bbetalist[ij], Cov_All_Models[0], marker='.', color='r')

					
				plt.plot(bbetalist, Cov_bbeta, color='k')
				plt.plot(bbeta, 0, marker='o', color='g')
				#plt.text(4, mean(Cov_bbeta), '$cov_initial$ = ('+str(Cov_c)+')', color='m')
				#plt.text(2, mean(Cov_bbeta), '$\\gamma$ = ('+str(beta_pliot)+')', color='b')
				
				plt.xlabel('$\\beta$', fontsize=14)
				plt.ylabel('$\\rho$', fontsize=14)
				
				CorrCoefficient = array(Cov_bbeta)*1.0
				X, Y = [], []
				
				if KS == 'K-S_OK':
				  
					DD = []
					
				elif KS == 'K-S_NO':
				    
					DD = []
				      
				else:
					print "You must to determine the condition of the Kolmogrov-Smirnov test as 'K-S_OK' or 'K-S_NO'. ", "\n"

				
				#=======================================================================================================
				
				#======================================### Secondly for stretch ###=======================================
				X, Y = [], []
				
				if KS == 'K-S_OK':
				  
					DD = []
					
				elif KS == 'K-S_NO':
				    
					DD = []
				      
				else:
					print "You must to determine the condition of the Kolmogrov-Smirnov test as 'K-S_OK' or 'K-S_NO'. ", "\n"

				return bbetalist, CorrCoefficient, X, Y, DD


				
				#=======================================================================================================
				
			elif Variables == 'stretch_only':
				All_models = []
				# regroupement de deux valeurs (alfa, bbeta) pour représenter un magnitude absolue M
				for i in range(len(alfalist)):
				    for j in range(len(bbetalist)):
					    All_models.append((alfalist[i], bbetalist[j]))

				def func(a, z, bet):     # a = All_models;  z = redhsift;  bet = Beta
					aa=[]
					#beta = bet
					for i in range(len(a)):

						#Mlist_tild = m_app - ApparentMagnitude(z, lambda0_Model, omega0_Model, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, lambda0_Model, omega0_Model, 0,0,0,0,0)
						#Mlist_tild = M_0 - a[i][0]*(slist-1) + a[i][1]*clist
						
						ss = (( M_0 + bbeta*clist + ApparentMagnitude(z, a[i][0],a[i][1], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, a[i][0],a[i][1], 0,0,0,0,0) - m_app )/alfa ) + 1.0
						#ss = slist*1.0

						#ten_power_s_m = (10**((m_app + a[i][0]*(slist-1))/5.))**bet
						#ten_power_s_m = (10**((m_app - slist)/5.))**bet
						ten_power_s_m = (10**((a[i][1]*clist)/5.))**bet

						WeightingManifold = ApparentMagnitude(z, lambda0_Model, omega0_Model, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, bet)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
						
						WeightingManifold = WeightingManifold*ten_power_s_m
						
						WeightingSum = [b/sum(b) for b in WeightingManifold]
						DetermineBeta = [abs(max(b) - min(b)) for b in WeightingSum]
						beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
						#L_1_Manifold = [1.0 + (1./np.log(len(b))*sum(b*np.log(b))) for b in WeightingSum]
						#beta = Beta[L_1_Manifold.index(min(L_1_Manifold))][0]
			
						#ten_power_s_m = (10**((m_app - slist)/5.))**beta
						ten_power_s_m = (10**((a[i][1]*clist)/5.))**beta
						w_k = (ApparentMagnitude(z, lambda0_Model, omega0_Model, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta)*ten_power_s_m)/sum((ApparentMagnitude(z, lambda0_Model, omega0_Model, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta)*ten_power_s_m))
							
						COVA_riance = sum(w_k*(ss - sum(w_k*ss))*(m_app - sum(w_k*m_app)))
						#COR_elation = (1./((len(slist)/(len(slist) - 1.))*sqrt(sum(w_k*(slist - sum(w_k*slist))**2.0))*sqrt(sum(w_k*(m_app - sum(w_k*m_app))**2.0))))*COVA_riance

						COR_elation = COVA_riance

						aa.append(COR_elation)

						#aa.append((w_k, Mlist_tild))
						
						#~plt.figure(44)
						#z2 = z*1.0
						#z2.sort()
						#ww = ApparentMagnitude(z2, a[i][0],a[i][1], 0.0).Weightingfactor(z2, beta, 0.0)
						#wk = ww/sum(ww)
						#plt.plot(z2, len(z2)*wk)
						#~for hg in range(len(w_k)):
						#~	plt.plot(z[hg], w_k[hg], marker='+')
							
						"""
						plt.figure(5)
						L_1 = 1.0 + (1./np.log(len(wk))*sum(wk*np.log(wk)))
						plt.plot(a[i][0], L_1, marker='+', color='g')
						plt.plot(a[i][1], L_1, marker='+', color='r')
						
						plt.figure(6)
						V = Model(a[i][0], a[i][1]).curvature(z2, a[i][0],a[i][1], 0.0)[1]
						plt.plot(z2, V)
						"""
					return aa
					
				Cov_All_Models = func(All_models,redshift, Beta)	

				#WeightingNormalisation_All_Models = func(All_models,redshift, Beta)
				#Cov_All_Models = [(1./(sqrt(sum(b[0]*(b[1] - sum(b[0]*b[1]))**2.0))*sqrt(sum(b[0]*(m_app - sum(b[0]*m_app))**2.0))))*sum(b[0]*(b[1] - sum(b[0]*b[1]))*(m_app - sum(b[0]*m_app))) for b in WeightingNormalisation_All_Models]	# (b[0], b[1]) = (w_k, Mlist_tild)
			
				# compute of mean of all true models
				#MeanCov_All_Models = np.mean(Cov_All_Models)
				MeanCov_All_Models = Cov_s
				
				#print "Cov_All_Models = ", Cov_All_Models, "\n"
				print "len(Cov_All_Models) = ", len(Cov_All_Models), "\n"
				print "max(Cov_All_Models) = ", max(Cov_All_Models), "\n"
				print "min(Cov_All_Models) = ", min(Cov_All_Models), "\n"
				print "MeanCov_All_Models = ", MeanCov_All_Models, "\n"
				
				print "len(All_models) = ", len(All_models), "\n"
				
				All_models = np.asarray(All_models)
				Cov_All_Models = np.asarray(Cov_All_Models)
				
				Cov_All_Models = Cov_All_Models - MeanCov_All_Models
				
				#***********grid of covariance**********************
				# figure of 3D(surface)
				
				#All_models[:,0] is a liste of lamda
				#All_models[:,1] is a liste of omega
				
				fig = plt.figure(7)
				ax = fig.add_subplot(111)

				col = plt.scatter(All_models[:,1], All_models[:,0], marker='.', s=150, c=Cov_All_Models, linewidths=0.3, cmap=plt.cm.Spectral_r)
				#for vv, ww, dd in zip(All_models[:,1], All_models[:,0], Cov_All_Models):
				#	pl.text(vv, ww, '%.3f' % dd, ha='center', va='bottom')
				

				# Add a colourbar.
				#cax = fig.colorbar(col, orientation='vertical',format='%.30f')
				cax = fig.colorbar(col, orientation='vertical',format='%.4f')
				cax.set_label('Covariance')
				
				plt.plot([bbeta],[alfa], marker='o', color='m')
				plt.xlabel('$\\beta$', fontsize=14)
				plt.ylabel('$\\alpha$', fontsize=14)
				
				#===========For the save of the data====================
				AlfaOfGrid = All_models[:,1]*1.0
				BetaOfGrid = All_models[:,0]*1.0
				CorrCoefficient = Cov_All_Models*1.0
				#=======================================================

				
				#plt.figure(8)
				#for i in range(len(Cov_All_Models)):
				#	plt.plot(Cov_All_Models[i], All_models[i][0], marker='o', color='r')
				
				#n, bins, patches = plt.hist(Cov_All_Models,40,normed='True')
				#(mu_hist, sigma_hist) = norm.fit(array(Cov_All_Models))
				#y = mlab.normpdf( bins, mu_hist, sigma_hist)
				#plt.plot(bins, y, 'r--', linewidth=2)
				#plt.title(r'$\rho_{0}=%.7f,\ \sigma_{0}=%.3f$' %(mu_hist, sigma_hist), color='r')
				#plt.xlabel('Correlation coefficient', fontsize=14)
				
				#***************************************************************************************
						
				Covariance_Vert = []
				Covariance_Hor = []
				
				bbetaCovariance_Vert = []
				alfaCovariance_Vert = []

				bbetaCovariance_Hor = []
				alfaCovariance_Hor = []
				
				
				# make the array as matrix to search the null correlation curve on the grid!
				Cov_All_Models = Cov_All_Models.reshape(sqrt(len(All_models[:,1])), sqrt(len(All_models[:,0])))
				
				All_modelsalfa = All_models[:,0].reshape(sqrt(len(All_models[:,0])), sqrt(len(All_models[:,0]))) # For matrix of lambda
				All_modelsbbeta = All_models[:,1].reshape(sqrt(len(All_models[:,1])), sqrt(len(All_models[:,1]))) # For matrix of omega

				#============**** To make the contours of confidence levels ****=============  
				#CS = pl.contour(All_modelsomega, All_modelslambda, abs(Cov_All_Models), 3, linewidths=np.arange(.5, 4, .5), colors=('r', 'green', 'blue', (1,1,0), '#afeeee', '0.5'))
				#plt.clabel(CS, inline=1, fontsize=10)
				#============***************************************************=============  
				#plt.show()
				#search of covariance zero by HORIZONTAL interpolation
				for i in range(len(bbetalist)):

					#----------------------------------------------------
					
					Condition_1 = any(Cov_All_Models[i] == 0.0)
					if Condition_1 == True:
						rr = Cov_All_Models[i]
						bb = All_modelsbbeta[i] # j'insiste sur omegalist parce que les lignes dans la matrice lambdalist sont toujours divisé avec des valeurs égaux. (interpolation horizontale de la grille (c.a.d de la matrice))# pareil pour omegalist si je fait aussi l'interpolation vertical.
						
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
						bb = All_modelsbbeta[i]
						
						
					rr = array(rr)
					bb = array(bb)
					
					Condition_2 = all(sign(rr)> 0)
					Condition_3 = all(sign(rr)< 0)
					
					if Condition_1 == False and Condition_2 == False:
						Cov_0 = 0.0
						alfa0 = All_modelsalfa[i][0]# lambda est constant ici / En réalité, toujours lamda va etre le meme dans All_modelslambda[i] parce que depuis le début, j'ai fait les deux boucles en fixant lambda0 et faire tourner la deuxième boucle sur omega0.
						bbeta0 = Functions(0.7, 0.3).Nonlinear_Interpolation(bb, rr, Cov_0)
							
						Covariance_Hor.append(Cov_0)
						bbetaCovariance_Hor.append(bbeta0)
						alfaCovariance_Hor.append(alfa0)
					
					else:
						pass
						
					#------------------------------------------------------

				#search of covariance zero by VERTICAL interpolation
				for i in range(len(alfalist)):
					
					#----------------------------------------------------
					
					Condition_1 = any(Cov_All_Models[:,i] == 0.0)
					if Condition_1 == True:
						rr = Cov_All_Models[:,i]
						bb = All_modelsalfa[:,i] # j'insiste sur omegalist parce que les lignes dans la matrice lambdalist sont toujours divisé avec des valeurs égaux. (interpolation horizontale de la grille (c.a.d de la matrice))# pareil pour omegalist si je fait aussi l'interpolation vertical.
						
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
						bb = All_modelsalfa[:,i]
						
						
					rr = array(rr)
					bb = array(bb)
					
					Condition_2 = all(sign(rr)> 0)
					Condition_3 = all(sign(rr)< 0)
					
					if Condition_1 == False and Condition_2 == False:
						Cov_0 = 0.0
						alfa0 = Functions(0.7, 0.3).Nonlinear_Interpolation(bb, rr, Cov_0) # lambda est constant ici / En réalité, toujours lamda va etre le meme dans All_modelslambda[i] parce que depuis le début, j'ai fait les deux boucles en fixant lambda0 et faire tourner la deuxième boucle sur omega0.
						bbeta0 = All_modelsbbeta[:,i][0]
						
						Covariance_Vert.append(Cov_0)
						bbetaCovariance_Vert.append(bbeta0)
						alfaCovariance_Vert.append(alfa0)
					
					else:
						pass
						
					#------------------------------------------------------
					
				
				#plt.figure(9)
				
				Y_Fit, X_Fit = [], []
				alfaCovariance = alfaCovariance_Hor + alfaCovariance_Vert
				bbetaCovariance = bbetaCovariance_Hor + bbetaCovariance_Vert
				#"""
				alfaCovariance1 = list(array(alfaCovariance)*1.0)
				bbetaCovariance1 = list(array(bbetaCovariance)*1.0)
				
				print "alfaCovariance1 = ", alfaCovariance1
				print "bbetaCovariance1 = ", bbetaCovariance1
				
				for i in range(len(alfaCovariance_Vert)):
					plt.plot(bbetaCovariance_Vert[i], alfaCovariance_Vert[i], marker='.', color='r')
					
				for i in range(len(alfaCovariance_Hor)):
					plt.plot(bbetaCovariance_Hor[i], alfaCovariance_Hor[i], marker='.', color='g')
					
				for i in range(len(alfaCovariance)):	
					X_Fit.append(min(bbetaCovariance1))
					Y_Fit.append(alfaCovariance1[bbetaCovariance1.index(min(bbetaCovariance1))])
					
					alfaCovariance1.remove(alfaCovariance1[bbetaCovariance1.index(min(bbetaCovariance1))])
					bbetaCovariance1.remove(min(bbetaCovariance1))
				
				Y_Fit2, X_Fit2 = list(array(alfaCovariance)*1.0), list(array(bbetaCovariance)*1.0)
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
				
				plt.plot(X, Y, color=couleur)
				
				plt.text(0.7, 1.32, '($\\beta$, $\\alpha$) = ('+str(bbeta)+', '+str(alfa)+')', color='r')
				
				plt.plot([bbeta],[alfa], marker='o', color='r')
				plt.xlabel('$\\beta$', fontsize=16)
				plt.ylabel('$\\alpha$', fontsize=16)
				plt.title('The null correlation curve in the ($\\beta$, $\\alpha$) plan', color='b')
				grid(True)
						
				plt.xlim(min(bbetalist), max(bbetalist))                                                                                                                                                                                                                                                                                                                                                                                                               
				plt.ylim(min(alfalist), max(alfalist))      

				
				plt.figure(9)
				
				Y_Fit, X_Fit = [], []
				alfaCovariance = alfaCovariance_Hor + alfaCovariance_Vert
				bbetaCovariance = bbetaCovariance_Hor + bbetaCovariance_Vert
				#"""
				alfaCovariance1 = list(array(alfaCovariance)*1.0)
				bbetaCovariance1 = list(array(bbetaCovariance)*1.0)
				
				print "alfaCovariance1 = ", alfaCovariance1
				print "bbetaCovariance1 = ", bbetaCovariance1
				
				for i in range(len(alfaCovariance_Vert)):
					plt.plot(bbetaCovariance_Vert[i], alfaCovariance_Vert[i], marker='.', color='r')
					
				for i in range(len(alfaCovariance_Hor)):
					plt.plot(bbetaCovariance_Hor[i], alfaCovariance_Hor[i], marker='.', color='g')
					
				for i in range(len(alfaCovariance)):	
					X_Fit.append(min(bbetaCovariance1))
					Y_Fit.append(alfaCovariance1[bbetaCovariance1.index(min(bbetaCovariance1))])
					
					alfaCovariance1.remove(alfaCovariance1[bbetaCovariance1.index(min(bbetaCovariance1))])
					bbetaCovariance1.remove(min(bbetaCovariance1))
				
				Y_Fit2, X_Fit2 = list(array(alfaCovariance)*1.0), list(array(bbetaCovariance)*1.0)
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
				
				plt.plot(X, Y, color=couleur)
				
				plt.text(0.7, 1.32, '($\\beta$, $\\alpha$) = ('+str(bbeta)+', '+str(alfa)+')', color='r')
				
				plt.plot([bbeta],[alfa], marker='o', color='r')
				plt.xlabel('$\\beta$', fontsize=16)
				plt.ylabel('$\\alpha$', fontsize=16)
				plt.title('The null correlation curve in the ($\\beta$, $\\alpha$) plan', color='b')
				grid(True)
						
				plt.xlim(min(bbetalist), max(bbetalist))                                                                                                                                                                                                                                                                                                                                                                                                               
				plt.ylim(min(alfalist), max(alfalist))      
				#plt.show()
				
				if KS == 'K-S_OK':
				  
					DD = []
					
				elif KS == 'K-S_NO':
				    
					DD = []
				      
				else:
					print "You must to determine the condition of the Kolmogrov-Smirnov test as 'K-S_OK' or 'K-S_NO'. ", "\n"

				return BetaOfGrid, AlfaOfGrid, CorrCoefficient, X, Y, DD


			elif Variables == 'color_only':
				All_models = []
				# regroupement de deux valeurs (alfa, bbeta) pour représenter un magnitude absolue M
				for i in range(len(alfalist)):
				    for j in range(len(bbetalist)):
					    All_models.append((alfalist[i], bbetalist[j]))

				def func(a, z, bet):     # a = All_models;  z = redhsift;  bet = Beta
					aa=[]
					#beta = bet
					for i in range(len(a)):

						#Mlist_tild = m_app - ApparentMagnitude(z, lambda0_Model, omega0_Model, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, lambda0_Model, omega0_Model, 0,0,0,0,0)
						#Mlist_tild = M_0 - a[i][0]*(slist-1) + a[i][1]*clist

						cc = ( m_app - M_0 + a[i][0]*(slist -1) - ApparentMagnitude(z, lambda0_Model, omega0_Model, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(z, lambda0_Model, omega0_Model, 0,0,0,0,0) )/a[i][1]
						#cc = clist*1.0
						
						#ten_power_c_m = (10**((m_app - bbeta*clist)/5.))**bet
						#ten_power_c_m = (10**((m_app - clist)/5.))**bet
						ten_power_c_m = (10**((-a[i][0]*slist)/5.))**bet
						
						WeightingManifold = ApparentMagnitude(z, lambda0_Model, omega0_Model, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, bet)   #  a[i][0],a[i][1] = (lamda, omega) pour un modèle donné
						
						WeightingManifold = WeightingManifold*ten_power_c_m
						
						WeightingSum = [b/sum(b) for b in WeightingManifold]
						DetermineBeta = [abs(max(b) - min(b)) for b in WeightingSum]
						beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
						#L_1_Manifold = [1.0 + (1./np.log(len(b))*sum(b*np.log(b))) for b in WeightingSum]
						#beta = Beta[L_1_Manifold.index(min(L_1_Manifold))][0]
						
						#ten_power_c_m = (10**((m_app - clist)/5.))**beta
						ten_power_c_m = (10**((-a[i][0]*slist)/5.))**beta
						w_k = (ApparentMagnitude(z, lambda0_Model, omega0_Model, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta)*ten_power_c_m)/sum((ApparentMagnitude(z, lambda0_Model, omega0_Model, 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(z, beta)*ten_power_c_m))

						COVA_riance = sum(w_k*(cc - sum(w_k*cc ))*(m_app - sum(w_k*m_app)))
						#COR_elation = (1./((len(bbeta*clist)/(len(bbeta*clist) - 1.))*sqrt(sum(w_k*(bbeta*clist - sum(w_k*bbeta*clist))**2.0))*sqrt(sum(w_k*(m_app - sum(w_k*m_app))**2.0))))*COVA_riance

						COR_elation = COVA_riance

						aa.append(COR_elation)

						#aa.append((w_k, Mlist_tild))
						
						#~plt.figure(44)
						#z2 = z*1.0
						#z2.sort()
						#ww = ApparentMagnitude(z2, a[i][0],a[i][1], 0.0).Weightingfactor(z2, beta, 0.0)
						#wk = ww/sum(ww)
						#plt.plot(z2, len(z2)*wk)
						#~for hg in range(len(w_k)):
						#~	plt.plot(z[hg], w_k[hg], marker='+')
							
						"""
						plt.figure(5)
						L_1 = 1.0 + (1./np.log(len(wk))*sum(wk*np.log(wk)))
						plt.plot(a[i][0], L_1, marker='+', color='g')
						plt.plot(a[i][1], L_1, marker='+', color='r')
						
						plt.figure(6)
						V = Model(a[i][0], a[i][1]).curvature(z2, a[i][0],a[i][1], 0.0)[1]
						plt.plot(z2, V)
						"""
					return aa

				Cov_All_Models = func(All_models,redshift, Beta)
			
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
				
				All_models = np.asarray(All_models)
				Cov_All_Models = np.asarray(Cov_All_Models)
				
				Cov_All_Models = Cov_All_Models - MeanCov_All_Models
				
				#***********grid of covariance**********************
				# figure of 3D(surface)
				
				#All_models[:,0] is a liste of lamda
				#All_models[:,1] is a liste of omega
				
				fig = plt.figure(7)
				ax = fig.add_subplot(111)

				col = plt.scatter(All_models[:,1], All_models[:,0], marker='.', s=150, c=Cov_All_Models, linewidths=0.3, cmap=plt.cm.Spectral_r)
				#for vv, ww, dd in zip(All_models[:,1], All_models[:,0], Cov_All_Models):
				#	pl.text(vv, ww, '%.3f' % dd, ha='center', va='bottom')
				

				# Add a colourbar.
				#cax = fig.colorbar(col, orientation='vertical',format='%.30f')
				cax = fig.colorbar(col, orientation='vertical',format='%.4f')
				cax.set_label('Covariance')
				
				plt.plot([bbeta],[alfa], marker='o', color='m')
				plt.xlabel('$\\beta$', fontsize=14)
				plt.ylabel('$\\alpha$', fontsize=14)
				
				#===========For the save of the data====================
				AlfaOfGrid = All_models[:,1]*1.0
				BetaOfGrid = All_models[:,0]*1.0
				CorrCoefficient = Cov_All_Models*1.0
				#=======================================================

				
				#plt.figure(8)
				#for i in range(len(Cov_All_Models)):
				#	plt.plot(Cov_All_Models[i], All_models[i][0], marker='o', color='r')
				
				#n, bins, patches = plt.hist(Cov_All_Models,40,normed='True')
				#(mu_hist, sigma_hist) = norm.fit(array(Cov_All_Models))
				#y = mlab.normpdf( bins, mu_hist, sigma_hist)
				#plt.plot(bins, y, 'r--', linewidth=2)
				#plt.title(r'$\rho_{0}=%.7f,\ \sigma_{0}=%.3f$' %(mu_hist, sigma_hist), color='r')
				#plt.xlabel('Correlation coefficient', fontsize=14)
				
				#***************************************************************************************
						
				Covariance_Vert = []
				Covariance_Hor = []
				
				bbetaCovariance_Vert = []
				alfaCovariance_Vert = []

				bbetaCovariance_Hor = []
				alfaCovariance_Hor = []
				
				
				# make the array as matrix to search the null correlation curve on the grid!
				Cov_All_Models = Cov_All_Models.reshape(sqrt(len(All_models[:,1])), sqrt(len(All_models[:,0])))
				
				All_modelsalfa = All_models[:,0].reshape(sqrt(len(All_models[:,0])), sqrt(len(All_models[:,0]))) # For matrix of lambda
				All_modelsbbeta = All_models[:,1].reshape(sqrt(len(All_models[:,1])), sqrt(len(All_models[:,1]))) # For matrix of omega

				#============**** To make the contours of confidence levels ****=============  
				#CS = pl.contour(All_modelsomega, All_modelslambda, abs(Cov_All_Models), 3, linewidths=np.arange(.5, 4, .5), colors=('r', 'green', 'blue', (1,1,0), '#afeeee', '0.5'))
				#plt.clabel(CS, inline=1, fontsize=10)
				#============***************************************************=============  
				#plt.show()
				#search of covariance zero by HORIZONTAL interpolation
				for i in range(len(bbetalist)):

					#----------------------------------------------------
					
					Condition_1 = any(Cov_All_Models[i] == 0.0)
					if Condition_1 == True:
						rr = Cov_All_Models[i]
						bb = All_modelsbbeta[i] # j'insiste sur omegalist parce que les lignes dans la matrice lambdalist sont toujours divisé avec des valeurs égaux. (interpolation horizontale de la grille (c.a.d de la matrice))# pareil pour omegalist si je fait aussi l'interpolation vertical.
						
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
						bb = All_modelsbbeta[i]
						
						
					rr = array(rr)
					bb = array(bb)
					
					Condition_2 = all(sign(rr)> 0)
					Condition_3 = all(sign(rr)< 0)
					
					if Condition_1 == False and Condition_2 == False:
						Cov_0 = 0.0
						alfa0 = All_modelsalfa[i][0]# lambda est constant ici / En réalité, toujours lamda va etre le meme dans All_modelslambda[i] parce que depuis le début, j'ai fait les deux boucles en fixant lambda0 et faire tourner la deuxième boucle sur omega0.
						bbeta0 = Functions(0.7, 0.3).Nonlinear_Interpolation(bb, rr, Cov_0)
							
						Covariance_Hor.append(Cov_0)
						bbetaCovariance_Hor.append(bbeta0)
						alfaCovariance_Hor.append(alfa0)
					
					else:
						pass
						
					#------------------------------------------------------

				#search of covariance zero by VERTICAL interpolation
				for i in range(len(alfalist)):
					
					#----------------------------------------------------
					
					Condition_1 = any(Cov_All_Models[:,i] == 0.0)
					if Condition_1 == True:
						rr = Cov_All_Models[:,i]
						bb = All_modelsalfa[:,i] # j'insiste sur omegalist parce que les lignes dans la matrice lambdalist sont toujours divisé avec des valeurs égaux. (interpolation horizontale de la grille (c.a.d de la matrice))# pareil pour omegalist si je fait aussi l'interpolation vertical.
						
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
						bb = All_modelsalfa[:,i]
						
						
					rr = array(rr)
					bb = array(bb)
					
					Condition_2 = all(sign(rr)> 0)
					Condition_3 = all(sign(rr)< 0)
					
					if Condition_1 == False and Condition_2 == False:
						Cov_0 = 0.0
						alfa0 = Functions(0.7, 0.3).Nonlinear_Interpolation(bb, rr, Cov_0) # lambda est constant ici / En réalité, toujours lamda va etre le meme dans All_modelslambda[i] parce que depuis le début, j'ai fait les deux boucles en fixant lambda0 et faire tourner la deuxième boucle sur omega0.
						bbeta0 = All_modelsbbeta[:,i][0]
						
						Covariance_Vert.append(Cov_0)
						bbetaCovariance_Vert.append(bbeta0)
						alfaCovariance_Vert.append(alfa0)
					
					else:
						pass
						
					#------------------------------------------------------
					
				
				#plt.figure(9)
				
				Y_Fit, X_Fit = [], []
				alfaCovariance = alfaCovariance_Hor + alfaCovariance_Vert
				bbetaCovariance = bbetaCovariance_Hor + bbetaCovariance_Vert
				#"""
				alfaCovariance1 = list(array(alfaCovariance)*1.0)
				bbetaCovariance1 = list(array(bbetaCovariance)*1.0)
				
				print "alfaCovariance1 = ", alfaCovariance1
				print "bbetaCovariance1 = ", bbetaCovariance1
				
				for i in range(len(alfaCovariance_Vert)):
					plt.plot(bbetaCovariance_Vert[i], alfaCovariance_Vert[i], marker='.', color='r')
					
				for i in range(len(alfaCovariance_Hor)):
					plt.plot(bbetaCovariance_Hor[i], alfaCovariance_Hor[i], marker='.', color='g')
					
				for i in range(len(alfaCovariance)):	
					X_Fit.append(min(bbetaCovariance1))
					Y_Fit.append(alfaCovariance1[bbetaCovariance1.index(min(bbetaCovariance1))])
					
					alfaCovariance1.remove(alfaCovariance1[bbetaCovariance1.index(min(bbetaCovariance1))])
					bbetaCovariance1.remove(min(bbetaCovariance1))
				
				Y_Fit2, X_Fit2 = list(array(alfaCovariance)*1.0), list(array(bbetaCovariance)*1.0)
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
				
				plt.plot(X, Y, color=couleur)
				
				plt.text(0.7, 1.32, '($\\beta$, $\\alpha$) = ('+str(bbeta)+', '+str(alfa)+')', color='r')
				
				plt.plot([bbeta],[alfa], marker='o', color='r')
				plt.xlabel('$\\beta$', fontsize=16)
				plt.ylabel('$\\alpha$', fontsize=16)
				plt.title('The null correlation curve in the ($\\beta$, $\\alpha$) plan', color='b')
				grid(True)
						
				plt.xlim(min(bbetalist), max(bbetalist))                                                                                                                                                                                                                                                                                                                                                                                                               
				plt.ylim(min(alfalist), max(alfalist))      

				
				plt.figure(9)
				
				Y_Fit, X_Fit = [], []
				alfaCovariance = alfaCovariance_Hor + alfaCovariance_Vert
				bbetaCovariance = bbetaCovariance_Hor + bbetaCovariance_Vert
				#"""
				alfaCovariance1 = list(array(alfaCovariance)*1.0)
				bbetaCovariance1 = list(array(bbetaCovariance)*1.0)
				
				print "alfaCovariance1 = ", alfaCovariance1
				print "bbetaCovariance1 = ", bbetaCovariance1
				
				for i in range(len(alfaCovariance_Vert)):
					plt.plot(bbetaCovariance_Vert[i], alfaCovariance_Vert[i], marker='.', color='r')
					
				for i in range(len(alfaCovariance_Hor)):
					plt.plot(bbetaCovariance_Hor[i], alfaCovariance_Hor[i], marker='.', color='g')
					
				for i in range(len(alfaCovariance)):	
					X_Fit.append(min(bbetaCovariance1))
					Y_Fit.append(alfaCovariance1[bbetaCovariance1.index(min(bbetaCovariance1))])
					
					alfaCovariance1.remove(alfaCovariance1[bbetaCovariance1.index(min(bbetaCovariance1))])
					bbetaCovariance1.remove(min(bbetaCovariance1))
				
				Y_Fit2, X_Fit2 = list(array(alfaCovariance)*1.0), list(array(bbetaCovariance)*1.0)
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
				
				plt.plot(X, Y, color=couleur)
				
				plt.text(0.7, 1.32, '($\\beta$, $\\alpha$) = ('+str(bbeta)+', '+str(alfa)+')', color='r')
				
				plt.plot([bbeta],[alfa], marker='o', color='r')
				plt.xlabel('$\\beta$', fontsize=16)
				plt.ylabel('$\\alpha$', fontsize=16)
				plt.title('The null correlation curve in the ($\\beta$, $\\alpha$) plan', color='b')
				grid(True)
						
				plt.xlim(min(bbetalist), max(bbetalist))                                                                                                                                                                                                                                                                                                                                                                                                               
				plt.ylim(min(alfalist), max(alfalist))      
				#plt.show()
				
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
	  
		
