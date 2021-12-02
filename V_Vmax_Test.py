#!/usr/bin/python2.7.3
# -*- coding: latin-1 -*-

# Nom de fichier: V_Vmax_Test.py
#==============================================================================
#title           :V_Vmax_Test.py
#description     :This file is to applicate the V/Vmax Test.
#author          :Dyaa Chbib
#date            :2015_10_19
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
import socket

import Methods
from Methods import Functions
from Cosmological_Model import Model
from HubbleDiagram import ApparentMagnitude
from Savedata import Save

H_0 = np.float64(70.0) #Km.s-1.Mpc-1
LightVelocity = np.float64(2.99792458*10**5.0) # Km.s-1

class V_Vmax:
	"""Creating of the V_Vmax"""
	def __init__(self, lambda_0, omega_0):
	  
		self.param = Functions(lambda_0, omega_0)

	def Search_of_best_Parameters(self, zz, zform, m_lim, m_app, redshift, LambdaOfCorrNull, OmegaOfCorrNull, Realdata, Subsample, Limit_of_Subsample, m_lim_min, m_lim_max, z_min_ofsubsample, z_max_ofsubsample, degree_of_Polyn, Ncoefficients_of_KcorrPolyn, degree_of_Polyn_l, Ncoefficients_of_KcorrPolyn_l):
	
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
	
		
		zzz = np.linspace(10**-15.95, zform, 10**4.)
		#zzz = np.linspace(10**-15.95, zform, len(zz)*10)
		Ncoefficients_of_KcorrPolyn_lll = array([Fit_of_Kcorrection(List_of_Poly, redsh) for redsh in zzz])
		
		DD = []
		COV = []
		List_of_VV = []

		for i in range(len(LambdaOfCorrNull)):
			kappaOfCorrNull = Functions(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).kappa_0
			Mlist_tild = m_app - ApparentMagnitude(redshift, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).m_theor(redshift, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
			#============determine of beta (constant used in the weighting factor term)=======================
			"""
			Betalist = linspace(np.float64(0.01), np.float64(5), 500)
			Beta = array([array([a]) for a in Betalist])
			WeightingManifold = ApparentMagnitude(redshift, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, Beta)
			
			WeightingSum = [a/sum(a) for a in WeightingManifold]
			DetermineBeta = [abs(max(b) - min(b)) for b in WeightingSum]
			beta = Beta[DetermineBeta.index(min(DetermineBeta))][0]
			Weighting = ApparentMagnitude(redshift, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn).Weightingfactor(redshift, beta)

			SumWeighting = sum(Weighting)
			factorWeighting = array(Weighting)/SumWeighting
			"""
			#===================================================================================================
			ZETA = ApparentMagnitude(zz, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, degree_of_Polyn_l, Ncoefficients_of_KcorrPolyn_l).m_theor(zz, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
			vol_l = []
			vol = []
			
			if Subsample == 'Subsample_OK':
				if Limit_of_Subsample == 'mappLimit_OK_redshiftLimit_NO':
					if kappaOfCorrNull > 0:

						ZETA_long = ApparentMagnitude(zzz, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, 0, Ncoefficients_of_KcorrPolyn_lll).m_theor(zzz, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
						
						#Angular_distance = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zzz, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[0]*sqrt(kappaOfCorrNull)
						#z_star = Functions(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).Nonlinear_Interpolation(zzz, Angular_distance, pi/2)
						#Ncoefficients_of_KcorrPolyn_z_star = Fit_of_Kcorrection(List_of_Poly,z_star)
						#degree_of_Polyn_z_star = 0 # pas importante tant que je ne calcul pas le W_k
						#ZETA_star = ApparentMagnitude(z_star, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, degree_of_Polyn_z_star, Ncoefficients_of_KcorrPolyn_z_star).m_theor(z_star, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
						#M_star = m_lim - ZETA_star
						z_star = Functions(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).Nonlinear_Interpolation(zzz, ZETA_long, max(ZETA_long))
						#--------Pour m_lim_min--------
						M_star_m_lim_min = m_lim_min - max(ZETA_long)
						Ncoefficients_of_KcorrPolyn_zform = Fit_of_Kcorrection(List_of_Poly,zform)
						degree_of_Polyn_zform = 0 # pas importante tant que je ne calcul pas le W_k
						ZETA_form = ApparentMagnitude(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, degree_of_Polyn_zform, Ncoefficients_of_KcorrPolyn_zform).m_theor(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
						M_form_m_lim_min = m_lim_min - ZETA_form
						#------------------------------
						#--------Pour m_lim_max--------
						M_star_m_lim_max = m_lim_max - max(ZETA_long)
						Ncoefficients_of_KcorrPolyn_zform = Fit_of_Kcorrection(List_of_Poly,zform)
						degree_of_Polyn_zform = 0 # pas importante tant que je ne calcul pas le W_k
						ZETA_form = ApparentMagnitude(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, degree_of_Polyn_zform, Ncoefficients_of_KcorrPolyn_zform).m_theor(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
						M_form_m_lim_max = m_lim_max - ZETA_form
						#------------------------------
						
						for j in range(len(redshift)):
	
							if Mlist_tild[j]>M_star_m_lim_max and Mlist_tild[j]<M_form_m_lim_max:
							  
							#----recherche de valeurs de zl avec m_lim_max (premier limite en magnitude apparente)--------------		  
								mu_lim_max = m_lim_max - Mlist_tild[j]
								
								if Realdata == 'Realdata_OK':
								
									ZETA_jj = ZETA_long*1.0
									zz_jj = zzz*1.0
									z_l_jj = []
									#for jj in range(20): # A utliser quand je cherche les valeurs de z_l avec la k-correction sur le diagramme M-V
									for jj in range(10): # A utliser quand je cherche les valeurs de z_l sans la k-correction sur le diagramme M-V
										z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim_max)
										z_l_jj.append(z_ljj)
										index_of_mu_lim_max = list(abs(ZETA_jj - mu_lim_max)).index(min(abs(ZETA_jj - mu_lim_max)))
										ZETA_jj = list(ZETA_jj*1)
										zz_jj = list(zz_jj*1)
										ZETA_jj.remove(ZETA_jj[index_of_mu_lim_max])
										zz_jj.remove(zz_jj[index_of_mu_lim_max])
										ZETA_jj = array(ZETA_jj*1)
										zz_jj = array(zz_jj*1)
									
									z_l_jj = array(z_l_jj)
									z_l_jj.sort()
									# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
									z_l1 = min(z_l_jj)
									z_l2 = max(z_l_jj)
									z_l_m_lim_max = [z_l1, z_l2]
									
									V_l1_m_lim_max = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_l2_m_lim_max =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									
								else:
									ZETA_jj = ZETA_long*1.0
									zz_jj = zzz*1.0
									z_l_jj = []
									for jj in range(5):
										z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim_max)
										z_l_jj.append(z_ljj)
										index_of_mu_lim_max = list(abs(ZETA_jj - mu_lim_max)).index(min(abs(ZETA_jj - mu_lim_max)))
										ZETA_jj = list(ZETA_jj*1)
										zz_jj = list(zz_jj*1)
										ZETA_jj.remove(ZETA_jj[index_of_mu_lim_max])
										zz_jj.remove(zz_jj[index_of_mu_lim_max])
										ZETA_jj = array(ZETA_jj*1)
										zz_jj = array(zz_jj*1)
									
									z_l_jj = array(z_l_jj)
									z_l_jj.sort()
									# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
									z_l1 = min(z_l_jj)
									z_l2 = max(z_l_jj)
									z_l_m_lim_max = [z_l1, z_l2]
									
									V_l1_m_lim_max = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_l2_m_lim_max =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

							#----recherche de valeurs de zl avec m_lim_min (premier limite en magnitude apparente)--------------
							  
								mu_lim_min = m_lim_min - Mlist_tild[j]
								
								if Realdata == 'Realdata_OK':
								
									ZETA_jj = ZETA_long*1.0
									zz_jj = zzz*1.0
									z_l_jj = []
									#for jj in range(20): # A utliser quand je cherche les valeurs de z_l avec la k-correction sur le diagramme M-V
									for jj in range(10): # A utliser quand je cherche les valeurs de z_l sans la k-correction sur le diagramme M-V
										z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim_min)
										z_l_jj.append(z_ljj)
										index_of_mu_lim_min = list(abs(ZETA_jj - mu_lim_min)).index(min(abs(ZETA_jj - mu_lim_min)))
										ZETA_jj = list(ZETA_jj*1)
										zz_jj = list(zz_jj*1)
										ZETA_jj.remove(ZETA_jj[index_of_mu_lim_min])
										zz_jj.remove(zz_jj[index_of_mu_lim_min])
										ZETA_jj = array(ZETA_jj*1)
										zz_jj = array(zz_jj*1)
									
									z_l_jj = array(z_l_jj)
									z_l_jj.sort()
									# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
									z_l1 = min(z_l_jj)
									z_l2 = max(z_l_jj)
									z_l_m_lim_min = [z_l1, z_l2]
									
									V_l1_m_lim_min = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_l2_m_lim_min =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								else:
									ZETA_jj = ZETA_long*1.0
									zz_jj = zzz*1.0
									z_l_jj = []
									for jj in range(5):
										z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim_min)
										z_l_jj.append(z_ljj)
										index_of_mu_lim_min = list(abs(ZETA_jj - mu_lim_min)).index(min(abs(ZETA_jj - mu_lim_min)))
										ZETA_jj = list(ZETA_jj*1)
										zz_jj = list(zz_jj*1)
										ZETA_jj.remove(ZETA_jj[index_of_mu_lim_min])
										zz_jj.remove(zz_jj[index_of_mu_lim_min])
										ZETA_jj = array(ZETA_jj*1)
										zz_jj = array(zz_jj*1)
									
									z_l_jj = array(z_l_jj)
									z_l_jj.sort()
									# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
									z_l1 = min(z_l_jj)
									z_l2 = max(z_l_jj)
									z_l_m_lim_min = [z_l1, z_l2]
									
									V_l1_m_lim_min = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_l2_m_lim_min =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									
								V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_zform = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

								"""
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
										condition2 = (abs(float(Distanguish_of_solutions[i]) - float(Distanguish_of_solutions[i-1]))>= 0.1)
										if condition1 == False and condition2 == True:
											z_l.append(z_l_l[i])
										else:
											pass
									
									if len(z_l)==4:
										
										z_1, z_2, z_3, z_4 = z_l[0], z_l[1], z_l[2], z_l[3]
										
										V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_l1 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_1, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
										V_l2 =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_2, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_l3 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_3, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
										V_l4 =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_4, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_zform = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

										if redshift[j] <= z_1:
											
											V_0 = V_z
											vol.append(V_0)
											V_l = V_l1 + V_zform - V_l4 + V_l3 - V_l2
											vol_l.append(V_l)

												
										elif redshift[j] >= z_2 and redshift[j] <= z_3:

											V_0 = V_l1 + V_z - V_l2
											vol.append(V_0)
											V_l = V_l1 + V_zform - V_l4 + V_l3 - V_l2
											vol_l.append(V_l)
											
										elif redshift[j] >= z_4:

											V_0 = V_l1 + V_z - V_l4 + V_l3 - V_l2
											vol.append(V_0)
											V_l = V_l1 + V_zform - V_l4 + V_l3 - V_l2
											vol_l.append(V_l)
											
										else:
											print "le redshift de l'objet est entre z_l1 et z_l2 ou entre z_l3 et z_l4", "\n"
											print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
											
									elif len(z_l)==2:
										z_l1 = min(z_l_jj)
										z_l2 = max(z_l_jj)
										z_l = [z_l1, z_l2]
										
										V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_l1 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
										V_l2 =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_zform = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

										if redshift[j] <= z_l[0]:
											
											V_0 = V_z
											vol.append(V_0)
											V_l = V_l1 + V_zform - V_l2
											vol_l.append(V_l)

												
										elif redshift[j] >= z_l[1]:
											V_0 = V_l1 + V_z - V_l2
											vol.append(V_0)
											V_l = V_l1 + V_zform - V_l2
											vol_l.append(V_l)
										else:
											print "le redshift de l'objet est entre z_l1 et z_l2 ou entre z_l3 et z_l4", "\n"
											print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
									else:
										print "len(z_l) n'est pas égale à 2 ni à 4, len(z_l) = ", len(z_l), ", z_l = ", z_l, "\n"
										print "z_l_l = ", z_l_l, "\n"
										print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
										
								"""	

								if redshift[j] <= z_l_m_lim_max[0]:
									
									V_0 = V_z - V_l1_m_lim_min
									vol.append(V_0)
									V_l = V_l1_m_lim_max - V_l1_m_lim_min
									vol_l.append(V_l)

								
								elif redshift[j] >= z_l_m_lim_max[0]:
									V_0 = V_l1_m_lim_max - V_l1_m_lim_min + V_z - V_l2_m_lim_max
									vol.append(V_0)
									V_l = V_l1_m_lim_max - V_l1_m_lim_min + V_l2_m_lim_min - V_l2_m_lim_max
									vol_l.append(V_l)
								else:
									print "le redshift de l'objet est entre z_l_m_lim_max[0] et z_l_m_lim_max[1]", "\n"
									print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
								
							else:
								mu_lim_min = m_lim_min - Mlist_tild[j]
								mu_lim_max = m_lim_max - Mlist_tild[j]
								if Realdata == 'Realdata_OK':
									# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
									z_l_m_lim_min = Functions(0.7, 0.3).Nonlinear_Interpolation(zzz, ZETA_long, mu_lim_min) 
									z_l_m_lim_max = Functions(0.7, 0.3).Nonlinear_Interpolation(zzz, ZETA_long, mu_lim_max)
									V_min_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_max_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_0 = V_z - V_min_lim
									vol.append(V_0)
									V_l = V_max_lim - V_min_lim
									vol_l.append(V_l)
									
									"""
									ZETA_jj = ZETA_long*1.0
									zz_jj = zzz*1.0
									z_l_jj = []
									for jj in range(10):
										z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim)
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
										condition2 = (abs(float(Distanguish_of_solutions[i]) - float(Distanguish_of_solutions[i-1]))>= 0.1)
										if condition1 == False and condition2 == True:
											z_l.append(z_l_l[i])
										else:
											pass
									if len(z_l)==3:
										
										z_1, z_2, z_3 = z_l[0], z_l[1], z_l[2]
										
										V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_l1 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_1, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_l2 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_2, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_l3 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_3, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

										#print "z_1, z_2, z_3 = ", z_1, z_2, z_3, "\n"
										
										if redshift[j] <= z_1:
										
											V_0 = V_z
											vol.append(V_0)
											V_l = V_l1 + V_l3 - V_l2
											vol_l.append(V_l)

										elif redshift[j] >= z_2:
											V_0 = V_l1 + V_z - V_l2
											vol.append(V_0)
											V_l = V_l1 + V_l3 - V_l2
											vol_l.append(V_l)
										else:
											print "le redshift de l'objet est entre z_l1 et z_l2 de K-correction", "\n"
											print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
											
									elif len(z_l)==1:
										z_l = z_l[0]
										V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_0 = V_z
										vol.append(V_0)
										V_l =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										vol_l.append(V_l)
									else:
										print "len(z_l) n'est pas égale à 1 ni à 3, len(z_l) = ", len(z_l), ", z_l = ", z_l, "\n"
										print "z_l_l = ", z_l_l, "\n"
										print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
									"""
								else:
									z_l_m_lim_min = Functions(0.7, 0.3).Nonlinear_Interpolation(zzz, ZETA_long, mu_lim_min) 
									z_l_m_lim_max = Functions(0.7, 0.3).Nonlinear_Interpolation(zzz, ZETA_long, mu_lim_max)
									V_min_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_max_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_0 = V_z - V_min_lim
									vol.append(V_0)
									V_l = V_max_lim - V_min_lim
									vol_l.append(V_l)

					else:
						for j in range(len(redshift)):
						
							mu_lim_min = m_lim_min - Mlist_tild[j]
							mu_lim_max = m_lim_max - Mlist_tild[j]
							z_l_m_lim_min = Functions(0.7, 0.3).Nonlinear_Interpolation(zz, ZETA, mu_lim_min) 
							z_l_m_lim_max = Functions(0.7, 0.3).Nonlinear_Interpolation(zz, ZETA, mu_lim_max)
							V_min_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
							V_max_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
							V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
							V_0 = V_z - V_min_lim
							vol.append(V_0)
							V_l = V_max_lim - V_min_lim
							vol_l.append(V_l)
							
				elif Limit_of_Subsample == 'mappLimit_OK_redshiftLimit_OK':
					if kappaOfCorrNull > 0:

						ZETA_long = ApparentMagnitude(zzz, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, 0, Ncoefficients_of_KcorrPolyn_lll).m_theor(zzz, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
						
						#Angular_distance = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zzz, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[0]*sqrt(kappaOfCorrNull)
						#z_star = Functions(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).Nonlinear_Interpolation(zzz, Angular_distance, pi/2)
						#Ncoefficients_of_KcorrPolyn_z_star = Fit_of_Kcorrection(List_of_Poly,z_star)
						#degree_of_Polyn_z_star = 0 # pas importante tant que je ne calcul pas le W_k
						#ZETA_star = ApparentMagnitude(z_star, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, degree_of_Polyn_z_star, Ncoefficients_of_KcorrPolyn_z_star).m_theor(z_star, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
						#M_star = m_lim - ZETA_star
						z_star = Functions(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).Nonlinear_Interpolation(zzz, ZETA_long, max(ZETA_long))
						#--------Pour m_lim_min--------
						M_star_m_lim_min = m_lim_min - max(ZETA_long)
						Ncoefficients_of_KcorrPolyn_zform = Fit_of_Kcorrection(List_of_Poly,zform)
						degree_of_Polyn_zform = 0 # pas importante tant que je ne calcul pas le W_k
						ZETA_form = ApparentMagnitude(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, degree_of_Polyn_zform, Ncoefficients_of_KcorrPolyn_zform).m_theor(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
						M_form_m_lim_min = m_lim_min - ZETA_form
						#------------------------------
						#--------Pour m_lim_max--------
						M_star_m_lim_max = m_lim_max - max(ZETA_long)
						Ncoefficients_of_KcorrPolyn_zform = Fit_of_Kcorrection(List_of_Poly,zform)
						degree_of_Polyn_zform = 0 # pas importante tant que je ne calcul pas le W_k
						ZETA_form = ApparentMagnitude(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, degree_of_Polyn_zform, Ncoefficients_of_KcorrPolyn_zform).m_theor(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
						M_form_m_lim_max = m_lim_max - ZETA_form
						#------------------------------
						
						for j in range(len(redshift)):
	
							if Mlist_tild[j]>M_star_m_lim_max and Mlist_tild[j]<M_form_m_lim_max:
							  
							#----recherche de valeurs de zl avec m_lim_max (premier limite en magnitude apparente)--------------		  
								mu_lim_max = m_lim_max - Mlist_tild[j]
								
								if Realdata == 'Realdata_OK':
								
									ZETA_jj = ZETA_long*1.0
									zz_jj = zzz*1.0
									z_l_jj = []
									#for jj in range(20): # A utliser quand je cherche les valeurs de z_l avec la k-correction sur le diagramme M-V
									for jj in range(10): # A utliser quand je cherche les valeurs de z_l sans la k-correction sur le diagramme M-V
										z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim_max)
										z_l_jj.append(z_ljj)
										index_of_mu_lim_max = list(abs(ZETA_jj - mu_lim_max)).index(min(abs(ZETA_jj - mu_lim_max)))
										ZETA_jj = list(ZETA_jj*1)
										zz_jj = list(zz_jj*1)
										ZETA_jj.remove(ZETA_jj[index_of_mu_lim_max])
										zz_jj.remove(zz_jj[index_of_mu_lim_max])
										ZETA_jj = array(ZETA_jj*1)
										zz_jj = array(zz_jj*1)
									
									z_l_jj = array(z_l_jj)
									z_l_jj.sort()
									# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
									z_l1 = min(z_l_jj)
									z_l2 = max(z_l_jj)
									z_l_m_lim_max = [z_l1, z_l2]
									
									V_l1_m_lim_max = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_l2_m_lim_max =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									
								else:
									ZETA_jj = ZETA_long*1.0
									zz_jj = zzz*1.0
									z_l_jj = []
									for jj in range(5):
										z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim_max)
										z_l_jj.append(z_ljj)
										index_of_mu_lim_max = list(abs(ZETA_jj - mu_lim_max)).index(min(abs(ZETA_jj - mu_lim_max)))
										ZETA_jj = list(ZETA_jj*1)
										zz_jj = list(zz_jj*1)
										ZETA_jj.remove(ZETA_jj[index_of_mu_lim_max])
										zz_jj.remove(zz_jj[index_of_mu_lim_max])
										ZETA_jj = array(ZETA_jj*1)
										zz_jj = array(zz_jj*1)
									
									z_l_jj = array(z_l_jj)
									z_l_jj.sort()
									# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
									z_l1 = min(z_l_jj)
									z_l2 = max(z_l_jj)
									z_l_m_lim_max = [z_l1, z_l2]
									
									V_l1_m_lim_max = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_l2_m_lim_max =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

							#----recherche de valeurs de zl avec m_lim_min (premier limite en magnitude apparente)--------------
							  
								mu_lim_min = m_lim_min - Mlist_tild[j]
								
								if Realdata == 'Realdata_OK':
								
									ZETA_jj = ZETA_long*1.0
									zz_jj = zzz*1.0
									z_l_jj = []
									#for jj in range(20): # A utliser quand je cherche les valeurs de z_l avec la k-correction sur le diagramme M-V
									for jj in range(10): # A utliser quand je cherche les valeurs de z_l sans la k-correction sur le diagramme M-V
										z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim_min)
										z_l_jj.append(z_ljj)
										index_of_mu_lim_min = list(abs(ZETA_jj - mu_lim_min)).index(min(abs(ZETA_jj - mu_lim_min)))
										ZETA_jj = list(ZETA_jj*1)
										zz_jj = list(zz_jj*1)
										ZETA_jj.remove(ZETA_jj[index_of_mu_lim_min])
										zz_jj.remove(zz_jj[index_of_mu_lim_min])
										ZETA_jj = array(ZETA_jj*1)
										zz_jj = array(zz_jj*1)
									
									z_l_jj = array(z_l_jj)
									z_l_jj.sort()
									# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
									z_l1 = min(z_l_jj)
									z_l2 = max(z_l_jj)
									z_l_m_lim_min = [z_l1, z_l2]
									
									V_l1_m_lim_min = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_l2_m_lim_min =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								else:
									ZETA_jj = ZETA_long*1.0
									zz_jj = zzz*1.0
									z_l_jj = []
									for jj in range(5):
										z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim_min)
										z_l_jj.append(z_ljj)
										index_of_mu_lim_min = list(abs(ZETA_jj - mu_lim_min)).index(min(abs(ZETA_jj - mu_lim_min)))
										ZETA_jj = list(ZETA_jj*1)
										zz_jj = list(zz_jj*1)
										ZETA_jj.remove(ZETA_jj[index_of_mu_lim_min])
										zz_jj.remove(zz_jj[index_of_mu_lim_min])
										ZETA_jj = array(ZETA_jj*1)
										zz_jj = array(zz_jj*1)
									
									z_l_jj = array(z_l_jj)
									z_l_jj.sort()
									# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
									z_l1 = min(z_l_jj)
									z_l2 = max(z_l_jj)
									z_l_m_lim_min = [z_l1, z_l2]
									
									V_l1_m_lim_min = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_l2_m_lim_min =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									
								V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_zform = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_max = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_max_ofsubsample, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_min = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_min_ofsubsample, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

								"""
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
										condition2 = (abs(float(Distanguish_of_solutions[i]) - float(Distanguish_of_solutions[i-1]))>= 0.1)
										if condition1 == False and condition2 == True:
											z_l.append(z_l_l[i])
										else:
											pass
									
									if len(z_l)==4:
										
										z_1, z_2, z_3, z_4 = z_l[0], z_l[1], z_l[2], z_l[3]
										
										V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_l1 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_1, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
										V_l2 =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_2, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_l3 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_3, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
										V_l4 =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_4, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_zform = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

										if redshift[j] <= z_1:
											
											V_0 = V_z
											vol.append(V_0)
											V_l = V_l1 + V_zform - V_l4 + V_l3 - V_l2
											vol_l.append(V_l)

												
										elif redshift[j] >= z_2 and redshift[j] <= z_3:

											V_0 = V_l1 + V_z - V_l2
											vol.append(V_0)
											V_l = V_l1 + V_zform - V_l4 + V_l3 - V_l2
											vol_l.append(V_l)
											
										elif redshift[j] >= z_4:

											V_0 = V_l1 + V_z - V_l4 + V_l3 - V_l2
											vol.append(V_0)
											V_l = V_l1 + V_zform - V_l4 + V_l3 - V_l2
											vol_l.append(V_l)
											
										else:
											print "le redshift de l'objet est entre z_l1 et z_l2 ou entre z_l3 et z_l4", "\n"
											print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
											
									elif len(z_l)==2:
										z_l1 = min(z_l_jj)
										z_l2 = max(z_l_jj)
										z_l = [z_l1, z_l2]
										
										V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_l1 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
										V_l2 =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_zform = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

										if redshift[j] <= z_l[0]:
											
											V_0 = V_z
											vol.append(V_0)
											V_l = V_l1 + V_zform - V_l2
											vol_l.append(V_l)

												
										elif redshift[j] >= z_l[1]:
											V_0 = V_l1 + V_z - V_l2
											vol.append(V_0)
											V_l = V_l1 + V_zform - V_l2
											vol_l.append(V_l)
										else:
											print "le redshift de l'objet est entre z_l1 et z_l2 ou entre z_l3 et z_l4", "\n"
											print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
									else:
										print "len(z_l) n'est pas égale à 2 ni à 4, len(z_l) = ", len(z_l), ", z_l = ", z_l, "\n"
										print "z_l_l = ", z_l_l, "\n"
										print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
										
								"""

								if redshift[j] <= z_l_m_lim_max[0]:
									
									V_0 = V_z - max(V_min, V_l1_m_lim_min)
									vol.append(V_0)
									V_l = min(V_max, V_l1_m_lim_max) - max(V_min, V_l1_m_lim_min)
									vol_l.append(V_l)

								
								elif redshift[j] >= z_l_m_lim_max[0]:
									V_0 = min(V_max, V_l1_m_lim_max) - max(V_min, V_l1_m_lim_min) + V_z - max(V_min, V_l2_m_lim_max)
									vol.append(V_0)
									V_l = min(V_max, V_l1_m_lim_max) - max(V_min, V_l1_m_lim_min) + min(V_max, V_l2_m_lim_max) - max(V_min, V_l2_m_lim_max)
									vol_l.append(V_l)
								else:
									print "le redshift de l'objet est entre z_l_m_lim_max[0] et z_l_m_lim_max[1]", "\n"
									print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
								
							else:
								mu_lim_min = m_lim_min - Mlist_tild[j]
								mu_lim_max = m_lim_max - Mlist_tild[j]
								V_max = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_max_ofsubsample, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_min = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_min_ofsubsample, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

								if Realdata == 'Realdata_OK':
									# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
									z_l_m_lim_min = Functions(0.7, 0.3).Nonlinear_Interpolation(zz, ZETA, mu_lim_min) 
									z_l_m_lim_max = Functions(0.7, 0.3).Nonlinear_Interpolation(zz, ZETA, mu_lim_max)
									V_min_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_max_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_0 = V_z - max(V_min, V_min_lim)
									vol.append(V_0)
									V_l = min(V_max, V_max_lim) - max(V_min, V_min_lim)
									vol_l.append(V_l)
									
									"""
									ZETA_jj = ZETA_long*1.0
									zz_jj = zzz*1.0
									z_l_jj = []
									for jj in range(10):
										z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim)
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
										condition2 = (abs(float(Distanguish_of_solutions[i]) - float(Distanguish_of_solutions[i-1]))>= 0.1)
										if condition1 == False and condition2 == True:
											z_l.append(z_l_l[i])
										else:
											pass
									if len(z_l)==3:
										
										z_1, z_2, z_3 = z_l[0], z_l[1], z_l[2]
										
										V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_l1 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_1, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_l2 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_2, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_l3 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_3, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

										#print "z_1, z_2, z_3 = ", z_1, z_2, z_3, "\n"
										
										if redshift[j] <= z_1:
										
											V_0 = V_z
											vol.append(V_0)
											V_l = V_l1 + V_l3 - V_l2
											vol_l.append(V_l)

										elif redshift[j] >= z_2:
											V_0 = V_l1 + V_z - V_l2
											vol.append(V_0)
											V_l = V_l1 + V_l3 - V_l2
											vol_l.append(V_l)
										else:
											print "le redshift de l'objet est entre z_l1 et z_l2 de K-correction", "\n"
											print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
											
									elif len(z_l)==1:
										z_l = z_l[0]
										V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										V_0 = V_z
										vol.append(V_0)
										V_l =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
										vol_l.append(V_l)
									else:
										print "len(z_l) n'est pas égale à 1 ni à 3, len(z_l) = ", len(z_l), ", z_l = ", z_l, "\n"
										print "z_l_l = ", z_l_l, "\n"
										print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
									"""
								else:
									z_l_m_lim_min = Functions(0.7, 0.3).Nonlinear_Interpolation(zz, ZETA, mu_lim_min) 
									z_l_m_lim_max = Functions(0.7, 0.3).Nonlinear_Interpolation(zz, ZETA, mu_lim_max)
									V_min_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_max_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_0 = V_z - max(V_min, V_min_lim)
									vol.append(V_0)
									V_l = min(V_max, V_max_lim) - max(V_min, V_min_lim)
									vol_l.append(V_l)

					else:
						for j in range(len(redshift)):
						
							mu_lim_min = m_lim_min - Mlist_tild[j]
							mu_lim_max = m_lim_max - Mlist_tild[j]
							V_max = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_max_ofsubsample, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
							V_min = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_min_ofsubsample, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

							if Realdata == 'Realdata_OK':
								# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
								z_l_m_lim_min = Functions(0.7, 0.3).Nonlinear_Interpolation(zz, ZETA, mu_lim_min) 
								z_l_m_lim_max = Functions(0.7, 0.3).Nonlinear_Interpolation(zz, ZETA, mu_lim_max)
								V_min_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
								V_max_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
								V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_0 = V_z - max(V_min, V_min_lim)
								vol.append(V_0)
								V_l = min(V_max, V_max_lim) - max(V_min, V_min_lim)
								vol_l.append(V_l)
									
							else:
								z_l_m_lim_min = Functions(0.7, 0.3).Nonlinear_Interpolation(zz, ZETA, mu_lim_min) 
								z_l_m_lim_max = Functions(0.7, 0.3).Nonlinear_Interpolation(zz, ZETA, mu_lim_max)
								V_min_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_min, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
								V_max_lim = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l_m_lim_max, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
								V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_0 = V_z - max(V_min, V_min_lim)
								vol.append(V_0)
								V_l = min(V_max, V_max_lim) - max(V_min, V_min_lim)
								vol_l.append(V_l)

				
				else:
					print "You must to determine the condition of the Subsample as 'mappLimit_OK_redshiftLimit_OK' or 'mappLimit_OK_redshiftLimit_NO'. ", "\n"
					
					
			elif Subsample == 'Subsample_NO':
				if kappaOfCorrNull > 0:
					
					ZETA_long = ApparentMagnitude(zzz, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, 0, Ncoefficients_of_KcorrPolyn_lll).m_theor(zzz, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
					
					#Angular_distance = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zzz, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[0]*sqrt(kappaOfCorrNull)
					#z_star = Functions(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).Nonlinear_Interpolation(zzz, Angular_distance, pi/2)
					#Ncoefficients_of_KcorrPolyn_z_star = Fit_of_Kcorrection(List_of_Poly,z_star)
					#degree_of_Polyn_z_star = 0 # pas importante tant que je ne calcul pas le W_k
					#ZETA_star = ApparentMagnitude(z_star, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, degree_of_Polyn_z_star, Ncoefficients_of_KcorrPolyn_z_star).m_theor(z_star, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
					#M_star = m_lim - ZETA_star
					z_star = Functions(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).Nonlinear_Interpolation(zzz, ZETA_long, max(ZETA_long))
					M_star = m_lim - max(ZETA_long)
					Ncoefficients_of_KcorrPolyn_zform = Fit_of_Kcorrection(List_of_Poly,zform)
					degree_of_Polyn_zform = 0 # pas importante tant que je ne calcul pas le W_k
					ZETA_form = ApparentMagnitude(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0, Realdata, degree_of_Polyn_zform, Ncoefficients_of_KcorrPolyn_zform).m_theor(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)
					M_form = m_lim - ZETA_form
					
					for j in range(len(redshift)):

						if Mlist_tild[j]>M_star and Mlist_tild[j]<M_form:
							
							#print "je suis entre M_star et M_form", "\n"
							#print "Mlist_tild[j] = ", Mlist_tild[j], "\n"
							#print "M_star, M_form = ", M_star, M_form, "\n"
							#print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
							mu_lim = m_lim - Mlist_tild[j]
							if Realdata == 'Realdata_OK':
								
								ZETA_jj = ZETA_long*1.0
								zz_jj = zzz*1.0
								z_l_jj = []
								#for jj in range(20): # A utliser quand je cherche les valeurs de z_l avec la k-correction sur le diagramme M-V
								for jj in range(10): # A utliser quand je cherche les valeurs de z_l sans la k-correction sur le diagramme M-V
									z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim)
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
								# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
								z_l1 = min(z_l_jj)
								z_l2 = max(z_l_jj)
								z_l = [z_l1, z_l2]
								
								V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_l1 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
								V_l2 =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_zform = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

								if redshift[j] <= z_l[0]:
									
									V_0 = V_z
									vol.append(V_0)
									V_l = V_l1 + V_zform - V_l2
									vol_l.append(V_l)

										
								elif redshift[j] >= z_l[1]:
									V_0 = V_l1 + V_z - V_l2
									vol.append(V_0)
									V_l = V_l1 + V_zform - V_l2
									vol_l.append(V_l)
								else:
									print "le redshift de l'objet est entre z_l1 et z_l2", "\n"
									print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
								
								"""
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
									condition2 = (abs(float(Distanguish_of_solutions[i]) - float(Distanguish_of_solutions[i-1]))>= 0.1)
									if condition1 == False and condition2 == True:
										z_l.append(z_l_l[i])
									else:
										pass
								
								if len(z_l)==4:
									
									z_1, z_2, z_3, z_4 = z_l[0], z_l[1], z_l[2], z_l[3]
									
									V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_l1 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_1, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_l2 =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_2, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_l3 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_3, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_l4 =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_4, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_zform = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

									if redshift[j] <= z_1:
										
										V_0 = V_z
										vol.append(V_0)
										V_l = V_l1 + V_zform - V_l4 + V_l3 - V_l2
										vol_l.append(V_l)

											
									elif redshift[j] >= z_2 and redshift[j] <= z_3:

										V_0 = V_l1 + V_z - V_l2
										vol.append(V_0)
										V_l = V_l1 + V_zform - V_l4 + V_l3 - V_l2
										vol_l.append(V_l)
										
									elif redshift[j] >= z_4:

										V_0 = V_l1 + V_z - V_l4 + V_l3 - V_l2
										vol.append(V_0)
										V_l = V_l1 + V_zform - V_l4 + V_l3 - V_l2
										vol_l.append(V_l)
										
									else:
										print "le redshift de l'objet est entre z_l1 et z_l2 ou entre z_l3 et z_l4", "\n"
										print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
										
								elif len(z_l)==2:
									z_l1 = min(z_l_jj)
									z_l2 = max(z_l_jj)
									z_l = [z_l1, z_l2]
									
									V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_l1 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
									V_l2 =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_zform = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

									if redshift[j] <= z_l[0]:
										
										V_0 = V_z
										vol.append(V_0)
										V_l = V_l1 + V_zform - V_l2
										vol_l.append(V_l)

											
									elif redshift[j] >= z_l[1]:
										V_0 = V_l1 + V_z - V_l2
										vol.append(V_0)
										V_l = V_l1 + V_zform - V_l2
										vol_l.append(V_l)
									else:
										print "le redshift de l'objet est entre z_l1 et z_l2 ou entre z_l3 et z_l4", "\n"
										print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
								else:
									print "len(z_l) n'est pas égale à 2 ni à 4, len(z_l) = ", len(z_l), ", z_l = ", z_l, "\n"
									print "z_l_l = ", z_l_l, "\n"
									print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
									
								"""	
							else:
								ZETA_jj = ZETA_long*1.0
								zz_jj = zzz*1.0
								z_l_jj = []
								for jj in range(5):
									z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim)
									z_l_jj.append(z_ljj)
									index_of_mu_lim = list(abs(ZETA_jj - mu_lim)).index(min(abs(ZETA_jj - mu_lim)))
									ZETA_jj = list(ZETA_jj*1)
									zz_jj = list(zz_jj*1)
									ZETA_jj.remove(ZETA_jj[index_of_mu_lim])
									zz_jj.remove(zz_jj[index_of_mu_lim])
									ZETA_jj = array(ZETA_jj*1)
									zz_jj = array(zz_jj*1)
								
								z_l_jj = array(z_l_jj)
								z_l1 = min(z_l_jj)
								z_l2 = max(z_l_jj)
								z_l = [z_l1, z_l2]
								
								V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_l1 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l[0], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1] 
								V_l2 =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l[1], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_zform = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(zform, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

								if redshift[j] <= z_l[0]:
									
									V_0 = V_z
									vol.append(V_0)
									V_l = V_l1 + V_zform - V_l2
									vol_l.append(V_l)

										
								elif redshift[j] >= z_l[1]:
									V_0 = V_l1 + V_z - V_l2
									vol.append(V_0)
									V_l = V_l1 + V_zform - V_l2
									vol_l.append(V_l)
								else:
									print "le redshift de l'objet est entre z_l1 et z_l2", "\n"
									print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
							
						else:
						
							mu_lim = m_lim - Mlist_tild[j]
							if Realdata == 'Realdata_OK':
								# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
								z_l = Functions(0.7, 0.3).Nonlinear_Interpolation(zzz, ZETA_long, mu_lim) 
								V_0 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								vol.append(V_0)
								if z_l > zform:
									print "z_l > zform, z_l = ", z_l, "\n"
								else:
									pass
								V_l =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(min(zform, z_l), LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								vol_l.append(V_l)
								
								"""
								ZETA_jj = ZETA_long*1.0
								zz_jj = zzz*1.0
								z_l_jj = []
								for jj in range(10):
									z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim)
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
									condition2 = (abs(float(Distanguish_of_solutions[i]) - float(Distanguish_of_solutions[i-1]))>= 0.1)
									if condition1 == False and condition2 == True:
										z_l.append(z_l_l[i])
									else:
										pass
								if len(z_l)==3:
									
									z_1, z_2, z_3 = z_l[0], z_l[1], z_l[2]
									
									V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_l1 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_1, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_l2 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_2, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_l3 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_3, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

									#print "z_1, z_2, z_3 = ", z_1, z_2, z_3, "\n"
									
									if redshift[j] <= z_1:
									
										V_0 = V_z
										vol.append(V_0)
										V_l = V_l1 + V_l3 - V_l2
										vol_l.append(V_l)

									elif redshift[j] >= z_2:
										V_0 = V_l1 + V_z - V_l2
										vol.append(V_0)
										V_l = V_l1 + V_l3 - V_l2
										vol_l.append(V_l)
									else:
										print "le redshift de l'objet est entre z_l1 et z_l2 de K-correction", "\n"
										print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
										
								elif len(z_l)==1:
									z_l = z_l[0]
									V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									V_0 = V_z
									vol.append(V_0)
									V_l =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
									vol_l.append(V_l)
								else:
									print "len(z_l) n'est pas égale à 1 ni à 3, len(z_l) = ", len(z_l), ", z_l = ", z_l, "\n"
									print "z_l_l = ", z_l_l, "\n"
									print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
								"""
							else:
								z_l = Functions(0.7, 0.3).Nonlinear_Interpolation(zzz, ZETA_long, mu_lim) 
								V_0 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								vol.append(V_0)
								if z_l > zform:
									print "z_l > zform, z_l = ", z_l, "\n"
								else:
									pass
								V_l =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(min(zform, z_l), LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								vol_l.append(V_l)

				else:
					for j in range(len(redshift)):
					
						mu_lim = m_lim - Mlist_tild[j]
						if Realdata == 'Realdata_OK':
							
							# NOTE: on fait comme si la modification de ZETA due à la contribution de la k-correction n'éxiste pas, alors c'est un lissage pour les parties qui sont modifié sur le diagram M-V. (on fais ca car la méthode de récupération de des 4 valeurs ou de deux valeurs de z_l dans le cas de M_star<M<M_form et avec k-correction prends beaucoup du temps [en faisant ZETA de taille 10**6])
							z_l = Functions(0.7, 0.3).Nonlinear_Interpolation(zzz, ZETA_long, mu_lim) 
							V_0 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
							vol.append(V_0)
							if z_l > zform:
								print "z_l > zform, z_l = ", z_l, "\n"
							else:
								pass
							V_l =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(min(zform, z_l), LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
							vol_l.append(V_l)
							
							"""	
							ZETA_jj = ZETA*1.0
							zz_jj = zz*1.0
							z_l_jj = []
							for jj in range(10):
								z_ljj = Functions(0.7, 0.3).Nonlinear_Interpolation(zz_jj, ZETA_jj, mu_lim)
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
								condition2 = (abs(float(Distanguish_of_solutions[i]) - float(Distanguish_of_solutions[i-1]))>= 0.1)
								if condition1 == False and condition2 == True:
									z_l.append(z_l_l[i])
								else:
									pass

							if len(z_l)==3:
								
								z_1, z_2, z_3 = z_l[0], z_l[1], z_l[2]
								
								V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_l1 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_1, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_l2 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_2, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_l3 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_3, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]

								#print "z_1, z_2, z_3 = ", z_1, z_2, z_3, "\n"
								
								if redshift[j] <= z_1:
								
									V_0 = V_z
									vol.append(V_0)
									V_l = V_l1 + V_l3 - V_l2
									vol_l.append(V_l)

								elif redshift[j] >= z_2:
									V_0 = V_l1 + V_z - V_l2
									vol.append(V_0)
									V_l = V_l1 + V_l3 - V_l2
									vol_l.append(V_l)
								else:
									print "le redshift de l'objet est entre z_l1 et z_l2 de K-correction", "\n"
									print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
									
							elif len(z_l)==1:
								z_l = z_l[0]
								V_z = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								V_0 = V_z
								vol.append(V_0)
								V_l =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(z_l, LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
								vol_l.append(V_l)
							else:
								print "len(z_l) n'est pas égale à 1 ni à 3, len(z_l) = ", len(z_l), ", z_l = ", z_l, "\n"
								print "z_l_l = ", z_l_l, "\n"
								print "LambdaOfCorrNull, OmegaOfCorrNull = ", LambdaOfCorrNull[i], OmegaOfCorrNull[i], "\n"
							"""

						else:			
							z_l = Functions(0.7, 0.3).Nonlinear_Interpolation(zz, ZETA, mu_lim)  
							V_0 = Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(redshift[j], LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
							vol.append(V_0)
							if z_l > zform:
								print "z_l > zform, z_l = ", z_l, "\n"
							else:
								pass
							V_l =  Model(LambdaOfCorrNull[i], OmegaOfCorrNull[i]).curvature(min(zform, z_l), LambdaOfCorrNull[i], OmegaOfCorrNull[i], 0.0)[1]
							vol_l.append(V_l)
					
			else:
				print "You must to determine the condition of the Subsample as 'Subsample_OK' or 'Subsample_NO'. ", "\n"
						
		
			vol_l = array(vol_l)
			vol = array(vol)
			vv = vol/vol_l
			#Cov = sum(factorWeighting*((Mlist_tild - sum(factorWeighting*Mlist_tild))*(vv - sum(factorWeighting*vv))))
			Cov = 0
			#COV.append(Cov)
			vv.sort()
			VV = vv
			D = []
			F_V = []
			D_Proba = []
			#H = HH[bb]
			#for k in range(0, len(VV)+1):
			
			H = []   #valeurs de probability density function 
			FH = []  #valeurs de fonction de répartition
			#on va comparer avec une distribution uniform
			for bb in range(1,len(VV)+1):                                                                                                               
			      H.append(bb/float(len(VV)))

			H=array(H)
			for k in range(len(VV)):                                                                                                                   
				x = H[k]
				Heavisidelist = Functions(0.7,0.3).Heavisidefuction(x - H)
				F_x = sum(Heavisidelist)/float(len(VV)) # #valeurs de fonction de répartition théorique qu'on a pris uniform 
		
				Heavisidelist = Functions(0.7,0.3).Heavisidefuction(x - VV) 
				F_H_x = sum(Heavisidelist)/len(VV) # #valeurs de fonction de répartition en question que je compatre avec celle de la théorie 

				D.append(abs(F_H_x - F_x))
				
			
			#DD.append((LambdaOfCorrNull[i], OmegaOfCorrNull[i], max(D), sum(factorWeighting*VV), mean(VV), (max(VV) - min(VV))/2., Cov))
			#List_of_VV.append((LambdaOfCorrNull[i], OmegaOfCorrNull[i], max(D), sum(factorWeighting*VV), mean(VV), (max(VV) - min(VV))/2., Cov, VV, F_V))
			DD.append((LambdaOfCorrNull[i], OmegaOfCorrNull[i], max(D), mean(VV), (max(VV) - min(VV))/2., Cov))
			List_of_VV.append((LambdaOfCorrNull[i], OmegaOfCorrNull[i], max(D), mean(VV), (max(VV) - min(VV))/2., Cov, VV, F_V))

		#plt.plot(OmegaOfCorrNull, array(DD)+(LambdaOfCorrNull[DD.index(min(DD))] - min(DD)), color='m')
		#plt.plot(array(DD)+(OmegaOfCorrNull[DD.index(min(DD))] - min(DD)), LambdaOfCorrNull, color='y')
		#print "LambdaOfCorrNull[DD.index(min(DD))] = ", LambdaOfCorrNull[DD.index(min(DD))], "\n"
		#print "OmegaOfCorrNull[DD.index(min(DD))] = ", OmegaOfCorrNull[DD.index(min(DD))], "\n"
		
		DDD = list(array(DD)[:,2])
		plt.plot(OmegaOfCorrNull, array(DDD)+(LambdaOfCorrNull[DDD.index(min(DDD))] - min(DDD)), color='m')
		plt.plot(array(DDD)+(OmegaOfCorrNull[DDD.index(min(DDD))] - min(DDD)), LambdaOfCorrNull, color='y')
		print "LambdaOfCorrNull[DDD.index(min(DDD))] = ", LambdaOfCorrNull[DDD.index(min(DDD))], "\n"
		print "OmegaOfCorrNull[DDD.index(min(DDD))] = ", OmegaOfCorrNull[DDD.index(min(DDD))], "\n"
		print "DD = ", array(DD), "\n"
		
		return min(DDD)







def y4(xx):
	#v = (4.*xx**2 + 0.5*xx)*np.exp(-xx/2) + 0.3 -0.038961067977
	v = (4.*xx**2 + 1.51*xx)*np.exp(-(xx)/2) 
	return v


x = linspace(0.2,0.4,100)
x_inv = linspace(0.2,0.1,100)

X = array(list(x) + list(x_inv))
X = linspace(0.051,0.6, 100)
Y = y4(X)

print y4(0.3)

plt.figure()
N = 0
XY_Point = []
while N < 400: 
#for i in range(500):
        
	x_a = np.random.uniform(low=0.1, high=0.5, size=(1))[0]
	x_b = np.random.uniform(low=0.1, high=0.5, size=(1))[0]
	X = linspace(min(x_a, x_b), max(x_a, x_b), 100)
	
	eps_x = np.random.uniform(low=0.0, high=0.05, size=(1))[0]
	eps_y = np.random.uniform(low=0.0, high=0.05, size=(1))[0]

	if 0.2999 < np.mean(X) < 0.3001:
		if 0.695 < y4(np.mean(X)) < 0.7051:
			xp = X*1.0
			yp = Y*1.0

		
			#plt.plot(xp+eps_x-0.01, yp+eps_y-0.05)
			plt.plot(xp+eps_x, yp+eps_y)
			N += 1
			#XY_Point.append((xp+eps_x, yp+eps_y))
			#XY_Point.append(y4(np.mean(X))+eps_y)
			XY_Point.append(np.random.normal(0.3, 0.3, 1)[0])

	else:
		pass



X = linspace(0.1, 0.5, 20)
Y = linspace(0.01, 1.5, 20)
XY = []
for i in range(len(X)):
	for j in range(len(Y)):
		XY.append((X[i], Y[j]))
		
XY = array(XY)	
"""
densite_de_point = []
densite_de_point.append(0)
for i in range(len(XY)-1):
	larg = XY[i+1][0] - XY[i][0]
	longu = XY[i+1][1] - XY[i][1]
	nn = 0
	for k in range(len(XY_Point)):
		for j in range(100):
			yy = XY_Point[k][1][j]
			xx = XY_Point[k][0][j]
			if XY[i][1] < yy < XY[i+1][1]:
				nn += 1
			else:
				pass
		
	densite_de_point.append(nn)
"""
densite_de_point = XY_Point
densite_de_point = array(densite_de_point)

print "densite_de_point = ", densite_de_point

X_c = XY[:,0].reshape(20,20)
Y_c = XY[:,1].reshape(20,20)
densite_de_point = densite_de_point.reshape(20,20)

CS = pl.contour(X_c, Y_c, densite_de_point, 3, linewidths=np.arange(.5, 4, .5), colors=('r', 'green', 'blue', (1,1,0), '#afeeee', '0.5'))
plt.clabel(CS, inline=1, fontsize=10)
		