#!/usr/bin/python2.7.3
# -*- coding: latin-1 -*-

#Nom de fichier: HubbleDiagram.py
#==============================================================================
#title           :HubbleDiagram.py
#description     :This file gives us the apparent magnitude that allows us to build the Hubble diagram, it gives us also the weighting factor.
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
from scipy.integrate import quad
from scipy import stats
from scipy.stats import norm

import Methods
from Methods import Functions
from Cosmological_Model import Model

H_0 = 70.0 #Km.s-1.Mpc-1
LightVelocity = 2.99792458*10**5.0 # Km.s-1
d_H = (LightVelocity*10**6.0)/H_0  # Parsec

class ApparentMagnitude:
	"""Computing of theoretical apparent magnitude"""
	def __init__(self, z, lambda_0, omega_0, M, Realdata, degree_of_Polyn, Ncoefficients_of_KcorrPolyn):
		
		self.Realdata = Realdata
		self.degree_of_Polyn = degree_of_Polyn
		self.Ncoefficients_of_KcorrPolyn = Ncoefficients_of_KcorrPolyn
		
		self.param = Functions(lambda_0, omega_0)
		self.z = z
		self.M = M
		self.cosmoModel = Model(lambda_0, omega_0)
		self.d_L = self.cosmoModel.curvature(self.z, lambda_0, omega_0, M)[2] # calcule de distance de luminosité à partir de class "Model"
		self.comobileDistance = self.cosmoModel.curvature(self.z, lambda_0, omega_0, M)[0]
		#--------------------------------K-correction-------------------------------------
		if self.Realdata == 'Realdata_OK':
			#kcorr_em = self.param.PolyNdeg(self.degree_of_Polyn, self.Ncoefficients_of_KcorrPolyn, self.z)
			
			alfa_nu = -0.5
			kcorr_con = self.param.function_kcorr(alfa_nu, self.z)
			
			#self.kcorr = kcorr_con + kcorr_em
			self.kcorr = Ncoefficients_of_KcorrPolyn
			#----dérivé de kcorr------
			#d_kcorr_em = self.param.d_PolyNdeg(self.degree_of_Polyn, self.Ncoefficients_of_KcorrPolyn, self.z)
			d_kcorr_con = self.param.d_function_kcorr(alfa_nu, self.z)
			
			#self.d_kcorr = d_kcorr_con + d_kcorr_em
			self.d_kcorr = degree_of_Polyn
			#-------------------------
		elif self.Realdata == 'Realdata_NO': 
			pass
		      
		      
		else:
			print "You must to determine the condition of the Realdata as 'Realdata_OK' or 'Realdata_NO'. ", "\n"
		#----------------------------------------------------------------------------------

	def m_theor(self, z, lambda_0, omega_0, M, alfa, bbeta, strech, col):

		Condit = any(self.d_L < 0.0)
		if Condit == True:
			print "d_L is negatif '~' ", "\n"
		else:
			pass
		
		#print "d_L = ", d_L, "\n"
		if self.Realdata == 'Realdata_OK':
			
			m_theorie = self.M - alfa*(strech - 1) + bbeta*col + 5.*np.log10(self.d_L/10.0) + self.kcorr  # - 5.*np.log10((LightVelocity*10**6)/(10*H_0)) #d_L/10. --- (10. est en pc)
			
		elif self.Realdata == 'Realdata_NO': 
		
			m_theorie = self.M - alfa*(strech - 1) + bbeta*col + 5.*np.log10(self.d_L/10.0) # - 5.*np.log10((LightVelocity*10**6)/(10*H_0)) #d_L/10. --- (10. est en pc)

		else:
			print "You must to determine the condition of the Realdata as 'Realdata_OK' or 'Realdata_NO'. ", "\n"
				
		#print "m_theorie = ", m_theorie, "\t"
		return m_theorie

	
	#==============================Jacobian==============================================


	def Jacobian(self, z):
	  
		#NOTE: For the reason that the case of the QSOs doesn't need to make a modification on the weighting factors after doing the K-correction, we replace 'Realdata_OK' by 'Realdata_NO' 
		#self.Realdata = 'Realdata_NO'

		if self.Realdata == 'Realdata_OK':
			
			tau = self.comobileDistance
			jac_1 = 0.0
			
			if self.param.kappa_0 == 0.0:
				jac_1 = (1./(tau**2.))*( (1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (1./tau) + self.d_kcorr*sqrt(self.param.P(1./(1. + self.z)))*(1. + self.z)**2. )
				

			elif self.param.kappa_0 > 0.0:
				#In the case where kappa0 > 0, the sin function gives us a negatif value with angles > pi (redshift of atipodes) (angle changed by (2*pi - angle) or sin(x) by abs(sin(x)). If we do not do it, it dangerous to obtain a negative values ​​for d_L  (changelly of tan(x) by : abs(sin(x))/cos(x))
				
				epsilon = 2*tau*sqrt(self.param.kappa_0)
				Condition = any(epsilon < 10**-7.7)
				
				if Condition == True:
					if type(epsilon) is not np.ndarray:
						if epsilon < 10**-7.7:
							jac_1 = ((1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (sqrt(self.param.kappa_0)/tan(epsilon/2.0) ))/((epsilon**2.)/2.)
		
						else:
							jac_1 = ((1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (sqrt(self.param.kappa_0)/(abs(sin(epsilon/2.0))/cos(epsilon/2.0) )))/(1. - cos(epsilon))
					else:
						jac_1 = []
						for i in range(len(epsilon)):
							if epsilon[i] < 10**-7.7:
								jac1 = ((1. + self.z[i])*sqrt(self.param.P(1./(1. + self.z[i]))) + (sqrt(self.param.kappa_0)/tan(epsilon[i]/2.0) ))/((epsilon[i]**2.)/2.)
								jac_1.append(jac1)
							
							else:
								jac1 = ((1. + self.z[i])*sqrt(self.param.P(1./(1. + self.z[i]))) + (sqrt(self.param.kappa_0)/(abs(sin(epsilon[i]/2.0))/cos(epsilon[i]/2.0) ) ))/(1. - cos(epsilon[i]))
								jac_1.append(jac1)

						jac_1 = array(jac_1)
				else:
					jac_1 = ((1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (sqrt(self.param.kappa_0)/(abs(sin(epsilon/2.0))/cos(epsilon/2.0) ) ) + ((((1. + self.z)**2.)*sqrt(self.param.P(1./(1. + self.z))))*self.d_kcorr))/(1. - cos(epsilon))	

				
			
			elif self.param.kappa_0 < 0.0 :
				epsilon = 2.*tau*sqrt(abs(self.param.kappa_0))
				Condition = any(epsilon < 10**-7.7)
				if Condition == True:
					if type(epsilon) is not np.ndarray:
						if epsilon < 10**-7.7:
							jac_1 = ((1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (sqrt(abs(self.param.kappa_0))/tanh(epsilon/2.0)))/(((epsilon**2.)/2.))
						else:
							jac_1 = ((1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (sqrt(abs(self.param.kappa_0))/tanh(epsilon/2.0)))/(cosh(epsilon) - 1.)
					else:
						jac_1 = []
						for i in range(len(epsilon)):
							if epsilon[i] < 10**-7.7:
								jac1 = ((1. + self.z[i])*sqrt(self.param.P(1./(1. + self.z[i]))) + (sqrt(abs(self.param.kappa_0))/tanh(epsilon[i]/2.0)))/(((epsilon[i]**2.)/2.))
								jac_1.append(jac1)
							else:
								jac1 = ((1. + self.z[i])*sqrt(self.param.P(1./(1. + self.z[i]))) + (sqrt(abs(self.param.kappa_0))/tanh(epsilon[i]/2.0)))/(cosh(epsilon[i]) - 1.)
								jac_1.append(jac1)

						jac_1 = array(jac_1)
				else:
					jac_1 = ((1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (sqrt(abs(self.param.kappa_0))/tanh(epsilon/2.0)) + ((((1. + self.z)**2.)*sqrt(self.param.P(1./(1. + self.z))))*self.d_kcorr))/(cosh(epsilon) - 1.)
				
			else:
				pass
		  
		  
		elif self.Realdata == 'Realdata_NO':
		    
			#tau = self.cosmoModel.curvature(self.z, self.param.lambda_0, self.param.omega_0, M)[0] # M n'influe pas sur tau, comme si elle n'existe pas
			#dtau = (1 + self.z)*self.param.FunctiontoIntegrate(1./(1+z))
			tau = self.comobileDistance
			jac_1 = 0.0
			
			if self.param.kappa_0 == 0.0:
				jac_1 = (1./(tau**3.))*( tau*(1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + 1. )
				#jac_1 = 1.0/tau**3.0
				#print "kappa0 = ", self.param.kappa_0
				

			elif self.param.kappa_0 > 0.0:
				#In the case where kappa0 > 0, the sin function gives us a negatif value with angles > pi (redshift of atipodes) (angle changed by (2*pi - angle) or sin(x) by abs(sin(x)). If we do not do it, it dangerous to obtain a negative values ​​for d_L  (changelly of tan(x) by : abs(sin(x))/cos(x))
				
				epsilon = 2*tau*sqrt(self.param.kappa_0)
				Condition = any(epsilon < 10**-7.7)
				
				if Condition == True:
					if type(epsilon) is not np.ndarray:
						if epsilon < 10**-7.7:
							jac_1 = ((1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (sqrt(self.param.kappa_0)/tan(epsilon/2.0) ))/((epsilon**2.)/2.)
		
						else:
							jac_1 = ((1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (sqrt(self.param.kappa_0)/(abs(sin(epsilon/2.0))/cos(epsilon/2.0) )))/(1. - cos(epsilon))
					else:
						jac_1 = []
						for i in range(len(epsilon)):
							if epsilon[i] < 10**-7.7:
								jac1 = ((1. + self.z[i])*sqrt(self.param.P(1./(1. + self.z[i]))) + (sqrt(self.param.kappa_0)/tan(epsilon[i]/2.0) ))/((epsilon[i]**2.)/2.)
								jac_1.append(jac1)
							
							else:
								jac1 = ((1. + self.z[i])*sqrt(self.param.P(1./(1. + self.z[i]))) + (sqrt(self.param.kappa_0)/(abs(sin(epsilon[i]/2.0))/cos(epsilon[i]/2.0) ) ))/(1. - cos(epsilon[i]))
								jac_1.append(jac1)

						jac_1 = array(jac_1)
				else:
					jac_1 = ((1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (sqrt(self.param.kappa_0)/(abs(sin(epsilon/2.0))/cos(epsilon/2.0) ) ))/(1. - cos(epsilon))	
				
				#(5.*self.param.kappa_0)/(2.*pi*np.log(10)) cette constante sera simplifié qd on divise par somme(w_k)
				#print "kappa0 = ", self.param.kappa_0, "; jac_1 = ", jac_1, "\n"
				
			
			elif self.param.kappa_0 < 0.0 :
				epsilon = 2.*tau*sqrt(abs(self.param.kappa_0))
				Condition = any(epsilon < 10**-7.7)
				if Condition == True:
					if type(epsilon) is not np.ndarray:
						if epsilon < 10**-7.7:
							jac_1 = ((1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (sqrt(abs(self.param.kappa_0))/tanh(epsilon/2.0)))/(((epsilon**2.)/2.))
						else:
							jac_1 = ((1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (sqrt(abs(self.param.kappa_0))/tanh(epsilon/2.0)))/(cosh(epsilon) - 1.)
					else:
						jac_1 = []
						for i in range(len(epsilon)):
							if epsilon[i] < 10**-7.7:
								jac1 = ((1. + self.z[i])*sqrt(self.param.P(1./(1. + self.z[i]))) + (sqrt(abs(self.param.kappa_0))/tanh(epsilon[i]/2.0)))/(((epsilon[i]**2.)/2.))
								jac_1.append(jac1)
							else:
								jac1 = ((1. + self.z[i])*sqrt(self.param.P(1./(1. + self.z[i]))) + (sqrt(abs(self.param.kappa_0))/tanh(epsilon[i]/2.0)))/(cosh(epsilon[i]) - 1.)
								jac_1.append(jac1)

						jac_1 = array(jac_1)
				else:
					jac_1 = ((1. + self.z)*sqrt(self.param.P(1./(1. + self.z))) + (sqrt(abs(self.param.kappa_0))/tanh(epsilon/2.0)))/(cosh(epsilon) - 1.)

				#print "kappa0 = ", self.param.kappa_0, "; jac_1 = ", jac_1, "\n"
				
			else:
				pass
		else:
			print "You must to determine the condition of the Realdata as 'Realdata_OK' or 'Realdata_NO'. ", "\n"
			
		return abs(jac_1)

	#==============================Weighting factor==============================================


	def Weightingfactor(self, z, beta):
		
		#NOTE: For the reason that the case of the QSOs doesn't need to make a modification on the weighting factors after doing the K-correction, we replace 'Realdata_OK' by 'Realdata_NO' 
		#self.Realdata = 'Realdata_NO'
		
		#tau = self.cosmoModel.curvature(self.z, self.param.lambda_0, self.param.omega_0, M)[0] # M n'influe pas sur tau, comme si elle n'existe pas
		#zeta = self.m_theor(z, self.param.lambda_0, self.param.omega_0, 0)
		tau = self.comobileDistance
		#d_L = self.cosmoModel.curvature(self.z, self.param.lambda_0, self.param.omega_0, 0.0)[2] # calcule de distance de luminosité à partir de class "Model"
		#weighing = abs(self.Jacobian(z, M))*( 10**(beta*zeta/5.) )#/( 10**(beta*zeta/5.) )
		
		if self.Realdata == 'Realdata_OK':
		  
			if self.param.kappa_0 == 0.0:
				weighing = abs(self.Jacobian(z))*( ((1. + z)*tau*10.**(self.kcorr))**beta )
					
			elif self.param.kappa_0 > 0.0:
				epsilon = tau*sqrt(self.param.kappa_0)	
				weighing = abs(self.Jacobian(z))*( ((1. + z)*abs(sin(epsilon))*10.**(self.kcorr)/sqrt(self.param.kappa_0))**beta )
				
			elif self.param.kappa_0 < 0.0 :
				epsilon = tau*sqrt(abs(self.param.kappa_0))
				weighing = abs(self.Jacobian(z))*( ((1. + z)*sinh(epsilon)*10.**(self.kcorr)/sqrt(abs(self.param.kappa_0)))**beta )
			else:
				pass
		
		
		elif self.Realdata == 'Realdata_NO':
		  
			if self.param.kappa_0 == 0.0:
				weighing = abs(self.Jacobian(z)) *( ((1. + z)*tau)**beta ) #/ ( ((1. + z)*tau)**beta )
				
			elif self.param.kappa_0 > 0.0:
				#weighing = ( self.param.kappa_0**(1. - beta/2.) )*abs(self.Jacobian(z, M))*( ((1. + z)*sin(tau*sqrt(self.param.kappa_0)))**beta ) #/ ( ((1. + z)*sin(tau*sqrt(self.param.kappa_0))/sqrt(self.param.kappa_0))**beta )
				#weighing = abs(self.Jacobian(z, M))*( ((1. + z)*sin(tau*sqrt(self.param.kappa_0)))**beta ) #/ ( ((1. + z)*sin(tau*sqrt(self.param.kappa_0))/sqrt(self.param.kappa_0))**beta )
				epsilon = tau*sqrt(self.param.kappa_0)	
				weighing = abs(self.Jacobian(z)) *( ((1. + z)*abs(sin(epsilon))/sqrt(self.param.kappa_0))**beta )	#/( ((1. + z)*abs(sin(epsilon))/sqrt(self.param.kappa_0))**beta )
				
			elif self.param.kappa_0 < 0.0 :
				#weighing = ( abs(self.param.kappa_0)**(1. - beta/2.) )*abs(self.Jacobian(z, M))*( ((1. + z)*sinh(tau*sqrt(abs(self.param.kappa_0))))**beta ) #/ ( ((1. + z)*sinh(tau*sqrt(abs(self.param.kappa_0)))/sqrt(abs(self.param.kappa_0)))**beta )
				#weighing = abs(self.Jacobian(z, M))*( ((1. + z)*sinh(tau*sqrt(abs(self.param.kappa_0))))**beta ) #/ ( ((1. + z)*sinh(tau*sqrt(abs(self.param.kappa_0)))/sqrt(abs(self.param.kappa_0)))**beta )
				epsilon = tau*sqrt(abs(self.param.kappa_0))
				weighing = abs(self.Jacobian(z)) *( ((1. + z)*sinh(epsilon)/sqrt(abs(self.param.kappa_0)))**beta )
			else:
				pass
			      
		else:
			print "You must to determine the condition of the Realdata as 'Realdata_OK' or 'Realdata_NO'. ", "\n"
		 
		return weighing
		
		
		
	def VOLUME(self, M, zz, ZETA, m_l):		# m_l : magnitude apparente limite.
		
		z_l = self.param.Nonlinear_Interpolation(zz, ZETA, m_l - M)	#z_l : redshift limite (z_max). ATTENTION: zz et ZETA doivent etre calculé par le bon model cosmologique, (modèle de simulation.)
		
		if self.param.kappa_0 == 0.0:
			mu = M + 5*log10(((LightVelocity*10**6.0)/(10.*H_0)))
			if self.Realdata == 'Realdata_OK':
				V = (4.*pi/3.)*((1/(1. + z_l))**3.)*10**(3.*(m_l - mu - self.kcorr)/5.)
			
			elif self.Realdata == 'Realdata_NO':
				V = (4.*pi/3.)*((1/(1. + z_l))**3.)*10**(3.*(m_l - mu)/5.)
			else:
				print "You must to determine the condition of the Realdata as 'Realdata_OK' or 'Realdata_NO'. ", "\n"
			
		elif self.param.kappa_0 > 0.0:
			
			#alfa_l = ((10*H_0)/(LightVelocity*10**6.0))*(sqrt(self.param.kappa_0)/(1 + z_l))*(10**((m_l - M)/5.))
			mu = M + 5*log10(((LightVelocity*10**6.0)/(10.*H_0)))
						
			if self.Realdata == 'Realdata_OK':
				
				alfa_l = (sqrt(self.param.kappa_0)/(1 + z_l))*(10**((m_l - mu - self.kcorr)/5.))	# ((10*H_0)/(LightVelocity*10**6.0)) must be calculated in µ = M + 5*log10(((LightVelocity*10**6.0)/(10*H_0))) ce qui rend m - µ (au lieu de m - M) petit et du coup alfa_l < 1. (on fait ca parce que les formule de Volume sont déjà sans dimension).

			elif self.Realdata == 'Realdata_NO':
				
				alfa_l = (sqrt(self.param.kappa_0)/(1 + z_l))*(10**((m_l - mu)/5.))	# ((10*H_0)/(LightVelocity*10**6.0)) must be calculated in µ = M + 5*log10(((LightVelocity*10**6.0)/(10*H_0))) ce qui rend m - µ (au lieu de m - M) petit et du coup alfa_l < 1. (on fait ca parce que les formule de Volume sont déjà sans dimension).
			else:
				print "You must to determine the condition of the Realdata as 'Realdata_OK' or 'Realdata_NO'. ", "\n"
			
			V = (pi/self.param.kappa_0**(3./2))*(2*arcsin(alfa_l) - 2*alfa_l*sqrt(1 - alfa_l**2.)) 
			
		elif self.param.kappa_0 < 0.0 :
		  
			mu = M + 5*log10(((LightVelocity*10**6.0)/(10.*H_0)))
			alfa_l = (sqrt(abs(self.param.kappa_0))/(1 + z_l))*(10**((m_l - mu)/5.))
			
			if self.Realdata == 'Realdata_OK':
				
				alfa_l = (sqrt(abs(self.param.kappa_0))/(1 + z_l))*(10**((m_l - mu - self.kcorr)/5.))	# ((10*H_0)/(LightVelocity*10**6.0)) must be calculated in µ = M + 5*log10(((LightVelocity*10**6.0)/(10*H_0))) ce qui rend m - µ (au lieu de m - M) petit et du coup alfa_l < 1. (on fait ca parce que les formule de Volume sont déjà sans dimension).

			elif self.Realdata == 'Realdata_NO':
				
				alfa_l = (sqrt(abs(self.param.kappa_0))/(1 + z_l))*(10**((m_l - mu)/5.))	# ((10*H_0)/(LightVelocity*10**6.0)) must be calculated in µ = M + 5*log10(((LightVelocity*10**6.0)/(10*H_0))) ce qui rend m - µ (au lieu de m - M) petit et du coup alfa_l < 1. (on fait ca parce que les formule de Volume sont déjà sans dimension).
			else:
				print "You must to determine the condition of the Realdata as 'Realdata_OK' or 'Realdata_NO'. ", "\n"
		
			V = (pi/abs(self.param.kappa_0)**(3./2))*(2*alfa_l*sqrt(1 + alfa_l**2.) - 2*arcsinh(alfa_l))
			
		else:    
			pass
		  
		return V	