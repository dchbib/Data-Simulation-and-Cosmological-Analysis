#!/usr/bin/python2.7.3
# -*- coding: latin-1 -*-

# Nom de fichier: Methods.py
#==============================================================================
#title           :Methods.py
#description     :This file contain all methods of calculation and interpolation that we need in our simulation.
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
from math import gamma
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


H_0 = 70.0 #Km.s-1.Mpc-1
LightVelocity = 2.99792458*10**5.0 # Km.s-1
alfa_0 = 10**(-5.)

class Functions:
	"""All the functions and methods that I need in my simulation"""
	def __init__(self,lambda_0, omega_0):
		self.lambda_0 = lambda_0
		self.omega_0 = omega_0
		self.kappa_0 = self.lambda_0 + self.omega_0 - 1.

	def P(self, a):
	
		self.a = a
		value = self.lambda_0*a**4. - self.kappa_0*a**2. + self.omega_0*a 	# + alfa_0
		#print "P(de class) = ", value, "\n"
		return value

	def FunctiontoIntegrate(self, a):
	
		try:
			value = 1./sqrt(self.P(a))
		except RuntimeWarning:
			print "lambda0, omega0 = ", self.lambda_0, self.omega_0
			
		return value
		
		
	def FunctiontoIntegrate2(self, a):
		
		value = a/sqrt(self.P(a))
		return value	


	def Functiontodmth_lambda0(self, a):

		value = (a**2. - a**4.)/(self.P(a))**(3/2.)
		return value

	def Functiontodmth_omega0(self, a):

		value = (a**2 - a)/(self.P(a))**(3/2.)
		return value

#================derivative of the inverse function is used in the method of Newton==============

	def LowerBound(self, z):

		return 1./(1. + z)


	def dIntegral(self, z):
	 		
		a = self.LowerBound(z)
		#print "dIntegral = ", (a**2)*self.FunctiontoIntegrate(a), "\n"
		return (a**2)*self.FunctiontoIntegrate(a)

#===============================simpson method================================================

	def simpsons_rule(self, f, a, b):

		self.f = f
		c = (a+b) / 2.0
		h = abs(b-a) / 6.0
		valuee = h*(self.f(a) + 4.0*self.f(c) + self.f(b))
		#print "simpsons_rule = ", valuee, "\n"
		return valuee

	
	def Trapeze(self, f, a, b):

		self.f = f
		#n = 100000
		#k = linspace(1,n-1,100000)
		if type(a) is np.ndarray:	
			def intgrl1(LowerLimit, UpperLimit):
				return quad(self.f, LowerLimit, UpperLimit)
				
			vintgrl1 = vectorize(intgrl1) # ce qui permet d'utiliser quad avec des bornes de type array
			
			valuee = vintgrl1(a, b)[0]
			#valuee = array([quad(self.f, r, b)[0] for r in a])
			#valuee = array([((b-r)/n)*((self.f(r)/2.) + sum(self.f(r + k*((b-r)/n))) + (self.f(b)/2.)) for r in a])
		else:
			"""
			h = (b-a) / n
			k = a + k*h
			c = sum(self.f(k))
			valuee = h*(self.f(a)/2. + c + self.f(b)/2.)
			"""
			valuee = quad(self.f, a, b)[0]

		return valuee
		
	def Trapeze_3D(self, f, a, b):

		self.f = f
		#n = 100000
		#k = linspace(1,n-1,100000)
		if type(a) is np.ndarray:	
			def intgrl1(LowerLimit, UpperLimit):
				return quad(self.f, LowerLimit, UpperLimit)
				
			vintgrl1 = vectorize(intgrl1) # ce qui permet d'utiliser quad avec des bornes de type array
			
			valuee = vintgrl1(a, b)[0]
			#valuee = array([quad(self.f, r, b)[0] for r in a])
			#valuee = array([((b-r)/n)*((self.f(r)/2.) + sum(self.f(r + k*((b-r)/n))) + (self.f(b)/2.)) for r in a])
		else:
			"""
			h = (b-a) / n
			k = a + k*h
			c = sum(self.f(k))
			valuee = h*(self.f(a)/2. + c + self.f(b)/2.)
			"""
			valuee = quad(self.f, a, b)[0]

		return valuee
	
#=============================computing chi2=================================

	def Chisquare(self, f_obs, f_exp, sigma):
 
		chisq = 0.
		for i in range(len(f_obs)):
			chisq += ((f_obs[i] - f_exp[i])/sigma)**2
		print "chisq = ", chisq
		return chisq

#========================================================================================

# Methode de Newton Raphson pour trouver racines d'une fonction de dimension 1
	def Newton_Raphson(self, f, df, x0, eps, Nmax, meanM, meanMsquare):
		j=0
		x=x0
		listofx = []
		listoffunctionx = []
		while (j<Nmax and abs(f(x, meanM, meanMsquare))>eps):
			x = x - (f(x, meanM, meanMsquare)/df(x, meanM, meanMsquare))
        		#print x, "\n"
			#ùplt.plot(x, f(x, meanM, meanMsquare), marker = '+', color = 'r') #*
			listofx.append(x)
			listoffunctionx.append(f(x, meanM, meanMsquare))
			j = j+1
		#plt.ylim(-50, 50)
		#plt.xlim(-22, -18)
		#ùplt.show() #*
		#c = listoffunctionx.index(min(listoffunctionx))
		#print "M0 = ", listofx[c], "\n"
		#return listofx[c];
		#print "f(x, meanM, meanMsquare) = ", f(x, meanM, meanMsquare), "\n"
		#print "-------------------------------------------", "\n"
		return x


#==============================bounced limit================================================

	def NoBigBang(self, omega0, condition):
		
		if condition == 'permission':
			a = (27./4.)*(self.omega_0**2.)*self.lambda_0
			#b = (lambda0 + omega0 + alfa_0 - 1.)**3.
			b = self.kappa_0**3.
			Condition = a>b
			c=0
			if Condition == True:
				c=1		# that means that the parameters are acceptable
			else:
				c=0		# that means that the parameters are not acceptable, we don't have a singularity of BigBang.
				
			return  c
		elif condition == 'limit':
		
			if omega0 < 0.5:
				lambda0_lim = ( 4.0*omega0*(cosh((1./3.)*arccosh((1.0 - omega0)/omega0)))**3.0 ) #- 0.01

			elif omega0 > 0.5:
				lambda0_lim = ( 4.0*omega0*(cos((1./3.)*arccos((1.0 - omega0)/omega0)))**3.0 ) # - 0.01

			else:	# omega0 == 0.5
				lambda0_lim = 2.0


			return lambda0_lim
		
		else:
			print "You must to determine the condition of the NoBigBang's function as 'permission' or 'limit'. ", "\n"
			
			
#==============================Heaviside fuction================================================

	def Heavisidefuction(self, t):
		Condition1 = type(t) is np.ndarray
		Condition2 = type(t) is list
		if Condition1 == True or Condition2 == True :
			Heavisidelist = []
			for i in range(len(t)):
				c = 0.0
				if t[i] >= 0.0:
					c = 1
				else:
					c = 0.0
				Heavisidelist.append(c)

			b = np.asarray(Heavisidelist)
			
		else:
			c = 0.0
			if t >= 0.0:
				c = 1
			else:
				c = 0.0
			b = c

		return b


#==============================quadraticsolver================================================

	def quadraticsolver(self, a,b,c):
		d = b**2. - 4.*a*c
		if d < 0: 
		
			#print "This equation has no real solution"
			return 0
		elif d > 0:
			s0 = (-b - sqrt(d)) / (2.*a)
			s1 = (-b + sqrt(d)) / (2.*a)
			
		elif d == 0:
			s0 = -b / (2.*a)
			s1 = s0 
		else:
			pass
		      
		      
		return s0, s1


#==============================Nonlinear interpolation================================================


	def Nonlinear_Interpolation(self, XInitial_list, YInitial_list, Searching_value):	# XInitial_list, YInitial_list, Searching_value(c'est le True_yvalue)
		
		True_xvalue = 0.0
		
		if Searching_value == 0.0: 
		      
			#~YInitial_list_negative = filter(lambda x:x<0,YInitial_list)
			#~YInitial_list_positive = filter(lambda x:x>0,YInitial_list)
			# **si la liste commenece de négative au positive*** // on suppose que (comme dans le reste de la fonction) que si la ligne est divisé entre deux partie positive et négative, alors cette ligne est bien organisé et les valeurs de signe différents ne sont pas mélangé, exemple: [-1,-2,-0.5,-6,-12,-0.1, 0.5,2,4,15,3,2,10]			
			#~indice = len(YInitial_list_negative)-1  #; or indice =  len(YInitial_list) -  len(YInitial_list_positive)
			# **si la liste commenece de positive au négative ***
			# OR indice = len(YInitial_list_positive) - 1 #; or indice =  len(YInitial_list) -  len(YInitial_list_negative)
			
			indice = list(abs(YInitial_list)).index(min(abs(YInitial_list)))
		
		else:
			indice = list(abs(YInitial_list - Searching_value)).index(min(abs(YInitial_list - Searching_value)))
		
		
		Y_k = YInitial_list[indice]
		X_k = XInitial_list[indice]
		
		if indice == 0 or indice == len(YInitial_list)-1 or indice == 1 or indice == len(YInitial_list)-2 or indice == 2 or indice == len(YInitial_list)-3:
			print "The indice is one of three indices condition", "\n"
			True_xvalue = X_k
      
		else:	 
			#----------changing of method------------------------
			Difference = len(XInitial_list) - indice

			Y_k_moins1 = YInitial_list[indice-1]
			Y_k_plus1 = YInitial_list[indice+1]

			#-----------------------------------
			X_k_moins1 = XInitial_list[indice-1]
			X_k_plus1 = XInitial_list[indice+1]

			try:
				# solve of the following system A.x = b
				A=matrix([[X_k_moins1**2.0, X_k_moins1, 1.],[X_k**2., X_k, 1.], [X_k_plus1**2., X_k_plus1, 1.]])
				b=matrix([[Y_k_moins1],[Y_k], [Y_k_plus1]])

				solution = array(linalg.solve(A,b)) # gives the three constants needed.
				
				# refind the true x-value correspond to the Searching_value (y-value)
				Solve_poly2deg = self.quadraticsolver(solution[0][0], solution[1][0], solution[2][0]-Searching_value)
				ii=1
				
				while type(Solve_poly2deg) is int and ii < min(Difference, indice)-1:
					#print "ii = ", ii, "\n"
					#print "Difference = ", Difference, "\n"
					#print "min(Difference, indice) = ", min(Difference, indice), "\n"
					#print "indice = ",indice, "\n"
					ii = ii+1

					Y_k_moins1 = YInitial_list[indice-ii]
					Y_k_plus1 = YInitial_list[indice+ii]
					#-----------------------------------
					
					X_k_moins1 = XInitial_list[indice-ii]
					X_k_plus1 = XInitial_list[indice+ii]

					try:
						# solve of the following system A.x = b
						A=matrix([[X_k_moins1**2.0, X_k_moins1, 1.],[X_k**2., X_k, 1.], [X_k_plus1**2., X_k_plus1, 1.]])
						b=matrix([[Y_k_moins1],[Y_k], [Y_k_plus1]])

						solution = array(linalg.solve(A,b)) # gives the three constants needed.
						
						# refind the true x-value correspond to the Searching_value (y-value)
						Solve_poly2deg = self.quadraticsolver(solution[0][0], solution[1][0], solution[2][0]-Searching_value)
					except LinAlgError:
						True_xvalue = X_k
						
				if type(Solve_poly2deg) is int:
                                        #print "BaaaBaaacoucou!!"
					True_xvalue = X_k
				else:
					#print "ii = ", ii, "\n"
					for i in range(2):
						if Solve_poly2deg[i] < min(XInitial_list) or Solve_poly2deg[i] > max(XInitial_list):
							continue
						elif Solve_poly2deg[i] > min(XInitial_list) and Solve_poly2deg[i] < max(XInitial_list):
                                                        #print "coucouBaaa!!"
							True_xvalue = Solve_poly2deg[list(abs(array(Solve_poly2deg) - X_k)).index(min(abs(array(Solve_poly2deg) - X_k)))]
							#True_xvalue = Solve_poly2deg[list(abs(Solve_poly2deg - X_k)).index(min(abs(Solve_poly2deg - X_k)))]
							#True_xvalue = Solve_poly2deg[np.array([abs(Solve_poly2deg[0] - X_k),abs(Solve_poly2deg[1] - X_k)]).argmin()]
							continue
						else:
							pass
					      
			except LinAlgError:
				True_xvalue = X_k
				
		return True_xvalue
		
#========================================================================================

#==============================Nonlinear interpolation================================================

	def Nonlinear_Interpolation_one_by_one(self, XInitial_list, YInitial_list, Searching_value):	# XInitial_list, YInitial_list, Searching_value(c'est le True_yvalue)
		
		True_xvalue = 0.0
		
		indice = 1
		Y_k = YInitial_list[indice]
		X_k = XInitial_list[indice]
		
		#-----------------------------------
		Y_k_moins1 = YInitial_list[indice-1]
		Y_k_plus1 = YInitial_list[indice+1]

		#-----------------------------------
		X_k_moins1 = XInitial_list[indice-1]
		X_k_plus1 = XInitial_list[indice+1]
		
		try:
			# solve of the following system A.x = b
			A=matrix([[X_k_moins1**2.0, X_k_moins1, 1.],[X_k**2., X_k, 1.], [X_k_plus1**2., X_k_plus1, 1.]])
			b=matrix([[Y_k_moins1],[Y_k], [Y_k_plus1]])

			solution = array(linalg.solve(A,b)) # gives the three constants needed.
					
			# refind the true x-value correspond to the Searching_value (y-value)
			Solve_poly2deg = self.quadraticsolver(solution[0][0], solution[1][0], solution[2][0]-Searching_value)
				
			if type(Solve_poly2deg) is int:
                                #print "BaaaBaaacoucou!!"
				True_xvalue = X_k
			else:
				#print "ii = ", ii, "\n"
				for i in range(2):
					if Solve_poly2deg[i] < min(XInitial_list) or Solve_poly2deg[i] > max(XInitial_list):
						continue
					elif Solve_poly2deg[i] > min(XInitial_list) and Solve_poly2deg[i] < max(XInitial_list):
                                                #print "coucouBaaa!!"
                                                True_xvalue = Solve_poly2deg[list(abs(array(Solve_poly2deg) - X_k)).index(min(abs(array(Solve_poly2deg) - X_k)))]
						#True_xvalue = Solve_poly2deg[np.array([abs(Solve_poly2deg[0] - X_k),abs(Solve_poly2deg[1] - X_k)]).argmin()]
						continue
					else:
						pass
		except LinAlgError:
			True_xvalue = X_k
				
		return True_xvalue
#========================================================================================


	def gaussiane(self, M, M0, sigma):
		v = (1./sqrt(2.*pi*sigma**2.))*(np.exp(-0.5*((M - M0)/sigma)**2.))
		return v
#========================================================================================

	def G_M_s_c(self, Mx, sy, cz, M0, s0, c0, sigmaM, sigmas, sigmac,alfa, bbeta):
		v = np.exp(-(Mx - alfa*(sy-1) + bbeta*cz -M0 -s0 -c0)**2./(2*(sigmaM**2. + sigmas**2. + sigmac**2.)))/sqrt(2*pi*(sigmaM**2. + sigmas**2. + sigmac**2.))
		return v
		
#========================================================================================

	#def G_M_s_c_2D(self, sy, cz, M0, s0, c0, sigmaM, sigmas, sigmac,alfa, bbeta):
	#	v = np.exp(-(- alfa*(sy-1) + bbeta*cz - s0 - c0)**2./(2*(sigmaM**2. + sigmas**2. + sigmac**2.)))/sqrt(2*pi*(sigmaM**2. + sigmas**2. + sigmac**2.))
	#	return v
		
	def G_M_s_c_2D(self, sy, cz, M0, s0, c0, sigmaM, sigmas, sigmac,alfa, bbeta):
                    v = np.exp(-0.5*(((sy - s0)/sigmas)**2. + ((cz - c0)/sigmac)**2. - 2*(alfa/bbeta)*((sy - s0)*(cz - c0)/(sigmas*sigmac))))/(2*pi*sigmas*sigmac*sqrt(1-(alfa/bbeta)**2.))
                    return v
	
#========================================================================================

	def G_M_s_c_sum_gauss(self, M, s, c, M0, s0, c0, sigmaM, sigmas, sigmac):
		v = self.gaussiane(M, M0, sigmaM) + self.gaussiane(s, s0, sigmas) + self.gaussiane(c, c0, sigmac)    
		return v
		
#========================================================================================

	def FHI_m(self, m, m_lim, beta):
		v = 1. +  self.Heavisidefuction(m_lim - m)*(10**(-beta*(m_lim - m)/5.) - 1)
		return v
#========================================================================================

	def Gumbel(self, M, M0, sigma):
		v = (1./sigma)*np.exp(-(M - M0)/sigma)*np.exp(-np.exp(-(M - M0)/sigma))
		return v
		
#========================================================================================

	def Exponential(self, M, a):
		v = np.exp(a*M)
		return v
#========================================================================================

	def Integral_Exponential(self, M, a):
		v = np.exp(a*M)/a
		return v
		
#========================================================================================

	def Gamma(self, M, shape, scale):
		v = M**(shape-1)*(np.exp(-M/scale)/(scale**(shape)*gamma(shape)))
		return v
		
#========================================================================================
	def ksi(self, M):                                
		v = 0.065*M +1.7034                                
		return v

	def Luminosity_function(self, M, M0, sigma, beta_over_five):		#Luminosity_function_w
		#const = self.ksi(M)
		const = -0.05
		v = ((((1./sigma)*(1 + const*((M-M0)/sigma))**((-1./const)-1))*np.exp(-(1 + const*((M-M0)/sigma))**(-1./const)))*self.Heavisidefuction((-24.5) - M) + (((1./sigma)*(1 + const*((M-M0)/sigma))**((-1./const)-1))*np.exp(-(1 + const*((M-M0)/sigma))**(-1./const))*10**(M/10.) + 0.078)*self.Heavisidefuction(M - (-24.5)) )
		return v
		
	def Luminosity_function_time_Tenpower(self, M, M0, sigma, beta_over_five):
		const = self.ksi(M)
		v = ((((1./sigma)*(1 + const*((M-M0)/sigma))**((-1./const)-1))*np.exp(-(1 + const*((M-M0)/sigma))**(-1./const)))*self.Heavisidefuction((-24.5) - M) + (((1./sigma)*(1 + const*((M-M0)/sigma))**((-1./const)-1))*np.exp(-(1 + const*((M-M0)/sigma))**(-1./const))*10**(M/10.) + 0.078)*self.Heavisidefuction(M - (-24.5)) )*10**(beta_over_five*M)
		return v
#========================================================================================

#========================================================================================
	def Pro_DensityFunction_QSO(self, M, z_l, m_l, M0, sigma):
		
		mu = M + 5*log10(((LightVelocity*10**6.0)/(10.*H_0)))
		
		if self.kappa_0 == 0.0:
			Value = (1./sqrt(2.*pi*sigma**2.))*(np.exp(-0.5*((M - M0)/sigma)**2.)) #*(4.*pi/3.)*((1/(1. + z_l))**3.)*10**(3.*(m_l - mu)/5.)  
			
		elif self.kappa_0 > 0.0:
			
			alfa_l = (sqrt(self.kappa_0)/(1 + z_l))*(10**((m_l - mu)/5.))
			Value = (1./sqrt(2.*pi*sigma**2.))*(np.exp(-0.5*((M - M0)/sigma)**2.)) #*(pi/self.kappa_0**(3./2))*(2*arcsin(alfa_l) - 2*alfa_l*sqrt(1 - alfa_l**2.))
			
		elif self.kappa_0 < 0.0 :
		  
			alfa_l = (sqrt(abs(self.kappa_0))/(1 + z_l))*(10**((m_l - mu)/5.))
			Value = (1./sqrt(2.*pi*sigma**2.))*(np.exp(-0.5*((M - M0)/sigma)**2.)) #*(pi/abs(self.kappa_0)**(3./2))*(2*alfa_l*sqrt(1 + alfa_l**2.) - 2*arcsinh(alfa_l))
		else:    
			pass
		  
		return Value


#========================================================================================
		
#========================================================================================
	def Pro_DensityFunction_SN(self, sy, cz, M0, s0, c0, rho, sigmas, sigmac,alfa, bbeta, z_l, m_l):
	  
		M = M0 - alfa*(sy - s0) + bbeta*(cz - c0)
		mu = M + 5.*log10(((LightVelocity*10**6.0)/(10.*H_0)))
		
		if self.kappa_0 == 0.0:

			v = (np.exp(-0.5*(((sy - s0)/sigmas)**2. + ((cz - c0)/sigmac)**2. - 2*rho*((sy - s0)*(cz - c0)/(sigmas*sigmac))))/(2*pi*sigmas*sigmac*sqrt(1-rho**2.))) #*(4.*pi/3.)*((1/(1. + z_l))**3.)*10**(3.*(m_l - mu)/5.) #*self.Heavisidefuction(m_l - m)   # rho = alfa/bbeta      
			
		elif self.kappa_0 > 0.0:
			
			alfa_l = (sqrt(self.kappa_0)/(1 + z_l))*(10**((m_l - mu)/5.))	# ((10*H_0)/(LightVelocity*10**6.0)) must be calculated in µ = M + 5*log10(((LightVelocity*10**6.0)/(10*H_0))) ce qui rend m - µ (au lieu de m - M) petit et du coup alfa_l < 1. (on fait ca parce que les formule de Volume sont déjà sans dimension).
			v = (np.exp(-0.5*(((sy - s0)/sigmas)**2. + ((cz - c0)/sigmac)**2. - 2*rho*((sy - s0)*(cz - c0)/(sigmas*sigmac))))/(2*pi*sigmas*sigmac*sqrt(1-rho**2.))) #**(pi/self.kappa_0**(3./2))*(2*arcsin(alfa_l) - 2*alfa_l*sqrt(1 - alfa_l**2.))    # rho = alfa/bbeta      
 
		elif self.kappa_0 < 0.0 :
		  
			alfa_l = (sqrt(abs(self.kappa_0))/(1 + z_l))*(10**((m_l - mu)/5.))
			v = (np.exp(-0.5*(((sy - s0)/sigmas)**2. + ((cz - c0)/sigmac)**2. - 2*rho*((sy - s0)*(cz - c0)/(sigmas*sigmac))))/(2*pi*sigmas*sigmac*sqrt(1-rho**2.))) #**(pi/abs(self.kappa_0)**(3./2))*(2*alfa_l*sqrt(1 + alfa_l**2.) - 2*arcsinh(alfa_l))    # rho = alfa/bbeta      
		else:    
			pass
		  
		return v


#========================================================================================
		
	def Poly1deg(self, a, b, x):
		v = a*x + b 
		return v
			
#========================================================================================
		
	def Poly2deg(self, a, b, c, x):
		v = a*x**2. + b*x + c
		return v
			
#========================================================================================

	def Poly4deg(self, a, b, c, d, e, x):
		v = a*x**4. + b*x**3. + c*x**2. + d*x + e
		return v
		
#========================================================================================

	def Poly8deg(self, a, b, c, d, e, f, g, h, i, x):
		v = a*x**8. + b*x**7. + c*x**6. + d*x**5. + e*x**4. + f*x**3. + g*x**2. + h*x + i
		return v
		
#========================================================================================

	def PolyNdeg(self, n, Ncoefficients, x):		# n = degree of polynom. Ncoefficients = list contain N=n+1 coefficents of this polynom
		v, pol = 0.0, ''
		for i in range(len(Ncoefficients)):
			if (n-i)>0:
				v += Ncoefficients[i]*x**(n-i)
				pol += str(Ncoefficients[i])+'x**'+str(n-i)+' + '
			elif (n-i)==0:
				v += Ncoefficients[i]*x**(n-i)
				pol += str(Ncoefficients[i])+'x**'+str(n-i)	
			elif (n-i)<0:
				pass
		#print pol
		return v
		
#========================================================================================

	def d_PolyNdeg(self, n, Ncoefficients, x):		# n = degree of polynom. Ncoefficients = list contain N=n+1 coefficents of this polynom
		v, pol = 0.0, ''
		for i in range(len(Ncoefficients)):
			if (n-i)>0:
				v += (n-i)*Ncoefficients[i]*x**((n-i)-1)
				pol += str(n-i)+'*'+str(Ncoefficients[i])+'x**'+str((n-i)-1)+' + '
			elif (n-i)==1:
				v += (n-i)*Ncoefficients[i]*x**((n-i)-1)
				pol += str(n-i)+'*'+str(Ncoefficients[i])+'x**'+str((n-i)-1)
			elif (n-i)==0:
				pass	
			elif (n-i)<0:
				print "STOP"
		#print pol
		return v
		
#========================================================================================

	def function_kcorr(self, alfa_nu, z):
		v = -2.5*(1 + alfa_nu)*np.log(1. + z) + 2.5*(1. + alfa_nu)*np.log(1. + 2.)
		return v
#========================================================================================

	def d_function_kcorr(self, alfa_nu, z):
		v = -2.5*(1. + alfa_nu)*(1./(1. + z))
		return v
#========================================================================================
	def Time_of_Run(self, NbOfGrid, NbOfSample):                                                     
		v = (NbOfGrid*NbOfSample*10)/(625*500)
		v2 = v/64.
		print "This run will take a time roughly: ", v, "min = ", v2, "h = ", v2/24., "days", "\n"

#========================================================================================
	def fact(self, n):
		i = 0
		nn = n - i
		result = 1
		while nn!=0:
			result = result*nn
			i = i + 1
			nn = n - i
		return result

#========================================================================================
#========================================================================================
	def sin_series(self, x, N):                 
		v=0
		for n in range(N):
			v = v + (((-1.)**n)/(self.fact(2*n + 1)))*x**(2*n + 1)
		return v

#========================================================================================
#========================================================================================
	def sinh_series(self, x, N):                 
		v=0
		for n in range(N):
			v = v + (1./(self.fact(2*n + 1)))*x**(2*n + 1)
		return v

#========================================================================================
