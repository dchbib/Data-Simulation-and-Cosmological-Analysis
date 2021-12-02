#!/usr/bin/python2.7.3
# -*- coding: latin-1 -*-

#Nom de fichier: Cosmological_Model.py
#==============================================================================
#title           :Cosmological_Model.py
#description     :This file gives us the comoving distance, volume and luminosity distance.
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


H_0 = 70.0 #Km.s-1.Mpc-1
LightVelocity = 2.99792458*10**5.0 # Km.s-1
toto = 0.0
#bbeta = (3*np.log(10))/5.

class Model:
	"""Cosmological Model correspond to given parameters"""
	def __init__(self, lambda_0, omega_0):

		self.param = Functions(lambda_0, omega_0)
		
	def curvature(self, z, lambda_0, omega_0, M):

		zz = self.param.LowerBound(z)
		comobileDistance = self.param.Trapeze(self.param.FunctiontoIntegrate, zz, 1.0)
		lookback_time = self.param.Trapeze(self.param.FunctiontoIntegrate2, zz, 1.0)
		if self.param.kappa_0 > 0.0:

			epsilon = comobileDistance*sqrt(self.param.kappa_0)
			volume = pi*(2.0*epsilon - sin(2.0*epsilon))/(self.param.kappa_0**(3./2.))    # dimensionless volume
			#volume = volume*(((LightVelocity*10**6.0)/H_0)**3.0)	# volume en [pc]**3
			#In the case where kappa0 > 0, we use this condition (angle changed by (2*pi - angle) or sin(x) by abs(sin(x)) because the sin function gives us a negatif value with angles > pi (redshift of atipodes). If we do not do it, it dangerous to obtain a negative values ​​for D_L  (changing of tan(x) by : abs(sin(x))/cos(x))
			d_L = (LightVelocity*10**6.0)*(1. + z)*abs(sin(epsilon))/(H_0*sqrt(self.param.kappa_0)) #je multiplie par 10**6 parce que H_0 est en Km.s-1.Mpc-1
				

		elif self.param.kappa_0 < 0.0:

			epsilon = comobileDistance*sqrt(-self.param.kappa_0)
			volume = pi*(sinh(2.0*epsilon) - 2.0*epsilon)/((-self.param.kappa_0)**(3./2.)) # dimensionless volume
			#volume = volume*(((LightVelocity*10**6.0)/H_0)**3.0)	# volume en [pc]**3
						
			d_L = (LightVelocity*10**6.0)*(1. + z)*sinh(epsilon)/(H_0*sqrt(-self.param.kappa_0)) #je multiplie par 10**6 parce que H_0 est en Km.s-1.Mpc-1
			#d_L = (1. + z)*sinh(epsilon)/(sqrt(-self.param.kappa_0)) #je multiplie par 10**6 parce que H_0 est en Km.s-1.Mpc-1

		elif self.param.kappa_0 == 0.0:

			volume = 4.*pi*(comobileDistance**3.0)/3.  # dimensionless volume
			#volume = volume*(((LightVelocity*10**6.0)/H_0)**3.0)	# volume en [pc]**3
			
			d_L = (LightVelocity*10**6.0)*(1. + z)*comobileDistance/H_0  #je multiplie par 10**6 parce que H_0 est en Km.s-1.Mpc-1

		else:

			pass
		      
		     
		return comobileDistance, volume, d_L, lookback_time
		
	""" 
	def Volume_series(self, redshift):
	  
		zz = self.param.LowerBound(redshift)
		tau = self.param.Trapeze(self.param.FunctiontoIntegrate, zz, 1.0)
		
		if self.param.kappa_0 > 0.0:
			VV = (pi/(self.param.kappa_0**1.5))*(2.*tau*sqrt(self.param.kappa_0) - self.param.sin_series(2.*tau*sqrt(self.param.kappa_0), 30))
		elif self.param.kappa_0 < 0.0:
			VV = (pi/(abs(self.param.kappa_0)**1.5))*(2.*tau*sqrt(abs(self.param.kappa_0)) - self.param.sinh_series(2.*tau*sqrt(abs(self.param.kappa_0), 30)))   
		elif self.param.kappa_0 == 0.0:
			
			VV = (4.*pi*tau**3.)/3.   
		else:
			pass
		      
   .    	return VV
	"""
	def Pro_DensityFunction(self, M, z_l, m_l, M0, sigma):
	  
		zz = self.param.LowerBound(z_l)
		comobileDistance = self.param.Trapeze(self.param.FunctiontoIntegrate, zz, 1.0)
	  
		if self.param.kappa_0 == 0.0:
			Value = (1./sqrt(2.*pi*sigma**2.))*(np.exp(-0.5*((M - M0)/sigma)**2.))*(4.*pi/3.)*(comobileDistance**3.0)
			
		elif self.param.kappa_0 > 0.0:
			
			epsilon = comobileDistance*sqrt(self.param.kappa_0)
			      
			Value = (1./sqrt(2.*pi*sigma**2.))*(np.exp(-0.5*((M - M0)/sigma)**2.))*(pi/self.param.kappa_0**(3./2))*(2.0*epsilon - abs(sin(2.0*epsilon)))
			
		elif self.param.kappa_0 < 0.0 :
		  
			epsilon = comobileDistance*sqrt(abs(self.param.kappa_0))
			
			Value = (1./sqrt(2.*pi*sigma**2.))*(np.exp(-0.5*((M - M0)/sigma)**2.))*(pi/abs(self.param.kappa_0)**(3./2))*(sinh(2.0*epsilon) - 2.0*epsilon)
		else:    
			pass
		  
		return Value


