#!/usr/bin/python2.7.3
# -*- coding: latin-1 -*-

# Nom de fichier: Savedata.py

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
from scipy.stats import norm
import socket


import Methods
from Methods import Functions
from Cosmological_Model import Model
from HubbleDiagram import ApparentMagnitude
from GradientMethod import Gradient



class Save:
	"""Registration of all lists and arrays in one file"""
	def __init__(self, listofarrays, stringlistofarrays, filename): # filename has a forme: 'filename'

		# listofarrays is a list of several lists, stringlistofarrays is the list where contain a name for each list in listofarrays
		
		self.listofarrays = listofarrays
		self.stringlistofarrays = stringlistofarrays
		self.filename = filename

	def ReadWrite(self, header, private_pc):

		if private_pc == 'local':
			
			#path = '/home/dchbib/Substitute/Mixte_Plot/Image_Simulation/'+str(self.filename)
			#path = Name_of_machine+'PC_Dyaamanal_avant_formatage/toutes_dossier/Manal/Fit_Doudou/File/'+str(self.filename)
			path = '/home/mchbib/Bureau/Links/securite_site/images_data/'+str(self.filename)
			
		elif private_pc == 'work':
		  
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
			
			#path = Name_of_machine+'Zip_file/Python/Simulations/POO/OO/NullCorrelation/Mixte_Plot/Image_Simulation/'+str(self.filename)
			#path = Name_of_machine+'PC_Dyaamanal_avant_formatage/toutes_dossier/Manal/Fit_Doudou/File/'+str(self.filename)
			path = Name_of_machine+'Zip_file/Python/Simulations/POO/OO/SNIa/NullCorrelation/NullCorrelation_Simulation/NullCorrelation_Simulation_lifetime_code/NullCorrelation_Simulation_Like_QSOs_Code/Stability/Contours/'+str(self.filename)
			

				
		else:
			print "#################################*********#####################################", "\n"
			print "     YOU SHOULD DETRMINE THE TYPE OF YOUR COMPUTER AS: 'local' or 'work'", "\n"
			print "#################################*********#####################################", "\n"
		
		print "je suis la!, je commence à enregistrer les données! ", "\n"
		
		f = open(path, 'w')
		if header == True:

			#f.write("#")
			for l in range(len(self.stringlistofarrays)):

				listname = self.stringlistofarrays[l]
				listname = str(listname)
				f.write(listname)
				f.write("                ")

			f.write("\n")

		elif header == False:

			pass

		length = [len(a) for a in self.listofarrays]
		for i in range(max(length)):

			for l in range(len(self.listofarrays)):

				if i < len(self.listofarrays[l]):
					element = self.listofarrays[l][i]
					element = str(element)
					f.write(element)
					f.write("        ")
				else:
					f.write("Nan")
					f.write("        ")
					
			f.write("\n")			

				
		f.close()
		print "Voilà!, j'ai terminé d'enregistrer les données! Voici le chemin: ", "\n"
		print path, "\n"
		return path

