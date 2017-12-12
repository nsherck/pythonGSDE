#******************General Imports**********************************
import os
import sys
import subprocess
import numpy as np
import scipy as sp
import shutil
import datetime
import math
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.integrate import quad

# The definition of the integrand for fluctuations
def integrand(x, m): 
	return x**2*math.log((2/(x*math.pi/m)**2)*(x*math.pi/m - 1 + np.exp(-1*x*math.pi/m)))
	
def chk():
	x = np.linspace(1000000,100000000,100) # separation
	m = np.linspace(1,1000,1000) 	# Integers
	Rg = 1
	
	
	sumTot = []
	hList  = []
	sumReg = []
	sum1 = 0
	for i in x:
		preFactFluct = math.pi**2
		#I2 = quad(integrand, 1, np.inf, args=(i))
		I2 = [0]
		#print "I2"
		#print I2[0]
		
		for j in m:
			sum1 = sum1 + integrand(j,i) 
			#print sum1
			
		sumReg.append(preFactFluct*(sum1 - I2[0]))	
		sumTot.append(preFactFluct*sum1)
		hList.append(i)

	
	hList = np.asarray(hList)
	sumTot  = np.asarray(sumTot)
	sumReg  = np.asarray(sumReg)
	
	AlignedSumTot = np.column_stack((hList,sumTot))
	AlignedSumReg = np.column_stack((hList,sumReg))

	np.savetxt("SumTotChk.txt", AlignedSumTot)
	np.savetxt("SumRegChk.txt", AlignedSumReg)
	
	plt.plot(hList,sumTot,'-')
	plt.title('GSDE')
	plt.ylabel('W/kbT/Rg/Cb/N')
	plt.xlabel("h/Rg")
	plt.savefig('SumTotChk')
	plt.close()
	
	plt.plot(hList,sumReg,'-')
	plt.title('GSDE')
	plt.ylabel('W/kbT/Rg/Cb/N')
	plt.xlabel("h/Rg")
	plt.savefig('SumRegChk')
	plt.close()
	
	