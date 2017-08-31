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

def pythonGSDE():
	"""pythonGSDE sovles equations (8), (29) and (33) from Shvets and Semenov J. Chem. Phys. 2013 paper. It is the first correction to the Ground-State Dominance Theory
	for finite chain effects, i.e. (N < infinity). GSDE (Ground-State Dominance Extended) is valid when h >> correlation length, h = plate separation distance, |deltaEnd| << Rg
    and |deltaSeg| << Rg:	
	
	USER INPUTS: There are two main places for users to input parameters, the first immediately follows the code header and is below. The second starts around lines 56 and deals
	with the separation lengths to be simulated. Need the minimum to be below the maximum peak position.
	
	"""
	
	#***************************  USER DEFINED INPUTS  ******************************#
	#********************************************************************************#
	lsegment  = 1.0			# The length of a statistical segment / math.sqrt(6)
	thetabulk = 1.0			# The volume fraction of the bulk	(Cb*lsegment**3)
	excvol    =	3.9 		# The excluded volume parameter     (v = (1-2*ChiPS)*lsegment**3)
	wparam    = 1.8			# The w parameter 					(w = lsegment**6)
	N         =	1.0			# Chain Length
	
	#************************** END USER DEFINED INPUTS  ****************************#
	#********************************************************************************#
	alpha     = 2. + (3.*excvol*(lsegment)**3)/(wparam*thetabulk)	# (2+3*v/w*Cb)
	print "alpha"
	print alpha
	A 		  = (4*(1+1/alpha))
	print"A"
	print A
	# The End Corrected Bulk Correlation Length. The Edwards Corrected Correlation Lenght from GSD.
	CorrLength= lsegment/math.sqrt(2*(excvol*thetabulk/lsegment**3 + wparam*thetabulk**2/lsegment**6 + 1/N))
	print "Correlation Length"
	print CorrLength
	Kappa     = 0.61617  
	# The I(alpha) function defined in A.A. Shvets and A.N. Semenov J. Chem. Phys. 2013
	fnx1      = (1-(1/math.sqrt(1+alpha)))*np.log((1+math.sqrt(1+alpha))/math.sqrt(alpha))-(np.log(2)/math.sqrt(1+alpha))+(2/math.sqrt(1+alpha)*np.log(math.sqrt(1+alpha)+1-math.sqrt(alpha)))
	deltaEnd  = (math.sqrt(6)*lsegment**4/(thetabulk*math.sqrt(wparam)))*fnx1
	print "deltaEnd"
	print deltaEnd
	deltaSeg  = (-1*math.sqrt(6)*lsegment**4)/(thetabulk*math.sqrt(wparam))*np.log((1+math.sqrt(1+alpha))/math.sqrt(alpha))
	print "deltaSeg"
	print deltaSeg
	
	#***************************  USER DEFINED INPUTS  ******************************#
	#********************************************************************************#
	# Separation Distance is h, Maximum Lenght is L
	refRg = lsegment*N**(0.5)
	L = 8
	# L2 Controls the maximum separation distance
	L2 = L/refRg
	# L1 Controls the minimum separation distance
	L1 = 0.4
	#************************** END USER DEFINED INPUTS  ****************************#
	#********************************************************************************#
	#Checks on validity
	check1 = deltaEnd/refRg
	check2 = deltaSeg/refRg
	print "CHECK 1: |deltaEnd/Rg| << 1"
	print check1
	print "CHECK 2: |deltaSeg/Rg| << 1"
	print check2
		
	h = np.linspace(L1,L2,400)
	n = np.linspace(1,10000,10000)

	WTot = []
	WGSD = []
	WEND = []
	SumTot = []
	hList = []
	wEND = 0.
	
	for i in h: # Separation loop 
	
		wGSD = (-1*N**(0.5)/CorrLength)*thetabulk*A**2*np.exp(-1*i/CorrLength*lsegment*N**(0.5))
		preFactor = (4*thetabulk/N/lsegment**3)*(deltaEnd-(i*lsegment*N**0.5+2*deltaSeg)/2*np.log(1+(2*deltaEnd/(i*lsegment*N**0.5+2*deltaSeg))))
		
		for j in n: # Summation Loop
			y = 4*(math.pi**2)*(j**2)/(i**2)
			wEND = wEND + ((1-(1+y)*np.exp(-1*y))/(y-1+np.exp(-1*y)))
			
		SumTot.append(2*wEND)
		wEND = preFactor*(2*wEND - Kappa*i + 1)
		WTot.append((wGSD + wEND)*math.sqrt(N)*lsegment**3/thetabulk)
		hList.append(i)
		WGSD.append(wGSD)
		WEND.append(wEND)
		wEND = 0.0
		
		
	lenWTot = len(WTot)
	cnt = 0
	flag = 0
	maxpts = []
	hmaxpts = []
	#****************************  PMF MAXIMUM SEARCH  ******************************#
	#********************************************************************************#
	# Maximum Location and Height
	while cnt < len(WTot)-2:
		h1 = hList[cnt]
		d1 = WTot[cnt]
		h2 = hList[cnt+1]
		d2 = WTot[cnt+1]
		h3 = hList[cnt+2]
		d3 = WTot[cnt+2]
		
		if cnt > 0 and cnt < lenWTot-2 and d2 > d1 and d3 < d2: # the maximum
			maxpts = np.append(maxpts,d1)
			maxpts = np.append(maxpts,d2)
			maxpts = np.append(maxpts,d3)
			hmaxpts = np.append(hmaxpts,h1)
			hmaxpts = np.append(hmaxpts,h2)
			hmaxpts = np.append(hmaxpts,h3)
			flag = 1
			cnt = cnt + 1
		else:
			
			cnt = cnt + 1
	
	if flag == 1:
		coef = np.zeros(3)
		# Interpolation to find the maximum 
		MaxPeakPolyCoef = np.polyfit(hmaxpts[0:3],maxpts[0:3],2,full=True) # fit 2nd order polynomial  
		coef1 = MaxPeakPolyCoef[0]		# second order polynomial fit coefficients from x^0 + x^1 + x^2
		coef[2] = coef1[0]
		coef[1] = coef1[1]
		coef[0] = coef1[2]
		coef2 = np.polyder(coef1,1)		# first order derivative
		roots = np.roots(coef2) 		# finds the root of the polynomial
		MaxPeakPos = roots
		MaxPeak = np.polynomial.polynomial.polyval(roots,coef)	# finds the maximum
		
		print "MaxPeakPos"
		print MaxPeakPos
		print "MaxPeak"
		print MaxPeak
		check3 = MaxPeakPos/CorrLength
		print "CHECK 3: |MaxPeakPos/CorrLength| >> 1"
		print check3
	
	else:
		
		print"no maximum found"
	
	#********************************  PMF PLOTTING  ********************************#
	#********************************************************************************#
	hList = np.asarray(hList)
	WGSD  = np.asarray(WGSD)
	WEND  = np.asarray(WEND)
	WTot  = np.asarray(WTot)
	
	AlignedGSD = np.column_stack((hList,WGSD))
	AlignedEND = np.column_stack((hList,WEND))
	AlignedGSDE = np.column_stack((hList,WTot))
	AlignedSumTot = np.column_stack((hList,SumTot))
	
	np.savetxt("GSD.txt", AlignedGSD)
	np.savetxt("END.txt", AlignedEND)
	np.savetxt("GSDE.txt", AlignedGSDE)
	np.savetxt("SumTot.txt", AlignedSumTot)
	
	plt.plot(hList,WTot,'-')
	plt.title('GSDE')
	plt.ylabel('W/kbT/Rg/Cb/N')
	plt.xlabel("h/Rg")
	plt.savefig('PMF_GSDE')
	plt.close()
	
	plt.plot(hList,WGSD,'-')
	plt.title('GSDE')
	plt.ylabel('W/kbT/Rg/Cb/N')
	plt.xlabel("h/Rg")
	plt.savefig('PMF_GSD')
	plt.close()
	
	plt.plot(hList,WEND,'-')
	plt.title('WEND')
	plt.ylabel('W/kbT/Rg/Cb/N')
	plt.xlabel("h/Rg")
	plt.savefig('PMF_WEND')
	plt.close()
	
	