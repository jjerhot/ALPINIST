# modes/dp_brems_functions.py
#
# brems_cross_section written by Tom
# methods for DP brems production, should be documented

import numpy as np
from ALP_production.general import exotic_constants as c

class brems_cross_section:
	"""Utility class for implementation of modified Weizssacker Williams approximation"""
	def __init__(self): # set default values
		self.mprime = 0.5 # GeV
		self.M = 0.98 # GeV
		self.epsilon = 1.E-6
		self.Ep = 70 # GeV

	def setMPrime(self,val):
		self.mprime = val

	def setM(self,val):
		self.M = val
	    
	def setEpsilon(self,val):
		self.epsilon = val
	    
	def setMomentum(self,val):
		self.Ep = val

	def evaluate(self, x, p):#   returns d^2N/dE'dPt(z; pt, m)
		z = x[0] 
		pt = p[0] 
		mp = p[1] 
		self.setMPrime(mp) 
		sPrime = (self.M+self.Ep)*(self.M+self.Ep)*(1.-z) 
		s = 2.*self.M*self.Ep 
		return (self.wba(pt,z)*self.SigmaPP(sPrime)/self.SigmaPP(s)/self.Ep) 

	def evaluate2d(self, x, p):#   returns d^2N/dE'dPt(z,pt; m)
		z = x[0] 
		pt = x[1] 
		mp = p
		self.setMPrime(mp) 
		sPrime = (self.M+self.Ep)*(self.M+self.Ep)*(1.-z) 
		s = 2.*self.M*self.Ep 
		return (self.wba(pt,z)*self.SigmaPP(sPrime)/self.SigmaPP(s)/self.Ep) #dimension of E^-2 

	def evaluateIntProb(self, x, p):#   return the total integrated probability of A' emission in the range [zmin,zmax] x [0,ptmax] = f(zmax; ptmax, mp, zmin)
		zmax  = x[0] 
		zbins   = 100
		ptmax = p[0] 
		mp = p[1] 
		zmin = p[2] 

		self.setMPrime(mp) 
		s = 2.*self.M*self.Ep 
		zwid = (zmax-zmin)/zbins 
		wtotal = 0. 
		for i in range(zbins+1):
			z = zmin + zwid*i 
			sPrime = (self.M+self.Ep)*(self.M+self.Ep)*(1.-z) 
			wtotal += (self.wbaInt(ptmax,z)*self.SigmaPP(sPrime)/self.SigmaPP(s)/self.Ep) # dimension of Energy^-1
		return wtotal*zwid*self.Ep 

	def evaluateInt(self, x, p):#   return dN/dE', i.e., the A' emission probability integrated up to a given ptmax = f(z; ptmax, mp)
		z = x[0] 
		ptmax = p[0] 
		mp = p[1] 
		self.setMPrime(mp) 
		sPrime = (self.M+self.Ep)*(self.M+self.Ep)*(1.-z) 
		s = 2.*self.M*self.Ep 
		return (self.wbaInt(ptmax,z)*self.SigmaPP(sPrime)/self.SigmaPP(s)/self.Ep) # dimension of Energy^-1

	def evaluateXSInt(self, x, p):#   return dsigma/dz, i.e., the A' emission diff. xs integrated up to a given ptmax = f(z; ptmax, mp)
		z = x[0] 
		ptmax = p[0] 
		mp = p[1] 
		self.setMPrime(mp) 
		sPrime = (self.M+self.Ep)*(self.M+self.Ep)*(1.-z) 
		s = 2.*self.M*self.Ep 
		return (self.wbaInt(ptmax,z)*self.SigmaPP(sPrime)/self.SigmaPP(s)*7000.) # dimension of microbarn, multiply by the elastic pp xs (??)

	def evaluateSigmaPP(self, x, p):#   return inelastic xs as f(s')
		sPrime = x[0] 
		return self.SigmaPP(sPrime) 

	def evaluateWba(self, x, p):#   returns wba(z; pt, mp)
		z = x[0] 
		pt = p[0] 
		mp = p[1] 
		self.setMPrime(mp)
		return self.wba(pt,z) 

	def evaluateWbaInt(self, x, p):#   returns integral of wba as f(z; ptmax, mp)
		z = x[0] 
		ptmax = p[0] 
		mp = p[1] 
		self.setMPrime(mp) 
		return self.wbaInt(ptmax,z) 

	def H(self, pt, z):
		return pt*pt + (1.-z)*self.mprime*self.mprime + z*z*self.M*self.M 

	def SigmaPP(self, sPrime):
		Z = 35.45 # mb
		B = 0.308 # mb
		Y1= 42.53 # mb
		Y2= 33.34 # mb
		s0 = 5.38*5.38 # GeV^2
		s1 = 1.0  # GeV^2
		eta1 = 0.458 
		eta2 = 0.545 
		return Z + B*(np.log(sPrime/s0)*np.log(sPrime/s0)) + Y1*np.power(s1/sPrime,eta1) - Y2*np.power(s1/sPrime,eta2)

	def wba(self, pt, z):# return wba such that SigmaPP(s') dz dsigmapt gives the differential cross section 
		alpha = 7.297E-3 
		acca = self.H(pt,z) 
		return self.epsilon*alpha*0.5/(np.pi*acca)*(
					(1.+(1.-z)*(1.-z))/z
					- 2.*z*(1.-z)*((2.*self.M*self.M+self.mprime*self.mprime)/acca -2*np.power(z*self.M*self.M/acca,2))
					+ 2.*z*(1.-z)*(1+(1.-z)*(1.-z))*np.power(self.M*self.mprime/acca,2)
					+ 2.*z*(1.-z)*(1.-z)*np.power(self.mprime**2/acca,2)
				)*2.*pt 

	def wbaInt(self, ptmax, z):#   return integral of wba(z,pt^2) on p_t^2, needs to be integrated in dz
		alpha = 7.297E-3 
		A = (1.-z)*self.mprime*self.mprime+z*z*self.M*self.M 
		return self.epsilon*alpha*0.5/np.pi*(
					(1.+(1.-z)*(1.-z))/z*np.log(1.+ptmax*ptmax/A)
					- 2.*z*(1.-z)*(2.*self.M*self.M+self.mprime*self.mprime)*ptmax*ptmax/(A*(ptmax*ptmax+A))
					+ 2.*z*(1.-z)*(2.*np.power(z*self.M*self.M,2) + (1.+(1.-z)*(1.-z))*self.M*self.M*self.mprime*self.mprime + (1.-z)*np.power(self.mprime,4))*
					ptmax*ptmax*(ptmax*ptmax+2*A)/(2.*np.power(ptmax*ptmax+A,2)*A*A)
				)