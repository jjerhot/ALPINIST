#universal functions definitions

import os
import numpy as np
from sklearn.neighbors import KernelDensity
from . import exotic_constants as c
import random

def export(experiment,fileName,data):
	output_dir = os.getcwd()+'/../tab_prod/'
	outPath = output_dir + experiment + '/'

	np.savetxt(outPath + fileName,data,fmt='%.3e')
	print('\n[Info:] \t', 'File ' + fileName + ' saved to ' + outPath)
	return

def lambda_Kallen(xx,yy,zz): #Kallen function
    return (xx**2-(yy+zz)**2)*(xx**2-(yy-zz)**2)

def lorentz_transf(gam,bet): #returns lorentz transformation matrix for given gamma and beta
    return [[gam,gam*bet[0],gam*bet[1],gam*bet[2]],
            [gam*bet[0],1+(gam-1)*bet[0]**2/np.linalg.norm(bet)**2,(gam-1)*bet[0]*bet[1]/np.linalg.norm(bet)**2,(gam-1)*bet[0]*bet[2]/np.linalg.norm(bet)**2],
            [gam*bet[1],(gam-1)*bet[1]*bet[0]/np.linalg.norm(bet)**2,1+(gam-1)*bet[1]**2/np.linalg.norm(bet)**2,(gam-1)*bet[1]*bet[2]/np.linalg.norm(bet)**2],
            [gam*bet[2],(gam-1)*bet[2]*bet[0]/np.linalg.norm(bet)**2,(gam-1)*bet[2]*bet[1]/np.linalg.norm(bet)**2,1+(gam-1)*bet[2]**2/np.linalg.norm(bet)**2]]

def beta_lor(p,E):
	return p/E

def gamma_lor(beta):
	return 1/np.sqrt(1-beta**2)
        
def f_kde(kde,norm,e,th): #returns probability for given energy and angle
	val = 0
	if e>0 and th>0:
		val = norm * np.exp(kde.score_samples([[np.log(e),np.log(th)]]))/(e*th)
	return val

def decay_meson(p_cm,m_exo,lor): #returns energy and angle of exotic in random decay m_1 -> m_2 + exotic
	th_exo_cm = np.arccos(random.uniform(-1, 1))
	phi_exo_cm = random.uniform(0, 2*np.pi)

	p_lor_cm = [np.sqrt(p_cm**2+m_exo**2),p_cm*np.sin(th_exo_cm)*np.cos(phi_exo_cm),p_cm*np.sin(th_exo_cm)*np.sin(phi_exo_cm),p_cm*np.cos(th_exo_cm)]
	p_lor_exo = np.dot(lor,p_lor_cm)

	#return e_exo and th_exo
	return p_lor_exo[0], np.arccos(p_lor_exo[3]/np.sqrt(p_lor_exo[1]**2+p_lor_exo[2]**2+p_lor_exo[3]**2))