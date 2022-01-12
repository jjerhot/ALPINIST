#universal functions definitions

import numpy as np
import alp_constants as c

def alpha_s(m_a): #leading order alpha_s running
	nf = 6
	beta_0 = 11. - 2/3.*nf

	Lambda_QCD = 0.34

	if m_a <= 1:
		return 1
	elif m_a > 1 and m_a < 1.5: #cubic fit
		return 5.2377 * m_a**3 - 19.913 * m_a**2 + 24.1129 * m_a - 8.4376

	else:
		return 4*np.pi/(beta_0*np.log(m_a**2 /Lambda_QCD**2 ))

	print("Error: alpha_s function - something went wrong")
	return 0

def lambda_Kallen(xx,yy,zz): #Kallen function
    return (xx**2-(yy+zz)**2)*(xx**2-(yy-zz)**2)