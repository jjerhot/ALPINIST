
def dP(self,m_ds,z,pt2):
	gSNN = 1.2E-3
	return gSNN**2/(8*np.pi**2) * z * (c.m_p**2*(2-z)**2 + pt2)/(c.m_p**2*z**2+m_ds**2*(1-z)+pt2)**2

def dsigma(self,p_ds,th_ds,m_ds,e_beam,sigma):
	gSNN = 1.2E-3
	e_ds = np.sqrt(p_ds**2+m_ds**2)
	z_ds = p_ds/setup.p_beam[self.exp]
	pt = p_ds*np.sin(th_ds)
	if pt > 1.: return 0. # restriction from 2008.08108
	pp2 = c.m_p**2 + m_ds**2 -2*(e_ds*e_beam + p_ds*setup.p_beam[self.exp]*np.cos(th_ds))
	jacobian = e_ds * p_ds * np.sin(2*th_ds) / setup.p_beam[self.exp]
	return jacobian * sigma * gSNN**2 * setup.coupling_ref**setup.coupling_exp * z_ds * (c.m_p**2 * (2-z_ds)**2 + pt**2) / (8*np.pi**2 * (c.m_p**2 * z_ds**2 + m_ds**2*(1-z_ds) + pt**2)**2)
	# return jacobian * sigma * gSNN**2 * setup.coupling_ref**setup.coupling_exp * (self.form_factor_s(m_ds)*self.form_factor_ppstar(pp2))**2 * z_ds * (c.m_p**2 * (2-z_ds)**2 + pt**2) / (8*np.pi**2 * (c.m_p**2 * z_ds**2 + m_ds**2*(1-z_ds) + pt**2)**2)

def dsigma_ptz(self,z_ds,pt2_ds,m_ds):
	gSNN = 1.2E-3
	boundary1 = 4*setup.p_beam[self.exp]**2 * z_ds * (1-z_ds)**2
	boundary2 = 4*setup.p_beam[self.exp]**2 * z_ds * (1-z_ds)
	boundary3 = 4*setup.p_beam[self.exp]**2 * (1-z_ds)**2 / z_ds
	# boundary1 = 0.1*setup.p_beam[self.exp]**2 * z_ds * (1-z_ds)**2
	# boundary2 = 0.1*setup.p_beam[self.exp]**2 * z_ds * (1-z_ds)
	# boundary3 = 0.1*setup.p_beam[self.exp]**2 * (1-z_ds)**2 / z_ds
	if boundary1/100 < pt2_ds or boundary2/100 < m_ds**2 or boundary3/100 < c.m_p**2: return 0. # ensure <<

	boundary = (m_ds**2 * (1-z_ds) + c.m_p**2 * z_ds**2 + pt2_ds)/(4 * setup.p_beam[self.exp]**2 * z_ds * (1-z_ds)**2)
	if boundary >= 0.1: return 0.

	prob = gSNN**2 / (8 * np.pi**2) * z_ds * (c.m_p**2 * (2 - z_ds)**2 + pt2_ds) / (c.m_p**2 * z_ds**2 + m_ds**2 * (1 - z_ds) + pt2_ds)**2
	# prob = gSNN**2 / (8 * np.pi**2) * z_ds * (c.m_p**2 * (2 - z_ds)**2 + pt2_ds) / (c.m_p**2 * z_ds**2 + m_ds**2 * (1 - z_ds) + pt2_ds)**2
	s_prime = 2*c.m_p*setup.p_beam[self.exp]*(1-z_ds) + 2*c.m_p**2
	if prob < 0: return 0.
	return bf.brems_cross_section().SigmaPP(s_prime) * prob

def dsigma_eth(self,e_ds,th_ds,m_ds):
	if e_ds < m_ds: return 0.
	p_ds = np.sqrt(e_ds**2 - m_ds**2)
	z_ds = p_ds/setup.p_beam[self.exp]
	pt2_ds = (p_ds*np.sin(th_ds))**2
	return self.dsigma_ptz(z_ds,pt2_ds,m_ds)

def sigma_int(self,m_ds):
	pt2_list = np.linspace(0,50**2,1001).tolist()
	dpt2 = pt2_list[1]-pt2_list[0]
	# print('pt2_list[0]:',pt2_list[0],'pt2_list[-1]',pt2_list[-1],'dpt2',dpt2)
	z_list = np.linspace(0.1,0.9,91).tolist()
	dz = z_list[1]-z_list[0]
	# print('z_list[0]:',z_list[0],'z_list[-1]',z_list[-1],'dz',dz)
	binsize = dpt2*dz
	print('m_ds',m_ds,'binsize:',binsize)
	sig = 0
	for z in z_list:
		p = z*setup.p_beam[self.exp]
		for pt2 in pt2_list:
			if pt2 < p**2:
				sig += self.dsigma_ptz(z,pt2,m_ds)
	sig *= binsize
	return sig

def sigma_int_eth(self,m_ds):
    
	de = self.en_list[1]-self.en_list[0]
	# print('self.en_list[0]:',self.en_list[0],'self.en_list[-1]',self.en_list[-1],'de',de)
	dth = self.th_list[1]-self.th_list[0]
	# print('self.th_list[0]:',self.th_list[0],'self.th_list[-1]',self.th_list[-1],'dth',dth)
	binsize = de*dth
	print('m_ds',m_ds,'binsize:',binsize)
	sig = 0
	for e_ds in self.en_list:
		if e_ds>m_ds: p_ds = np.sqrt(e_ds**2-m_ds**2)
		z_ds = p_ds/setup.p_beam[self.exp]
		# if z_ds > 0.1 and z_ds < 0.9 and e_ds > m_ds and e_ds != self.en_list[0]:
		if e_ds > m_ds and e_ds != self.en_list[0]:
		# sigmaPP = bcsTot.SigmaPP(2*c.m_p*(setup.p_beam[self.exp]-p_ds+c.m_p))
		# data_list_sublist = []
		#skip first energy bin to satisfy WW approx
		# if e_ds > m_ds and e_ds != self.en_list[0]:
			for th_ds in self.th_list:
				sig += self.dsigma_eth(e_ds,th_ds,m_ds)
		# 		pt = p_ds*np.sin(th_ds)
		# 		pp2 = c.m_p**2 + m_ds**2 -2*(e_ds*e_beam + p_ds*setup.p_beam[self.exp]*np.cos(th_ds))
		# 		jacobian = 1 # e_ds * p_ds * np.sin(2*th_ds) / setup.p_beam[self.exp]
		# 		dsigma = jacobian * sigmaPP * gSNN**2 * setup.coupling_ref**setup.coupling_exp * (self.form_factor_scalar(m_ds)*self.form_factor_ppstar(pp2))**2 * z_ds * (c.m_p**2 * (2-z_ds)**2 + pt**2) / (8*np.pi**2 * (c.m_p**2 * z_ds**2 + m_ds**2*(1-z_ds) + pt**2)**2)
		# 		data_list_sublist.append([th_ds, e_ds, m_ds, dsigma])
		# else:
		# 	for th_ds in self.th_list:
		# 		data_list_sublist.append([th_ds, e_ds, m_ds, 0.])

	# for z in z_list:
	# 	p = z*setup.p_beam[self.exp]
	# 	for pt2 in pt2_list:
	# 		if pt2 < p**2:
	# 			sig += self.dsigma_ptz(z,pt2,m_ds)
	sig *= binsize
	return sig


def ds_brems_process(self,m_ds):
	bcsTot = bf.brems_cross_section()
	gSNN = 1.2E-3
	data_list = []
	p_ds = 0
	e_beam = np.sqrt(setup.p_beam[self.exp]**2 + c.m_p**2)
	for e_ds in self.en_list:
		if e_ds>m_ds: p_ds = np.sqrt(e_ds**2-m_ds**2)
		z_ds = p_ds/setup.p_beam[self.exp]
		sigmaPP = bcsTot.SigmaPP(2*c.m_p*(setup.p_beam[self.exp]-p_ds+c.m_p))
		data_list_sublist = []
		#skip first energy bin to satisfy WW approx
		if e_ds > m_ds and e_ds != self.en_list[0]:
			for th_ds in self.th_list:
				pt = p_ds*np.sin(th_ds)
				pp2 = c.m_p**2 + m_ds**2 -2*(e_ds*e_beam + p_ds*setup.p_beam[self.exp]*np.cos(th_ds))
				jacobian = 1 # e_ds * p_ds * np.sin(2*th_ds) / setup.p_beam[self.exp]
				dsigma = jacobian * sigmaPP * gSNN**2 * setup.coupling_ref**setup.coupling_exp * (self.form_factor_scalar(m_ds)*self.form_factor_ppstar(pp2))**2 * z_ds * (c.m_p**2 * (2-z_ds)**2 + pt**2) / (8*np.pi**2 * (c.m_p**2 * z_ds**2 + m_ds**2*(1-z_ds) + pt**2)**2)
				data_list_sublist.append([th_ds, e_ds, m_ds, dsigma])
		else:
			for th_ds in self.th_list:
				data_list_sublist.append([th_ds, e_ds, m_ds, 0.])

		data_list.append([data_list_sublist])

	return data_list
