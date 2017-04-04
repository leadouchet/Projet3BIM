# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import re




def readparam(param, setting) : 
	file = param
	data=pd.read_csv(file, sep=' ',header=0, index_col = 0)
	parameters = {}
	n = data.shape[0]
	for i in range(n) :
		c = data.iloc[i, setting]
		if (re.search('[a-zA-Z]',c) == None) : 
			   parameters[data.index[i]] = float(c)
		else : 
			c = re.split("\*",c)
			print c		
			if (len(c)== 1) :
				parameters[data.index[i]] = parameters[c[0]]
			else : 
				parameters[data.index[i]] = float(c[0]) * parameters[c[1]]




readparam('param.csv', 1)

def CircadianRythme(t, initial_conditions, param) :
	"""This function allows to simulate the cycle of proteins 
	during the circadian cycle. The inital_condition must be into a 
	in the order shows below, 
	param is a dictonnary containing as keys the name of the 
	parameter. dt and Run_time are numbers respectively describing 
	the integration step and the  Running time. The function return a list 
	with de dQuantites for one step

	The model has been built by Jean-Christophe Leloup, Albert Goldbeter
	in their article : Modeling the mammalian circadian clock :Sensitivity 
	analysis and multiplicity of oscillatory mechanisms.
	"""
#--------------------------
# Parameters importation : 
#--------------------------




#-----------------------
# Initial conditions  : 
#-----------------------

	# mRNAs of per, Cry and Bmal : 
	Mp = initial_conditions[0]
	Mc = initial_conditions [1]
	Mb = initial_conditions[2]

	# Phosporylated and non-phosphorylated proteins PER
	# and Cry in the cytosol : 

	Pc = initial_conditions[3]
	Bn = initial_conditions[4]
	Cc = initial_conditions[5]
	Pcp = initial_conditions[6]
	Ccp = initial_conditions[7]

	# Phosporylated and non-phosphorylated PER- Cry complexe
	# in the cytosol and nucleus : 

	PCc = initial_conditions[8]
	Pcn = initial_conditions[9]
	PCcp = initial_conditions[10]
	PCnp = initial_conditions[11]

	# Phosphorylated and non-phosphorylated protein BMAL1 in
	# the cytosol and nucleus : 

	Bc = initial_conditions[12]
	Bcp = initial_conditions[13]
	Bn = initial_conditions[14]
	Bnp = initial_conditions[15]

	# Inactive complex between PER-CRY and CLOCK-BMAL1 in 
	# nucleus : 

	In = initial_conditions[16]

#--------------
# Parameters : 
#--------------


	# Rate constants for modification : 

	k1 = param['k1']
	k2 = param['k2']
	k3 = param['k3']
	k4 = param['k4']
	k5 = param['k5']
	k6 = param['k6']
	k7 = param['k7']
	k8 = param['k8']

	# Activation constant

	Kap = param['Kap']
	Kac = param['Kac']
	Kib = param['Kib']

	#Nonspecific degradation rate constant

	kdmb = param['kdmb']
	kdmc = param['kdmc']
	kdmp = param['kdmp']
	kdnc = param['kdnc']
	kdn = param['kdn']

	# Michaelis constant : 

	Kd = param['Kd']
	Kdp = param['Kdp']
	Kp = param['Kp']
	Kmb = param['Kmb']
	Kmc = param['Kmc']
	Kmp = param['Kmp']

	#Rate constant for synthesis : 

	kstot = param['kstot']
	ksb = param['ksb']
	ksc = param['ksc']
	ksp = param['ksp']

	# Degree of cooperativity : 

	n = param['n']
	m = param['m']

	#Phosphorylation rate : 

	Vphos = param['Vphos']
	#Maximum Rate : 

	V1b = param['V1b']
	V1c = param['V1c']
	V1p = param['V1p']
	V1pc = param['V1pc']
	V2b = param['V2b']
	V2c = param['V2c']
	V2pc = param['V2pc']
	V3b = param['V3b']
	V3pc = param['V3pc']
	V4b = param['V4b']
	V4pc = param['V4pc']

	#Maximum rate of degradation

	Vdbc = param['Vdbc']
	Vdbn = param['Vdbn']
	Vdcc = param['Vdcc']
	Vdin = param['Vdin']
	Vdpc = param['Vdpc']
	Vdpcc = param['Vdpcc']
	Vdpcn = param['Vdpcn']
	Vmb = param['Vmb']
	Vmc = param['Vmc']
	Vmp = param['Vmp']

	# Maximum rate of synthesis/transcription : 

	Vstot = param['Vstot']
	Vsb = param['Vsb']
	Vsc = param['Vsc']
	Vsp = param['Vsp']

#--------------------------
# Kinetic equations : 
#--------------------------

	# mRNAs of per, Cry and Bmal : 

	dMp = Vsp * Bn**n/(Kap**n+Bn**n) - Vmp * Mp/(Kmp+Mp) - kdmp*Mp
	dMc = Vsc * Bn**n/(Kac**n+Bn**n) - Vmc * Mc/(Kmc + Mc) - kdmc*Mc
	dMb = Vsb * Kib**m/(Kib**m+Bn**m) - Vmb * Mb/(Kmb + Mb) - kdmb*Mb

	#Phosphorylated and non-phosphorylated proteins PER and CRY in the cytosol : 

	dPc = ksp * Mp - V1p*Pc/(Kp+Pc) + V2p * pcp/(Kdp + Pcp) + k4 * PCc - k3 * Pc * Cc - kdn * Pc
	dCc = ksc * Mc - V1c * Cc / (Kp +Cc) + V2c * Ccp/(Kdp + Ccp) + k4 * PCc - k3 * Pc * Cc - kdnc * Cc
	dPcp = V1p * Pc/(Kp + Pc) - V2p * Pcp/(Kdp + Pcp) - Vdpc * Pcp/(Kp+Pcp) - kdn * Pcp
	dCcp = V1c * Cc/(Kp+Cc) - V2c * Ccp/(Kdp + Ccp) - vdcc * Ccp/(Kd + Ccp) - kdnCcp

	# Phosphorylated and non-phosphorylated PER-CRY complex in cytosom and nucleus : 

	dPCc = -V1pc * PCc/(Kp+PCc) + V2pc * PCcp/(Kdp + PCcp) - k4 * PCc + k3 * Pc * Cc + k2 * Pcn - k1 * PCc - kdn * PCc 
	dPCn = -V3pc * Pcn/(Kp+Pcn) + V4pc * PCnp/(Kdp+PCnp) - k2*Pcn + k1*PCc - k7 * Bn * Pcn + k8 * In - kdn * Pcn
	dPCcp = V1pc * PCc/(Kp+PCc) - V2pc * PCcp/(Kdp + PCcp) - Vdpcc * PCcp/(Kd + PCcp) - kdn * PCcp
	dPCnp = V3pc * Pcn/(Kp+Pcn) - V4pc * PCnp/(Kdp + PCnp) - Vdpcn * PCnp/(Kd + PCnp) - kdn * PCnp

	#  Phosphorylated  and  non-phosphorylated  protein BMAL1 in the cytosol and nucleus
	dBc = kib * Mb - V1b * Bc/(Kp+Bc) + V2b * Bcp/(Kdp + Bcp) - k5*Bc + k6*Bc - kdn*Bc
	dBcp = V1b * Bc/(Kp + Bc) - V2b * Bcp/(Kdp + Bcp) - Vdbc * Bcp/(Kd + Bcp) - kdn*Bcp
	dBn = -V3b * Bn/(Kp+Bn) - V4b * Bnp/(Kdp+Bnp) + k5*Bc - k6 * Bn - k7 * Bn * Pcn + k8 * In - kdn*Bn
	dBnp = V3b*Bn/(Kp+Bn) - V4b * Bnp/(Kdp + Bnp) - Vdbn * Bnp/(Kd + Bnp) - kdn * Bnp

	#Inactive complex between PER–CRY and CLOCK–BMAL1 in nucleus :
	dIn = -k8 * In + k7 * Bn * Pcn -Vdin * In/(Kd + In) - kdn*In


	
	return [dMp, dMc, dMb, pPc, dCc, dPcp, dCcp, dPCc, dPCn, dPCcp, dPCnp, dBc, dBcp, dBn, dBnp, dIn]

