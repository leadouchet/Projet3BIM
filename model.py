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
			if (len(c)== 1) :
				parameters[data.index[i]] = parameters[c[0]]
			else : 
				parameters[data.index[i]] = float(c[0]) * parameters[c[1]]
	return parameters


def set_initials(x) :
	initial_condition = {}
	molecules = ['Mp', 'Mc', 'Mb', 'Pc', 'Cc', 'Pcp', 'Ccp', 'PCc', 'Pcn', 'PCcp', 'PCnp', 'Bc', 'Bcp', 'Bn', 'Bnp', 'In']
	for i in molecules : 
		initial_condition[i] = x
	return initial_condition


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

#-----------------------
# Initial conditions  : 
#-----------------------

	# mRNAs of per, Cry and Bmal : 
	Mp = initial_conditions['Mp']
	Mc = initial_conditions ['Mc']
	Mb = initial_conditions['Mb']

	# Phosporylated and non-phosphorylated proteins PER
	# and Cry in the cytosol : 

	Pc = initial_conditions['Pc']
	Bn = initial_conditions['Bn']
	Cc = initial_conditions['Cc']
	Pcp = initial_conditions['Pcp']
	Ccp = initial_conditions['Ccp']

	# Phosporylated and non-phosphorylated PER- Cry complexe
	# in the cytosol and nucleus : 

	PCc = initial_conditions['PCc']
	Pcn = initial_conditions['Pcn']
	PCcp = initial_conditions['PCcp']
	PCnp = initial_conditions['PCnp']

	# Phosphorylated and non-phosphorylated protein BMAL1 in
	# the cytosol and nucleus : 

	Bc = initial_conditions['Bc']
	Bcp = initial_conditions['Bcp']
	Bn = initial_conditions['Bn']
	Bnp = initial_conditions['Bnp']

	# Inactive complex between PER-CRY and CLOCK-BMAL1 in 
	# nucleus : 

	In = initial_conditions['In']

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

	KAP = param['KAP']
	KAC = param['KAC']
	KIB = param['KIB']

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
	KmB = param['KmB']
	KmC = param['KmC']
	KmP = param['KmP']

	#Rate constant for synthesis : 

	kstot = param['kstot']
	ksB = param['ksB']
	ksC = param['ksC']
	ksP = param['ksP']

	# Degree of cooperativity : 

	n = param['n']
	m = param['m']

	#Phosphorylation rate : 

	Vphos = param['Vphos']
	#Maximum Rate : 

	V1B = param['V1B']
	V1C = param['V1C']
	V1P = param['V1P']
	V1PC = param['V1PC']
	V2B = param['V2B']
	V2C = param['V2C']
	V2P = param['V2P']
	V2PC = param['V2PC']
	V3B = param['V3B']
	V3PC = param['V3PC']
	V4B = param['V4B']
	V4PC = param['V4PC']

	#Maximum rate of degradation

	vndBC = param['vndBC']
	vndBN = param['vndBN']
	vndCC = param['vndCC']
	vndIN = param['vndIN']
	vndPC = param['vndPC']
	vndPCC = param['vndPCC']
	vndPCN = param['vndPCN']
	vnmB = param['vnmB']
	vnmC = param['vnmC']
	vnmP = param['vnmP']

	# Maximum rate of synthesis/transcription : 

	vnsTot = param['vnsTot']
	vnsB = param['vnsB']
	vnsC = param['vnsC']
	vnsP = param['vnsP']

#--------------------------
# Kinetic equations : 
#--------------------------

	# mRNAs of per, Cry and Bmal : 

	dMp = vnsP * Bn**n/(KAP**n+Bn**n) - vnmP * Mp/(KmP+Mp) - kdmp*Mp
	dMc = vnsC * Bn**n/(KAC**n+Bn**n) - vnmC * Mc/(KmC + Mc) - kdmc*Mc
	dMb = vnsB * KIB**m/(KIB**m+Bn**m) - vnmB * Mb/(KmB + Mb) - kdmb*Mb

	#Phosphorylated and non-phosphorylated proteins PER and CRY in the cytosol : 

	dPc = ksP * Mp - V1P*Pc/(Kp+Pc) + V2P * Pcp/(Kdp + Pcp) + k4 * PCc - k3 * Pc * Cc - kdn * Pc
	dCc = ksC * Mc - V1C * Cc / (Kp +Cc) + V2C * Ccp/(Kdp + Ccp) + k4 * PCc - k3 * Pc * Cc - kdnc * Cc
	dPcp = V1P * Pc/(Kp + Pc) - V2P * Pcp/(Kdp + Pcp) - vndPC * Pcp/(Kp+Pcp) - kdn * Pcp
	dCcp = V1C * Cc/(Kp+Cc) - V2C * Ccp/(Kdp + Ccp) - vndCC * Ccp/(Kd + Ccp) - kdn * Ccp

	# Phosphorylated and non-phosphorylated PER-CRY complex in cytosom and nucleus : 

	dPCc = -V1PC * PCc/(Kp+PCc) + V2PC * PCcp/(Kdp + PCcp) - k4 * PCc + k3 * Pc * Cc + k2 * Pcn - k1 * PCc - kdn * PCc 
	dPCn = -V3PC * Pcn/(Kp+Pcn) + V4PC * PCnp/(Kdp+PCnp) - k2*Pcn + k1*PCc - k7 * Bn * Pcn + k8 * In - kdn * Pcn
	dPCcp = V1PC * PCc/(Kp+PCc) - V2PC * PCcp/(Kdp + PCcp) - vndPCC * PCcp/(Kd + PCcp) - kdn * PCcp
	dPCnp = V3PC * Pcn/(Kp+Pcn) - V4PC * PCnp/(Kdp + PCnp) - vndPCN * PCnp/(Kd + PCnp) - kdn * PCnp

	#  Phosphorylated  and  non-phosphorylated  protein BMAL1 in the cytosol and nucleus
	dBc = KIB * Mb - V1B * Bc/(Kp+Bc) + V2B * Bcp/(Kdp + Bcp) - k5*Bc + k6*Bc - kdn*Bc
	dBcp = V1B * Bc/(Kp + Bc) - V2B * Bcp/(Kdp + Bcp) - vndBC * Bcp/(Kd + Bcp) - kdn*Bcp
	dBn = -V3B * Bn/(Kp+Bn) - V4B * Bnp/(Kdp+Bnp) + k5*Bc - k6 * Bn - k7 * Bn * Pcn + k8 * In - kdn*Bn
	dBnp = V3B*Bn/(Kp+Bn) - V4B * Bnp/(Kdp + Bnp) - vndBN * Bnp/(Kd + Bnp) - kdn * Bnp

	#Inactive complex between PER–CRY and CLOCK–BMAL1 in nucleus :
	dIn = -k8 * In + k7 * Bn * Pcn -vndIN * In/(Kd + In) - kdn*In

	return [dMp, dMc, dMb, dPc, dCc, dPcp, dCcp, dPCc, dPCn, dPCcp, dPCnp, dBc, dBcp, dBn, dBnp, dIn]


param = readparam('param.csv', 1)
print param
init = set_initials(1)
print CircadianRythme(10, init, param)








