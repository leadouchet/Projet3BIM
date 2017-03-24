def CircadianRythme(t, initial_conditions, param) :

"""This function allows to simulate the cycle of proteins 
during the circadian cycle. The inital_condition is a 
dictionnary containing as keys the names of the molecule, 
param is also a dictonnary containing as keys the name of the 
parameter. dt and Run_time are numbers respectively describing 
the integration step and the  Running time """

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