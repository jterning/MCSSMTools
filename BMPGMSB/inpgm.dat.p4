# INPUT FILE FOR NMSSMTools
# BASED ON SUSY LES HOUCHES ACCORD II

BLOCK MODSEL
	3	1		# NMSSM particle content
	1	2		# IMOD (0=general NMSSM, 1=mSUGRA, 2=GMSB)
	10	0		# ISCAN (0=no scan, 1=grid scan, 2=random scan)

BLOCK SMINPUTS
	1	127.92D0	# ALPHA_EM^-1(MZ)
	2	1.16639D-5	# GF
	3	.1172D0		# ALPHA_S(MZ)
	4	91.187D0	# MZ
	5	4.214D0		# MB(MB) (running mass)
	6	171.4D0		# MTOP (pole mass)
	7	1.777D0		# MTAU

BLOCK MINPAR
	1	5.D4		# MSUSYEFF = m^2/MMESS
	2	7.5D6		# MMESS
	3	1.9D0		# TANB AT MSUSY
	4	1.D0		# SIGMU
	5	2.D0		# N5 = number of messenger 5-plets

BLOCK EXTPAR
	61	.6D0		# LAMBDA AT MSUSY
	63	0.D0		# ALAMBDA at MMESS
	66      0.D0		# XIF in GeV^2 at MMESS
	67      0.D0		# XIS in GeV^2 at MMESS (If MS^2 is not an input)
	68	0.D0		# MU' at MMESS
	69	0.D0		# MS'^2 at MMESS
	71	0.D0		# DH, Shift in mH at MMESS, |DH|<1
