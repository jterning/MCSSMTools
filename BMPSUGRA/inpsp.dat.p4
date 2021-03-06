# Input file for NMSSMTools
# Based on SUSY LES HOUCHES ACCORD II

BLOCK MODSEL
	3	1		# NMSSM particle content
	1	1		# IMOD (0=general NMSSM, 1=mSUGRA, 2=GMSB)
	10	0		# ISCAN (0=no scan, 1=grid scan, 2=random scan)
	9	1		# Call micrOmegas (default 0=no, 1=relic density only,
#				  2=dir. det. rate, 3=indir. det. rate, 4=both det. rates)

BLOCK SMINPUTS
	1	127.92D0	# ALPHA_EM^-1(MZ)
	2	1.16639D-5	# GF
	3	.1172D0		# ALPHA_S(MZ)
	4	91.187D0	# MZ
	5	4.214D0		# MB(MB) (running mass)
	6	171.4D0		# MTOP (pole mass)
	7	1.777D0		# MTAU

BLOCK MINPAR
	1	780.D0		# M0
	2	775.D0		# M12
	3	3.05D0		# TANB at MSUSY
	4	-1.D0		# SIGMU
	5	-2250.D0	# A0  

BLOCK EXTPAR
	21	7.25D5		# MHD^2 at MGUT
	22	4.65D6		# MHU^2 at MGUT
	61	.49D0		# LAMBDA at MSUSY
	64	-965.D0		# AKAPPA at MGUT

