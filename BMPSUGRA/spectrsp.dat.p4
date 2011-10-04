# NMSSMTools OUTPUT IN SLHA FORMAT
# Info about spectrum calculator
BLOCK SPINFO   # Program information
     1   NMSSMTools # Spectrum calculator
     2   3.0.1      # Version number
     8   0          # Higgs mass precision
     3   # Muon magn. mom. more than 2 sigma away
# Input parameters
BLOCK MODSEL
    3     1         # NMSSM PARTICLE CONTENT
BLOCK SMINPUTS
     1     1.27920000E+02   # ALPHA_EM^-1(MZ)
     2     1.16639000E-05   # GF
     3     1.17200000E-01   # ALPHA_S(MZ)
     4     9.11870000E+01   # MZ
     5     4.21400000E+00   # MB(MB)
     6     1.71400000E+02   # MTOP (POLE MASS)
     7     1.77700000E+00   # MTAU
# SMINPUTS Beyond SLHA:
# MW:     0.80420000E+02
# MS:     0.19000000E+00
# MC:     0.14000000E+01
# VUS:     0.22000000E+00
# VCB:     0.40000000E-01
# VUB:     0.40000000E-02
BLOCK MINPAR
     4    -1.00000000E+00   # SIGMU
     3     3.05000000E+00   # TANBETA(MZ)
     1     7.80000000E+02   # M0(MGUT)
     2     7.75000000E+02   # M12(MGUT)
     5    -2.25000000E+03   # A0(MGUT)
BLOCK EXTPAR
    21     7.25000000E+05   # MHD^2 AT THE GUT SCALE
    22     4.65000000E+06   # MHU^2 AT THE GUT SCALE
    61     4.90000000E-01   # LAMBDA AT THE SUSY SCALE
    64    -9.65000000E+02   # AKAPPA AT THE GUT SCALE
    66     0.00000000E+00   # XIF AT THE GUT SCALE
    67     0.00000000E+00   # XIS AT THE GUT SCALE
    68     0.00000000E+00   # MU' AT THE GUT SCALE 
    69     0.00000000E+00   # MS'^2 AT THE GUT SCALE 
    69     0.00000000E+00   # M3H^2 AT THE GUT SCALE 
# 
BLOCK MASS   # Mass spectrum 
#  PDG Code     mass             particle 
        25     4.06285454E+01   # lightest neutral scalar
        35     1.20772449E+02   # second neutral scalar
        45     6.83318210E+02   # third neutral scalar
        36     1.14654574E+02   # lightest pseudoscalar
        46     6.82975010E+02   # second pseudoscalar
        37     6.78441838E+02   # charged Higgs
   1000001     1.74230726E+03   #  ~d_L
   2000001     1.67388885E+03   #  ~d_R
   1000002     1.74089743E+03   #  ~u_L
   2000002     1.73674821E+03   #  ~u_R
   1000003     1.74230726E+03   #  ~s_L
   2000003     1.67388885E+03   #  ~s_R
   1000004     1.74089743E+03   #  ~c_L
   2000004     1.73674821E+03   #  ~c_R
   1000005     1.26381972E+03   #  ~b_1
   2000005     1.62512398E+03   #  ~b_2
   1000006     4.13282989E+02   #  ~t_1
   2000006     1.28743520E+03   #  ~t_2
   1000011     9.80803188E+02   #  ~e_L
   2000011     7.03205083E+02   #  ~e_R
   1000012     9.78291958E+02   #  ~nue_L
   1000013     9.80803188E+02   #  ~mu_L
   2000013     7.03205083E+02   #  ~mu_R
   1000014     9.78291958E+02   #  ~numu_L
   1000015     6.99984270E+02   #  ~tau_1
   2000015     9.79679540E+02   #  ~tau_2
   1000016     9.77151419E+02   #  ~nutau_L
   1000021     1.75889131E+03   #  ~g
   1000022    -6.07124336E+01   # neutralino(1)
   1000023     2.25840252E+02   # neutralino(2)
   1000025    -2.27491310E+02   # neutralino(3)
   1000035     3.34836971E+02   # neutralino(4)
   1000045     6.40537919E+02   # neutralino(5)
   1000024     2.13215865E+02   # chargino(1)
   1000037    -6.40512215E+02   # chargino(2)
# 
# Low energy observables
BLOCK LOWEN
# Exp. 2 Sigma: 3.03E-04 < BR(b -> s gamma) < 4.01E-04:
         1     3.89763079E-04   # BR(b -> s gamma)
        11     4.35529141E-04   # (BR(b -> s gamma)+Theor.Err.)
        12     3.24668396E-04   # (BR(b -> s gamma)-Theor.Err.)
# Exp. 2 Sigma: 4.99E-01 < Delta M_d < 5.15E-01:
         2     6.29352938E-01   # Delta M_d in ps^-1
        21     1.10721245E+00   # Delta M_d +Theor.Err.
        22     1.67558211E-01   # Delta M_d -Theor.Err.
# Exp. 2 Sigma: 1.753E+01 < Delta Ms < 1.801E+01:
         3     2.18042691E+01   # Delta M_s in ps^-1
        31     2.89865426E+01   # Delta M_s +Theor.Err.
        32     1.48560157E+01   # Delta M_s -Theor.Err.
# Exp. 95% C.L.: BR(Bs->mu+mu-) < 1.1E-08:
         4     3.54043574E-09   # BR(Bs -> mu+mu-)
        41     6.01336345E-09   # BR(Bs -> mu+mu-)+Theor.Err.
        42     1.71864510E-09   # BR(Bs -> mu+mu-)-Theor.Err.
# Exp. 2 Sigma: 8.90E-05 < BR(B+ > tau+ + nu_tau) < 2.45E-04:
         5     1.31652145E-04   # BR(B+ -> tau+ + nu_tau)
        51     2.63355788E-04   # BR(B+ -> tau+ + nu_tau) + Theor.Err.
        52     5.68179067E-05   # BR(B+ -> tau+ + nu_tau) - Theor.Err.
# BSM contr. to the muon anomalous magn. moment:
         6    -8.44883547E-11   # Del_a_mu
        61     1.97371802E-10   # Del_a_mu + Theor.Err.
        62    -3.66348511E-10   # Del_a_mu - Theor.Err.
# Minimal Exp.-SM (2 sigma):  8.77306222E-10
# Maximal Exp.-SM (2 sigma):  4.61144414E-09
# 
# Omega h^2 (allowed: 0.094 < Omega h^2 < 0.136):
        10     1.04541764E-01   # Omega h^2
# 
BLOCK HMIX Q=  8.14768987E+02 # (STOP/SBOTTOM MASSES)
     1    -2.06172995E+02   # MUEFF
     2     3.04128389E+00   # TAN(BETA)
     3     2.42973521E+02   # V(Q)
     4     4.65389357E+05   # MA^2
     5     1.86786342E+04   # MP^2
# 
# 3*3 Higgs mixing
BLOCK NMHMIX
  1  1    -9.95669073E-02   # S_(1,1)
  1  2     1.82323891E-02   # S_(1,2)
  1  3     9.94863815E-01   # S_(1,3)
  2  1     3.16460670E-01   # S_(2,1)
  2  2     9.48498007E-01   # S_(2,2)
  2  3     1.42890165E-02   # S_(2,3)
  3  1     9.43365823E-01   # S_(3,1)
  3  2    -3.16257983E-01   # S_(3,2)
  3  3     1.00208847E-01   # S_(3,3)
# 
# 3*3 Pseudoscalar Higgs mixing
BLOCK NMAMIX
  1  1     1.05748082E-01   # P_(1,1)
  1  2     3.46715023E-02   # P_(1,2)
  1  3     9.93788323E-01   # P_(1,3)
  2  1     9.44327022E-01   # P_(2,1)
  2  2     3.09615417E-01   # P_(2,2)
  2  3    -1.11286881E-01   # P_(2,3)
# 
# 3rd generation sfermion mixing
BLOCK STOPMIX  # Stop mixing matrix
  1  1     1.71797137E-01   # Rst_(1,1)
  1  2     9.85132348E-01   # Rst_(1,2)
  2  1    -9.85132348E-01   # Rst_(2,1)
  2  2     1.71797137E-01   # Rst_(2,2)
BLOCK SBOTMIX  # Sbottom mixing matrix
  1  1     9.99972305E-01   # Rsb_(1,1)
  1  2     7.44237630E-03   # Rsb_(1,2)
  2  1    -7.44237630E-03   # Rsb_(2,1)
  2  2     9.99972305E-01   # Rsb_(2,2)
BLOCK STAUMIX  # Stau mixing matrix
  1  1     7.17858930E-03   # Rsl_(1,1)
  1  2     9.99974234E-01   # Rsl_(1,2)
  2  1    -9.99974234E-01   # Rsl_(2,1)
  2  2     7.17858930E-03   # Rsl_(2,2)
# 
# Gaugino-Higgsino mixing
BLOCK NMNMIX  # 5*5 Neutralino Mixing Matrix
  1  1    -3.46165710E-02   # N_(1,1)
  1  2     3.64947853E-02   # N_(1,2)
  1  3     3.92864993E-02   # N_(1,3)
  1  4     3.43379942E-01   # N_(1,4)
  1  5     9.37025405E-01   # N_(1,5)
  2  1    -1.58780486E-01   # N_(2,1)
  2  2     7.57626138E-02   # N_(2,2)
  2  3     7.03960671E-01   # N_(2,3)
  2  4     6.32820758E-01   # N_(2,4)
  2  5    -2.70233316E-01   # N_(2,5)
  3  1    -6.57932930E-02   # N_(3,1)
  3  2     7.97970426E-02   # N_(3,2)
  3  3    -7.05908709E-01   # N_(3,3)
  3  4     6.65356390E-01   # N_(3,4)
  3  5    -2.19766796E-01   # N_(3,5)
  4  1     9.84385754E-01   # N_(4,1)
  4  2     3.46655176E-02   # N_(4,2)
  4  3     6.77363404E-02   # N_(4,3)
  4  4     1.56707343E-01   # N_(4,4)
  4  5    -2.52504715E-02   # N_(4,5)
  5  1    -1.56965134E-02   # N_(5,1)
  5  2     9.92652452E-01   # N_(5,2)
  5  3    -7.92161292E-04   # N_(5,3)
  5  4    -1.19882387E-01   # N_(5,4)
  5  5     4.72381081E-03   # N_(5,5)
# 
BLOCK UMIX  # Chargino U Mixing Matrix
  1  1     9.65407873E-04   # U_(1,1)
  1  2     9.99999534E-01   # U_(1,2)
  2  1    -9.99999534E-01   # U_(2,1)
  2  2     9.65407873E-04   # U_(2,2)
# 
BLOCK VMIX  # Chargino V Mixing Matrix
  1  1     1.69042110E-01   # V_(1,1)
  1  2    -9.85608830E-01   # V_(1,2)
  2  1     9.85608830E-01   # V_(2,1)
  2  2     1.69042110E-01   # V_(2,2)
# 
# Higgs reduced couplings
# (as compared to a SM Higgs with same mass)
BLOCK REDCOUP
# H1
  1  1     1.91873524E-02   # U-type fermions
  1  2    -3.19584957E-01   # D-type fermions
  1  3    -1.36951819E-02   # W,Z bosons
  1  4     1.66778509E-01   # Gluons
  1  5     7.22519558E-02   # Photons
# H2
  2  1     9.98177772E-01   # U-type fermions
  2  2     1.01575988E+00   # D-type fermions
  2  3     9.99884359E-01   # W,Z bosons
  2  4     9.73260860E-01   # Gluons
  2  5     9.95316105E-01   # Photons
# H3
  3  1    -3.32822722E-01   # U-type fermions
  3  2     3.02796917E+00   # D-type fermions
  3  3    -6.61142434E-03   # W,Z bosons
  3  4     3.18960186E-01   # Gluons
  3  5     1.36547467E+00   # Photons
# A1
  4  1     3.64875020E-02   # U-type fermions
  4  2     3.39424987E-01   # D-type fermions
  4  3     0.00000000E+00   # W,Z bosons
  4  4     3.24760135E-02   # Gluons
  4  5     2.76038480E-01   # Photons
# A2
  5  1     3.25832237E-01   # U-type fermions
  5  2     3.03105438E+00   # D-type fermions
  5  3     0.00000000E+00   # W,Z bosons
  5  4     3.32382381E-01   # Gluons
  5  5     2.87568154E-01   # Photons
# 
# GAUGE AND YUKAWA COUPLINGS AT THE SUSY SCALE
BLOCK GAUGE Q=  1.66075041E+03 # (SUSY SCALE)
         1     3.64129060E-01   # g1(Q,DR_bar)
         2     6.42230144E-01   # g2(Q,DR_bar)
         3     1.03597648E+00   # g3(Q,DR_bar)
BLOCK YU Q=  1.66075041E+03 # (SUSY SCALE)
  3  3     8.82706126E-01   # HTOP(Q,DR_bar)
BLOCK YD Q=  1.66075041E+03 # (SUSY SCALE)
  3  3     4.38586231E-02   # HBOT(Q,DR_bar)
BLOCK YE Q=  1.66075041E+03 # (SUSY SCALE)
  3  3     3.18976709E-02   # HTAU(Q,DR_bar)
# 
# SOFT TRILINEAR COUPLINGS AT THE SUSY SCALE
BLOCK AU Q=  1.66075041E+03 # (SUSY SCALE)
  3  3    -1.76577142E+03   # ATOP
BLOCK AD Q=  1.66075041E+03 # (SUSY SCALE)
  3  3    -3.73049948E+03   # ABOT
BLOCK AE Q=  1.66075041E+03 # (SUSY SCALE)
  2  2    -2.53124732E+03   # AMUON
  3  3    -2.52476298E+03   # ATAU
# 
# SOFT MASSES AT THE SUSY SCALE
BLOCK MSOFT Q=  1.66075041E+03 # (SUSY SCALE)
         1     3.33803265E+02   # M1
         2     6.13911940E+02   # M2
         3     1.67034315E+03   # M3
        21     3.75581052E+05   # M_HD^2
        22     8.11178852E+04   # M_HU^2
        31     9.79927069E+02   # M_eL
        32     9.79927069E+02   # M_muL
        33     9.78788436E+02   # M_tauL
        34     7.02149068E+02   # M_eR
        35     7.02149068E+02   # M_muR
        36     6.98938449E+02   # M_tauR
        41     1.67918623E+03   # M_q1L
        42     1.67918623E+03   # M_q2L
        43     1.25768032E+03   # M_q3L
        44     1.67451313E+03   # M_uR
        45     1.67451313E+03   # M_cR
        46     5.27835643E+02   # M_tR
        47     1.60904962E+03   # M_dR
        48     1.60904962E+03   # M_sR
        49     1.60663232E+03   # M_bR
# 
# NMSSM SPECIFIC PARAMETERS THE SUSY SCALE
BLOCK NMSSMRUN Q=  1.66075041E+03 # (SUSY SCALE)
     1     4.90000000E-01   # LAMBDA(Q,DR_bar)
     2     5.54362462E-02   # KAPPA(Q,DR_bar)
     3    -6.84510509E+02   # ALAMBDA
     4     1.51726925E+02   # AKAPPA
     5    -2.07635596E+02   # MUEFF
     6     0.00000000E+00   # XIF
     7     0.00000000E+00   # XIS
     8     0.00000000E+00   # MU'
     9     0.00000000E+00   # MS'^2
    10     6.07350035E+03   # MS^2
    12     0.00000000E+00   # M3H^2
# 
# GAUGE AND YUKAWA COUPLINGS AT THE GUT SCALE
BLOCK GAUGE Q=  1.74527935E+16 # (GUT SCALE)
         1     7.09781105E-01   # g1(MGUT,DR_bar), GUT normalization
         2     7.09781104E-01   # g2(MGUT,DR_bar)
         3     6.99955531E-01   # g3(MGUT,DR_bar)
BLOCK YU Q=  1.74527935E+16 # (GUT SCALE)
  3  3     6.46595828E-01   # HTOP(MGUT,DR_bar)
BLOCK YD Q=  1.74527935E+16 # (GUT SCALE)
  3  3     1.91122391E-02   # HBOT(MGUT,DR_bar)
BLOCK YE Q=  1.74527935E+16 # (GUT SCALE)
  3  3     2.35545716E-02   # HTAU(MGUT,DR_bar)
# 
# SOFT TRILINEAR COUPLINGS AT THE GUT SCALE
BLOCK AU Q=  1.74527935E+16 # (GUT SCALE)
  3  3    -2.25000740E+03   # ATOP
BLOCK AD Q=  1.74527935E+16 # (GUT SCALE)
  3  3    -2.25000377E+03   # ABOT
BLOCK AE Q=  1.74527935E+16 # (GUT SCALE)
  2  2    -2.25000314E+03   # AMUON
  3  3    -2.25000321E+03   # ATAU
# 
# SOFT MASSES SQUARED AT THE GUT SCALE
BLOCK MSOFT Q=  1.74527935E+16 # (GUT SCALE)
         1     7.75000197E+02   # M1
         2     7.75000042E+02   # M2
         3     7.74999941E+02   # M3
        21     7.24938711E+05   # M_HD^2
        22     4.64992415E+06   # M_HU^2
        31     6.08399757E+05   # M_eL^2
        32     6.08399757E+05   # M_muL^2
        33     6.08399751E+05   # M_tauL^2
        34     6.08400011E+05   # M_eR^2
        35     6.08400011E+05   # M_muR^2
        36     6.08400000E+05   # M_tauR^2
        41     6.08399751E+05   # M_q1L^2
        42     6.08399751E+05   # M_q2L^2
        43     6.08395203E+05   # M_q3L^2
        44     6.08399826E+05   # M_uR^2
        45     6.08399826E+05   # M_cR^2
        46     6.08390692E+05   # M_tR^2
        47     6.08399931E+05   # M_dR^2
        48     6.08399931E+05   # M_sR^2
        49     6.08399924E+05   # M_bR^2
# 
# NMSSM SPECIFIC PARAMETERS AT THE GUT SCALE
BLOCK NMSSMRUN Q=  1.74527935E+16 # (GUT SCALE)
     1     6.53512128E-01   # LAMBDA(MGUT,DR_bar)
     2     8.12028152E-02   # KAPPA(MGUT,DR_bar)
     3    -2.25002310E+03   # ALAMBDA
     4    -9.65042871E+02   # AKAPPA
     6     0.00000000E+00   # XIF
     7     0.00000000E+00   # XIS
     8     0.00000000E+00   # MU'
     9     0.00000000E+00   # MS'^2
    10     1.41789949E+06   # MS^2
    12     0.00000000E+00   # M3H^2
