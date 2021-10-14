import sys,os
sys.path.append('/home/jlw2245/Pynterface')
import pynterface2 as pynt
from pyscf import gto

atoms     = '''
C 0.00000000000000 0.00000000000000 -0.64944465799361
O 0.00000000000000 0.00000000000000 0.48844465799361
'''
basis     = {'C': gto.parse('''
C    S
   4563.240                  0.00196665
    682.0240                 0.0152306
    154.9730                 0.0761269
     44.45530                0.2608010
     13.02900                0.6164620
      1.827730               0.2210060
C    SP
     20.96420                0.114660               0.0402487
      4.803310               0.919999               0.237594
      1.459330              -0.00303068             0.815854
C    SP
      0.4834560              1.000000               1.000000
C    SP
      0.1455850              1.000000               1.000000
C    SP
      0.0438000              1.0000000              1.0000000
C    D
      2.50400                1.000000
C    D
      0.62600                1.000000
C    D
      0.15650                1.000000
C    F
      0.8000000              1.0000000
      '''),
      'O': gto.parse('''
O    S
   8588.500                  0.00189515
   1297.230                  0.0143859
    299.2960                 0.0707320
     87.37710                0.2400010
     25.67890                0.5947970
      3.740040               0.2808020
O    SP
     42.11750                0.113889               0.0365114
      9.628370               0.920811               0.237153
      2.853320              -0.00327447             0.819702
O    SP
      0.905661               1.000000               1.000000
O    SP
      0.255611               1.000000               1.000000
O    SP
      0.0845000              1.0000000              1.0000000
O    D
      5.16000                1.000000
O    D
      1.29200                1.000000
O    D
      0.32250                1.000000
O    F
      1.4000000              1.0000000      
      ''')}   #6-311++G(3df,3pd)
chrg      = 0          # charge
spn       = 0          # number of unpaired electrons (e.g. 1 in ferrocenium NOT the multiplicity which is 2 for ferrocenium )
psp       = None       # whether or not to use pseudopotential
sym       = 'c1'       # use 'c1' if possible ( no symmetry ); insert here 1 to atomatically detect the symmetry group (Note:  CASSCF wavefunction much more compact, but energy will be higher...)
atomSym   = False      # custom CAS symmetrizing procedure (rarely use)
psit      = 'rhf'     # |psi_T> , the trial wavefunction (choices: rhf, uhf, rcas, ucas, uks); ph-AFQMC cannot use UCAS trial (so don't use)

# if the trial wavefunction is CASSCF
na_act    = 9         # number of active electrons with spin-up
nb_act    = 9          # same for spin-down BR
ncas_act  = 14         # number of active orbitals
tot_wt    = 0.995      # tot_wt = sum_i |c_i|^2.  When i goes over all CASSCF configurations, tot_wt = 1; for truncated CASSCF, < 1.
smallCASsteps = 0      # use only when want to converge on a solution "close" to the initial guess
casscf_stepsize = None   # use only if you want to change the stepsize (overwrites smallCASsteps option)
CASMaxmem = 64000      # MB (looks like f48-19 has 132 GB, but double check)
FCIMaxmem = 64000      # MB
CASReadchk = 0         # set to true if to read from CAS_chk_name the CAS typically for AVAS
CAS_chk_name = 'oldchk'        #put in name typically called oldchk if you want to read from old chkfile name
casMOLDENgenerate = 1 			#whether to produce a file that can be opened by Molden for the CAS orbitals; will not work for QZ and up
cas_cycles = 1000      #number of CASSCF cycles
auto_mp_noons_cas = False               # set to true to run an mp2 calculation on the
                                        #ROHF wfn and use the fractional NOONs from that 
uno_cas = ''           #set to either uhf or uks to use UNOs from these
                       #to choose active space
cas_list = [6,7,8,9,10,11,12,13]        #this an example only works if next option is 1; pick orbitals for CAS space, 1-based indices
pick_cas_list = 0  			#put 1 if you want to pick the CAS list
ao_labels = ['5 Fe 3d', '5 Fe 4d', 'C 2pz'] 	#labels needed by AVAS to pick orbitals
AVAS_threshold = 0 			# AVAS flag for the threshold; if positive then AVAS is activated
avas_basis_set = '' 		        # set to ano if desired (otherwise will use minao)
CASCI = 0			        #whether to use CASCI instead of CASSCF solver if rcas is selected
SHCI = 0                                # whether to use SHCISCF as CASSCF solver
shci_np = 1                             # number of processes in shci
shci_ndet = 1e6                         # number of determinants to print out from DICE
level_shift_cas = None #whether to use/change default value of level_shift for the CASSCF solver
NEVPT2_flag = 0 #whether to run NEVPT2 on top of the CASSCF calculation
PIOS_flag = 0 #whether to use PIOS (cannot be combined with AVAS yet)
PiAtomsList = [6,7,8,9,10,11,12,13] #these atoms compose the Pi system
number_pi_occupied = None #use entire pi space or use integer pi occupied orbitals (set both similarly)
number_pi_virtual = None #use entire pi space or use integer pi virtual orbitals (set both similarly)

# if the trial wavefunction is rks
fun       = None       # dft xc functional, e.g. "b3lyp"
read_in_qchem = False   # whether or not to read in MOs from QChem file 53.0
qchem_scratch_dir = './' #directory in which 53.0 is stored
#for ROHF, UHF and UKS:
init_with_rks = 0 #use ROKS to initialization otherwise use ROHF - doesn't work yet.
rhfMOLDENgenerate = 1 			#whether to produce a file that can be opened by Molden for the RHF orbitals ; will not work for QZ and up
x2c       = 0          # set 1 if doing scalar relativistic calculation, 0 o/w
conv      = 1e-10      # convergence of scf calculations
lvl_shft_rhf  = 0.0        # level shifting in scf calculations (may need to converge TMs)
lvl_shft_trial  = 0.0        # level shifting in scf calculations (may need to converge TMs)
lindep    = 0          # linear dependence threshold (need only when have very large / "redundant" basis set)
HFMaxmem  = 64000      # memory of scf calculation in MB
noise     = 0.0        # whether to add some noise to the scf calculations, to help going away from local minima - only for UHF and UKS
printh    = True       # whether to print the hamiltonian (needed to print QMC inputs)
dec       = 'cholesky' # decomposition of electron repulsion integrals (V_ijkl), choices are 'cholesky' or 'df'
initguess_rhf = ''     # atom, 1e, or hcore for rhf beginning #default is minao but don't need to write that
initguess_trial = ''     # atom, 1e, or hcore for the trial wavefunction #default is minao but don't need to write that
scf_cycles_rhf = 500      # how many cycles allowed for SCF in the RHF initial guess #recommended 500
scf_cycles_trial = 500      # how many cycles allowed for SCF in the trial #recommend 500 for UHF, 5000 for DFT/KS
nstability_iter_rhf = 2 		#how many times to do stability checks for the initial guess ; 2 is recommended (will also check original stability)
nstability_iter_trial = 2 		#how many times to do stability checks for the UHF trial ; 2 is recommended (will also check original stability)
newton_rhf = 0 				#whether to use Newton solver for the RHF initial guess (1 for yes)
newton_trial = 0 			#whether to use Newton solver for the UHF trial (1 for yes)
RHFReadchk = 0                  #change to 1 if we will read in the chk file for the ROHF part
RHFchk_name = '' 		#put in name of chk file for ROHF to start from - can project from small basis to large basis e.g. 'TZ.chk' to go to QZ from TZ
#rename the chkfile from rhf.chk so the code does not override it
UHForUKSReadchk = 0                  #change to 1 if we will read in the chk file for the ROHF part
UHForUKSReadchk_name = '' 		#put in name of chk file for UHF or UKS to start from - can project from small basis to large basis e.g. 'TZ.chk' to go to QZ from TZ
#rename the chkfile from uhf.chk or uks.chk so the code does not override it
uhf_noons = False
#ddcosmo dielectric settings
ddcosmoscrf = 0 #whether to turn dielectric on
dielectric = 78.3553 # dielectric value to use e.g. this one is for water

# if the decomposition is cholesky
storage   = 2          # 1 for small systems; 2 for large (less efficient, but necessary when 1 crashes, guess ~250bfs)
chtol     = 1e-4       # threshold of cholesky decomposition (1e-4 is fairly aggressive, 1e-5 is safest)
chmax     = None       # ratio between number of cholesky vectors and basis functions; default is 10

# if the decomposition is df, choose auxiliary basis
auxbasis  = 'weigend+etb'  #e.g. ccpvtzri

#----------------------- don't touch these ------------------------------#
thc       = False # whether to do the thc
eigtol    = 1e-4  # eigenvalue threshold for thc decomposition
props     = False # whether to compute properties (1rdm, dipole moment)
#------------------------------------------------------------------------#

# number of orbitals to freeze, highest orbital to keep (frozen core, downfolding procedures); fcbasis!=None if the basis comes from external input
#not sure these work yet - BR
frozen    = None
downfold  = None
fcbasis   = None

pynt.gafqmc_setup(
 atoms     = atoms,
 chrg      = chrg,
 spn       = spn,
 basis     = basis,
 sym       = sym,
 psp       = psp,
 frags     = None,
 psit      = psit,
 conv      = conv,
 lindep    = lindep,
 HFMaxmem  = HFMaxmem,
 noise     = noise,
 lvl_shft_rhf   = lvl_shft_rhf,
 lvl_shft_trial = lvl_shft_trial,
 fun            = fun,
 init_with_rks  = init_with_rks,
 na_act    = na_act,
 nb_act    = nb_act,
 ncas_act  = ncas_act,
 tot_wt    = tot_wt,
 smallCASsteps  = smallCASsteps,
 CASMaxmem = CASMaxmem,
 FCIMaxmem = FCIMaxmem,
 dec       = dec,
 chtol     = chtol,
 chmax     = chmax,
 storage   = storage,
 auxbasis  = auxbasis,
 thc       = thc,
 eigtol    = eigtol,
 props     = props,
 printh    = printh,
 frozen    = frozen,
 downfold  = downfold,
 fcbasis   = fcbasis,
 x2c       = x2c,
 atomSym   = atomSym,
 CASReadchk     = CASReadchk,
 CAS_chk_name   = CAS_chk_name,
 AVAS_threshold = AVAS_threshold,
 ao_labels = ao_labels,
 casMOLDENgenerate = casMOLDENgenerate,
 scf_cycles_trial  = scf_cycles_trial,
 avas_basis_set    = avas_basis_set,
 initguess_rhf     = initguess_rhf,
 initguess_trial   = initguess_trial,
 newton_rhf        = newton_rhf,
 newton_trial      = newton_trial,
 pick_cas_list     = pick_cas_list,
 cas_cycles        = cas_cycles,
 nstability_iter_rhf   = nstability_iter_rhf,
 nstability_iter_trial = nstability_iter_trial,
 cas_list          = cas_list,
 scf_cycles_rhf    = scf_cycles_rhf,
 rhfMOLDENgenerate = rhfMOLDENgenerate,
 CASCI             = CASCI,
 level_shift_cas   = level_shift_cas,
 RHFReadchk        = RHFReadchk,
 RHFchk_name       = RHFchk_name,
 UHForUKSReadchk   = UHForUKSReadchk,
 UHForUKSReadchk_name = UHForUKSReadchk_name,
 ddcosmoscrf       = ddcosmoscrf,
 dielectric        = dielectric,
 casscf_stepsize   = casscf_stepsize,
 NEVPT2_flag       = NEVPT2_flag,
 PIOS_flag         = PIOS_flag,
 PiAtomsList       = PiAtomsList,
 number_pi_occupied = number_pi_occupied,
 number_pi_virtual  = number_pi_virtual,
 SHCI = SHCI, 
 shci_np = shci_np, 
 shci_ndet = shci_ndet,
 auto_mp_noons_cas = auto_mp_noons_cas,
 uno_cas = uno_cas,
 uhf_noons = uhf_noons,
 read_in_qchem = read_in_qchem,
 qchem_scratch_dir = qchem_scratch_dir)
