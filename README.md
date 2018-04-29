############
# To install
#############

This module requires hookups inside the VASP program. These hookups have been added to
the current VASP 6 development version and will be released with it. In the meantime,
experienced licensed VASP users may obtain patch files detailing them from the authors.

Once the hookups are present, replace the dummy extpot.F file shipped with VASP with
the one present in this github and (re)compile VASP (see VASP manual).

###################
# TO RUN - EXAMPLE
###################

cd /run/path

cp -r $VASPEMB_PATH/examples/cl2_test/* ./

cp $VASPEMB_PATH/drivers/* ./

# open svasp.emb, search for "VASP_EXEC", change its value to: $INSTALL_PATH/vasp.5.3.3/vasp
# open optimize.py, search for "VASP_EXEC", change its value to: $INSTALL_PATH/vasp.5.3.3/vasp

# do reference calculation
cd ref
../svasp.emb -np 16 -walltime 1:00:00 -mem 60000mb ref
# this will submit the reference calculation to the SLURM system, wait until it's finished

# Split the DERINFO and POSCAR use split_DERINFO.py
# NOTE: for other systems, make sure you configured the cluster definition in split_DERINFO.py.
# NOTE: The cluster definition is configured by the following variables in split_DERINFO.py
# n_atom: total number of atoms, including cluster and environtment.
# cluster: the atom number list for cluster. ALERT: atom number start from 1, yes, even it's in python!
# na_cluster: a dictionary defining the number of atoms by elements within cluster, used to make POSCAR for cluster
#             key   - element labels;
#             value - numbers of this particular type of element

cd /run/path
./split_DERINFO.py # This should generate POSCAR.cluster, POSCAR.environ, DERINFO.cluster, and DERINFO.environ

# Run the actual optimization
# ALERT: must use same number of processors (in here, 16) as the reference calculation !!!
# ALERT: must turn off symmetry (set ISYM=0 in INCAR) in the actual optimization !!!
./soptimize.sh -np 16 -walltime 1:00:00 cl2_test

# Hope your embedding calculation works :-)


#################
# params in input
#################
This section describes the important inputs in the VASP embedding run

***** In the reference calculations, watch out the following parameters:

--- INCAR file ---
LEXTPOT = .True. # This turns on the embedding calculation, must be set!

FCORR = 0.70     # This defines the region for derivative corrections. The exact derivatives are computed
                 # within the radius of FCORR*PSDMAX (PSDMAX determined by POTCAR)
                 # Usually, 0.70 is fine, default is 0.67, unless see difficulties in optimization, there 
                 # is no need to use larger values.

NGXF = 100       # Defines the dimension of the grid usded to express density and potential.
NGYF = 100       # Note that for efficiecy consideration, VASP may change them to different values, always
NGZF = 100       # check the OUTCAR for the true dimensions that are used by VASP.
       

***** In the Vemb optimization, watch out the following parameters:

--- split_DERINFO.py ---
n_atom = 2       # total number of atoms, including cluster and environtment.
cluster = [1]    # the atom number list for cluster. ALERT: atom number start from 1, yes, even it's in python!
na_cluster = {'Cl':1} # a dictionary defining the number of atoms by elements within cluster, used to make POSCAR for cluster
                      # key   - element labels;
                      # value - numbers of this particular type of element

--- INCAR file ---
ISYM = 0         # Turn off the symmetry, MUST BE SET in the Vemb optimization! Not necessary in reference calculation
LEXTPOT = .True. # MUST SET!
FCORR = 0.70     # MUST be consistent with reference calculations
NGXF/NGYF/NGZF   # MUST be consistent with reference calculations

--- optimize.py ---
lambd = 1.0e-7   # Prefactors for penalty function, something you need to play with...
istart = 0       # 0 - start from scratch; 1 - restart from existing EXTPOT file

ALERT (again...): The number of processors used in Vemb optimization must be consistent with the reference calculation.

#########
# authors
#########

The original version of this code was developed by Kuang Yu in the group of Emily A Carter in Mechanical and Aerospace
Engineering at Princeton University. Please cite the original publication, DOI: 10.1063/1.4922260 by Yu, Libisch, Carter.

It was subsequently ported to the VASP6 development version and refactored by Johannes M Dieterich (Princeton University)
and Florian Libisch (Vienna University of Technology).

Extensions were developed and tests carried out by Eduardo Schiavo (Unina Naples) and Florian Libisch.

We are grateful to the VASP developers for accepting our code for inclusion in the regular VASP code base.

