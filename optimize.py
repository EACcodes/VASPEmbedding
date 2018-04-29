#!/usr/bin/python2.7
import sys
import os
import commands
from math import *
from numpy import reshape
from numpy import fft
from numpy import sum
from numpy import meshgrid
from scipy import dot
from scipy import array
from scipy import average
from scipy import optimize
from scipy import random

lambd = 1.0e-7 # coefficients for penalty
istart = 0
environ = 'slurm' # 'pbs'

# read information from running environment
if environ == 'pbs':
    nproc = int( commands.getoutput("wc -l $PBS_NODEFILE | awk '{print $1}'") )
    WD = commands.getoutput("echo $PBS_O_WORKDIR")
elif environ == 'slurm':
    nproc = int( commands.getoutput("echo $SLURM_NPROCS") )
    WD = commands.getoutput("echo $SLURM_SUBMIT_DIR")

ofile = file( WD+'/runvasp.sh', 'w' )
# Generate vasp run script
print >> ofile, """\
#!bin/bash
module purge
module load openmpi/intel-13.0/1.6.3/64
module load intel/13.0/64/13.0.1.117
module load intel-mkl/11.1/0/64
export VASP_VERSION=5.3.3

export VIADEV_USE_AFFINITY=0
VASP_EXEC=/INSTALLATION/PATH/vasp
if [ -d cluster ]; then cp cluster/CHGCAR ./CHGCAR.cluster; cp cluster/WAVECAR ./WAVECAR.cluster; rm -r cluster; fi
if [ -d environ ]; then cp environ/CHGCAR ./CHGCAR.environ; cp environ/WAVECAR ./WAVECAR.environ;  rm -r environ; fi
mkdir cluster
mkdir environ
cp POSCAR.cluster cluster/POSCAR
cp POSCAR.environ environ/POSCAR
cp -r DERINFO.cluster cluster/DERINFO
cp -r DERINFO.environ environ/DERINFO
cp CHGCAR.cluster cluster/CHGCAR
cp CHGCAR.environ environ/CHGCAR
cp WAVECAR.cluster cluster/WAVECAR
cp WAVECAR.environ environ/WAVECAR
cp INCAR KPOINTS EXTPOT POTCAR cluster
cp INCAR KPOINTS EXTPOT POTCAR environ
"""

if environ == 'slurm':
    print >> ofile, """\
# run VASP
cd cluster
srun -n $SLURM_NPROCS $VASP_EXEC > logfile
cd ../environ
srun -n $SLURM_NPROCS $VASP_EXEC > logfile
cd ..
"""
elif environ == 'pbs':
    print >> ofile, """\
cd cluster
srun -np `wc -l $PBS_NODEFILE | awk '{print $1}'` --mca btl ^tcp --bind-to-socket $VASP_EXEC > logfile
cd ../environ
mpiexec -np `wc -l $PBS_NODEFILE | awk '{print $1}'` --mca btl ^tcp --bind-to-socket $VASP_EXEC > logfile
cd ..
"""

ofile.close()

# Read grid dimensions
ifile = file( 'INCAR', 'r' )
NGXF = -1
NGYF = -1
NGZF = -1
LEXTPOT = '.False.'
ISYM = 1
#LAECHG = '.False.'
for line in ifile:
    words = line.split(';')
    for word in words:
        if 'NGXF' in word:
            NGXF = int(word.split('=')[1])
        elif 'NGYF' in word:
            NGYF = int(word.split('=')[1])
        elif 'NGZF' in word:
            NGZF = int(word.split('=')[1])
        elif 'LEXTPOT' in word:
            LEXTPOT = word.split('=')[1]
        elif 'ISYM' in word:
            ISYM = int(word.split('=')[1])
#        elif 'LAECHG' in word:
#            LAECHG = word.split('=')[1]
ifile.close()
# check the necessary settings in INCAR
if NGXF == -1 or NGYF == -1 or NGZF == -1:
    sys.exit('ERROR: Need to specify NGXF, NGYF and NGZF in INCAR, exit code.')
if LEXTPOT.split()[0][1] != 'T' and LEXTPOT.split()[0][1] != 't':
    sys.exit('ERROR: Need to set the LEXTPOT tag in INCAR.')
if ISYM != 0:
    sys.exit('ERROR: ISYM should always be set to 0')
#if LAECHG.split()[0][1] != 'T' and LAECHG.split()[0][1] != 't':
#    sys.exit('ERROR: Need to set the LAECHG tag in INCAR.')
ntot = NGXF*NGYF*NGZF
# read reference density
def read_CHGCAR( ifn ):
    ifile = file( ifn, 'r' )
    density = []
    iread = 0
    for line in ifile:
        words = line.split()
        if len(words) == 0 and iread == 0:
            iread = 1
            continue
        if iread > 0 and len(words) == 0:
            iread = 0
            continue
        if iread == 1:
            [ nx, ny, nz ] = [ int(words[i]) for i in range(3) ]
            if nx!= NGXF or ny!=NGYF or nz!=NGZF:
                sys.exit('ERROR in read_CHGCAR: Density grid dimension mismatch: file %s'%ifn)
            iread = 2
            continue
        elif iread == 2:
            for word in words:
                density.append(float(word))
    ifile.close()
    density = array(density)
    return density
rou_ref = read_CHGCAR( 'ref/DERIV' )

def print_grid( extpot, ofn ):
    itot = 0
    iprint = 0
    os.system( 'if [ -f %s ]; then rm %s; fi'%(ofn, ofn) )
    ofile = file( ofn, 'w' )
    print >>ofile, '%d %d %d'%(NGXF, NGYF, NGZF)
    for iz in range(NGZF):
        for iy in range(NGYF):
           for ix in range(NGXF):
                print >>ofile, '%15.8e'%extpot[itot],
                iprint += 1
                if iprint%8 == 0:
                    print >>ofile
                itot += 1
    ofile.close()

def laplacian3d(field, nx, ny, nz):
    kx = fft.fftfreq(nx,1.0)*2*pi
    ky = fft.fftfreq(ny,1.0)*2*pi
    kz = fft.fftfreq(nz,1.0)*2*pi
    KX, KY, KZ = meshgrid(kx, ky, kz, indexing='ij')
    return fft.ifft(-KX**2*fft.fft(field, axis = 0), axis = 0).real + \
        fft.ifft(-KY**2*fft.fft(field, axis = 1), axis = 1).real + \
        fft.ifft(-KZ**2*fft.fft(field, axis = 2), axis = 2).real

eval_counter = 0
def Lagrangian( params, *args ):
    global eval_counter
    extpot = params
    # set EXTPOT file
    print_grid( extpot, 'EXTPOT' )
    # run vasp
    os.system( 'bash runvasp.sh' )    
    # read E_cluster (eV)
    line = commands.getoutput( 'grep E0 cluster/logfile' )
    tmp = line.split('E0=')[1]
    E_cluster = float(tmp.split()[0])
    # read E_environ (eV)
    line = commands.getoutput( 'grep E0 environ/logfile' )
    tmp = line.split('E0=')[1]
    E_environ = float(tmp.split()[0])
    # check for unfinished VASP jobs
    ierr = int(commands.getoutput( 'grep "General timing and accounting informations" cluster/OUTCAR | wc' ).split()[0])
    if ierr == 0:
        sys.exit('ERROR: VASP jobs exit abnormally in cluster calculation!')
    ierr = int(commands.getoutput( 'grep "General timing and accounting informations" environ/OUTCAR | wc' ).split()[0])
    if ierr == 0:
        sys.exit('ERROR: VASP jobs exit abnormally in environ calculation!')
    ierr = int(commands.getoutput( 'grep "Internal error, mismatch with DERINFO" cluster/logfile | wc' ).split()[0])
    if ierr > 0:
        sys.exit('ERROR: DERINFO error in cluster calculation, please check for possible file corruption.')
    ierr = int(commands.getoutput( 'grep "Internal error, mismatch with DERINFO" environ/logfile | wc' ).split()[0])
    if ierr > 0:
        sys.exit('ERROR: DERINFO error in environ calculation, please check for possible file corruption.')
    # compute \int V(r)*(na+nb-nref) (eV)
    rou_cluster = read_CHGCAR( 'cluster/DERIV' )
    rou_environ = read_CHGCAR( 'environ/DERIV' )
    d_rou = rou_cluster + rou_environ - rou_ref
    # remove the normalization differences
#    d_rou = d_rou - average(d_rou)
    if len(args) == 0 or args[0] == 0:
        W = E_cluster + E_environ - dot( extpot, rou_ref )/ntot
    elif args[0] == 1:
        W = E_cluster
    elif args[0] == 2:
        W = E_environ
    else:
        W = - dot( extpot, rou_ref )/ntot
    # the library lbfgs algorithm is minimizing the function, while we need to maximize W, so add negative sign
    W = -W 
    if len(args) == 0 or args[0] == 0:
        grad = -d_rou/ntot
    elif args[0] == 1:
        grad = -rou_cluster/ntot
    elif args[0] == 2:
        grad = -rou_environ/ntot
    else:
        grad = rou_ref/ntot
    eval_counter += 1
    print '--------------------'
    print ' Evaluation No: %d'%eval_counter
    print ' Norm of gradient without penalty: %12.6e'%sqrt(dot(grad,grad))
    print ' Lagrangian value without penalty: %12.6e'%W
    # add Laplacian penalty function
    field = reshape(array(extpot),(NGXF,NGYF,NGZF),order='F')
    laplacian = laplacian3d(field,NGXF,NGYF,NGZF)
    W -= sum(field * laplacian) * lambd
    grad -= reshape(laplacian,(NGXF*NGYF*NGZF), order='F') * lambd*2.0
    # +1 in execution counter
    print ' Norm of gradient: %12.6e'%sqrt(dot(grad,grad))
    print ' Lagrangian value: %12.6e'%W
    sys.stdout.flush()
    print_grid(rou_cluster+rou_environ, 'DERIV.current')
    return W, grad

def read_pot( ifn ):
    ifile = file( ifn, 'r' )
    iread = 0
    extpot = []
    for line in ifile:
        words = line.split()
        if iread == 0:
            [ nx, ny, nz ] = [ int(words[i]) for i in range(3) ]
            if nx!= NGXF or ny!=NGYF or nz!=NGZF:
                sys.exit('ERROR in read_pot: Density grid dimension mismatch: file %s'%ifn)
            iread = 1
            continue
        elif iread == 1:
            for word in words:
                extpot.append(float(word))
    ifile.close()
    return array(extpot)

# Initial guess, allocating potential array
if istart == 0:
    extpot0 = array( [0.0 for i in range(ntot)] )
else:
    extpot0 = read_pot('EXTPOT')
#direct = read_pot('direct')

# identify direction
#W0, direct = Lagrangian( extpot0, (0) )
#direct = array([0.0 for i in range(ntot)])
#direct[405050] = 1.0
#
#W0, grad0 = Lagrangian( extpot0, (0) )
#
#delta = 1e0 # 1e4
#ofile = file( 'logfile', 'w' )
#print >> ofile, '# %15.8e'%dot(grad0, direct)
#for i in range(-2,3):
#    extpot = extpot0 + delta*i*direct
#    W, grad = Lagrangian( extpot, (0) )
#    print >> ofile, '%15.6e%20.8e'%(i*delta,W)
#ofile.close()

W0, direct = Lagrangian( extpot0 )
os.system('mv DERIV.current DERIV.init')
x,f,d = optimize.fmin_l_bfgs_b( Lagrangian, x0=extpot0, args=(), pgtol=1e-06 )

print_grid(x,'EXTPOT.final')
print d['warnflag']
if d['warnflag'] == 2:
    print d['task']
