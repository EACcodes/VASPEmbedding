#!/usr/bin/env python3.9
import sys
import os
import subprocess
import numpy
from subprocess import check_output
from math import *
from numpy import reshape
from numpy import fft
from numpy import sum
from numpy import meshgrid
from numpy import dot
from numpy import array
from scipy import average
from scipy import optimize
from scipy import random

lambd = 1.0e-7 # coefficients for penalty
istart = 0
environ = 'slurm' # 'pbs'

# read information from running environment
if environ == 'slurm':
    WD = subprocess.check_output(['pwd']).decode().rstrip("\n")

filepath = WD + '/runvasp.sh'
ofile = open(filepath, 'w')
# Generate vasp run script
print("""\
#!bin/bash

if [ -d cluster ]; then mv cluster/CHGCAR ./CHGCAR.cluster; mv cluster/WAVECAR ./WAVECAR.cluster; mv cluster/DERINFO ./DERINFO.cluster; rm -r cluster; fi
if [ -d environ ]; then mv environ/CHGCAR ./CHGCAR.environ; mv environ/WAVECAR ./WAVECAR.environ; mv environ/DERINFO ./DERINFO.environ; rm -r environ; fi

mkdir cluster
mkdir environ
cp POSCAR.cluster cluster/POSCAR
cp POSCAR.environ environ/POSCAR
cp POTCAR.cluster cluster/POTCAR
cp POTCAR.environ environ/POTCAR
mv DERINFO.cluster cluster/DERINFO
mv DERINFO.environ environ/DERINFO
mv CHGCAR.cluster cluster/CHGCAR
mv WAVECAR.cluster cluster/WAVECAR
mv WAVECAR.environ environ/WAVECAR
mv CHGCAR.environ environ/CHGCAR
cp INCAR.cluster cluster/INCAR 
cp KPOINTS.cluster cluster/KPOINTS
ln -s ../EXTPOT cluster/EXTPOT
cp INCAR.environ environ/INCAR
cp KPOINTS.environ environ/KPOINTS
ln -s ../EXTPOT environ/EXTPOT
""", file=ofile)

if environ == 'slurm':
    print("""\
# run VASP
cd cluster
srun -n $SLURM_NPROCS $VASP_EXEC > logfile  
cd ../environ
srun -n $SLURM_NPROCS $VASP_EXEC > logfile  
cd ..
""",file=ofile)

ofile.close()

# Read grid dimensions
ifile = open( 'INCAR.cluster', 'r' )
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
    ifile = open( ifn, 'r' )
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
    density = numpy.array(density)
    return density
rou_ref = read_CHGCAR( 'ref/DERIV' )

def print_grid( extpot, ofn ):
    itot = 0
    iprint = 0
    os.system( 'if [ -f %s ]; then rm %s; fi'%(ofn, ofn) )
    ofile = open( ofn, 'w' )
    print('%d %d %d'%(NGXF, NGYF, NGZF),file=ofile)
    for iz in range(NGZF):
        for iy in range(NGYF):
           for ix in range(NGXF):
                print('%15.8e'%extpot[itot],end=' ',file=ofile)
                iprint += 1
                if iprint%8 == 0:
                    print('',file=ofile)
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

def read_POSCAR( ifn ):
    ifile = open( ifn, 'r' )
    ifile.readline()
    line = ifile.readline()
    sfac = float( line )
    box = []
    for idim in range(3):
        line = ifile.readline()
        words = line.split()
        l = [ float(words[i]) for i in range(3) ]
        box.append(l)
    box = numpy.array(box)*sfac
    ifile.close()    
    return box

box = read_POSCAR( 'POSCAR.cluster' )
vtot=box[0,0]*box[1,1]*box[2,2]*(1/(0.52917721067**3))
   
print('Cell volume (bohr^3): %16.10e'%vtot)
print('Step      -W_p[eV]         |sum(drho_p)|[e]   rmsd_p(c)[a.u.]      -W[eV]           |sum(drho)|[e]    rmsd(c)[a.u.]       method    loopn')

eval_counter = 0

def Lagrangian( params, *args ):
    global eval_counter
    extpot = params
    # set EXTPOT file
    print_grid( extpot, 'EXTPOT' )
    # run vasp
    os.system( 'bash runvasp.sh' )    
    # read E_cluster (eV)
    line = subprocess.check_output(['grep E0 cluster/logfile'],shell = True).decode().rstrip("\n")
    tmp = line.split('E0=')[1]
    E_cluster = float(tmp.split()[0])
    # read E_environ (eV)
    line = subprocess.check_output(['grep E0 environ/logfile'],shell = True).decode().rstrip("\n")
    tmp = line.split('E0=')[1]
    E_environ = float(tmp.split()[0])
    # compute \int V(r)*(na+nb-nref) (eV)
    rou_cluster = read_CHGCAR( 'cluster/DERIV' )
    rou_environ = read_CHGCAR( 'environ/DERIV' )
    d_rou = rou_cluster + rou_environ - rou_ref
    rou_ref2 = rou_ref
    # remove the normalization differences
#    d_rou = d_rou - average(d_rou)
    if len(args) == 0 or args[0] == 0:
        W = E_cluster + E_environ - numpy.dot( extpot, rou_ref2 )/ntot
    elif args[0] == 1:
        W = E_cluster
    elif args[0] == 2:
        W = E_environ
    else:
        W = - numpy.dot( extpot, rou_ref2 )/ntot
    # the library lbfgs algorithm is minimizing the function, while we need to maximize W, so add negative sign
    W = -W 
    if len(args) == 0 or args[0] == 0:
        grad = -d_rou/ntot
    elif args[0] == 1:
        grad = -rou_cluster/ntot
    elif args[0] == 2:
        grad = -rou_environ/ntot
    else:
        grad = rou_ref2/ntot
    eval_counter += 1
    # without penalty
    W_o = W
    normgrad_o= sqrt(numpy.dot(grad,grad))
    rms_o=normgrad_o*(sqrt(ntot))/vtot
    # add Laplacian penalty function
    field = reshape(numpy.array(extpot),(NGXF,NGYF,NGZF),order='F')
    laplacian = laplacian3d(field,NGXF,NGYF,NGZF)
    W -= sum(field * laplacian) * lambd
    grad -= reshape(laplacian,(NGXF*NGYF*NGZF), order='F') * lambd*2.0
    # +1 in execution counter
    normgrad=sqrt(numpy.dot(grad,grad))
    rms=normgrad*(sqrt(ntot))/vtot
    print('%4d' '%19.10e' '%19.10e' '%19.10e' '%19.10e' '%19.10e' '%19.10e' '%8d' '%8d' %(eval_counter,W,normgrad,rms,W_o,normgrad_o,rms_o,method,loopn))
    sys.stdout.flush()
    print_grid(rou_cluster+rou_environ, 'DERIV.current')
    if method==1:
       return W, grad
    if method==2:
       return W

def read_pot( ifn ):
    ifile = open( ifn, 'r' )
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
    return numpy.array(extpot)


if istart == 0:
    extpot0 = numpy.array( [0.0 for i in range(ntot)] )
else:
    extpot0 = read_pot('EXTPOT')

loopn = 0
method = 1
W0, direct = Lagrangian( extpot0 )
os.system('mv DERIV.current DERIV.init')

rou_cluster = read_CHGCAR( 'cluster/DERIV' )
rou_environ = read_CHGCAR( 'environ/DERIV' )
d_rou = rou_cluster + rou_environ - rou_ref
file_name = 'drho.init'
print_grid(d_rou, file_name)

f_old = 0
while True:
    method = 1
    loopn += 1
    extpot0 = read_pot('EXTPOT')
    x,f,d = optimize.fmin_l_bfgs_b( Lagrangian, x0=extpot0, args=(), pgtol=1e-05)
    dW = abs(f-f_old)
    if dW <= 1e-5:
       print('------------------------------------------')
       print('converged abs(dW): %19.10e'%dW)
       print('after %6d LBFGS loop'%loopn)
       print('------------------------------------------')
       rou_cluster = read_CHGCAR( 'cluster/DERIV' )
       rou_environ = read_CHGCAR( 'environ/DERIV' )
       d_rou = rou_cluster + rou_environ - rou_ref
       file_name = 'drho.' + str(loopn)
       print_grid(d_rou, file_name)
       break
    f_old = f
    print_grid(x,'EXTPOT')

print_grid(x,'EXTPOT.final')
print(d['warnflag'])
if d['warnflag'] == 2:
    print(d['task'])
