#!/usr/bin/env python3.9
import sys
import subprocess
import os
n_atom = 100
cluster = [6,8,16,18,23,25,26,32,33,35,38,43,47,48,49]
na_cluster = {'Pt':15}
#please remove the line "Selective dynamics" from the POSCAR of ref if present

environ = list(range(1,n_atom+1))
for i in cluster:
    environ.remove(i)

def split_DERINFO( ifn, ofn_cluster, ofn_environ ):
    ifile = open( ifn, 'r' )
    i_atom_cluster = 0
    i_atom_environ = 0
    i_atom_prev = 0
    ofile1 = open( ofn_cluster, 'w' )
    ofile2 = open( ofn_environ, 'w' )
    outtag = 0
    for line in ifile:
        words = line.split()
        i_atom = int(line[:4])
#       i_atom = int(words[0])
        if i_atom != i_atom_prev: # new atom
            i_atom_prev = i_atom
            if i_atom in cluster: # new atom in cluster
                outtag = 1
                i_atom_cluster += 1
            elif i_atom in environ:
                outtag = 2
                i_atom_environ += 1
            else:
                sys.exit('Error!')
        if outtag == 1:
            print('%4d%s'%(i_atom_cluster,line[4:]),end="",file=ofile1)
        else:
            print('%4d%s'%(i_atom_environ,line[4:]),end="",file=ofile2)
    ifile.close()
    ofile1.close()
    ofile2.close()

os.system('mkdir DERINFO.cluster DERINFO.environ')
nproc = int( subprocess.check_output(['ls -l ref/DERINFO/* | wc | cut -c -10'],shell = True).decode().rstrip("\n"))
for iproc in range(1, nproc+1):
    split_DERINFO('ref/DERINFO/%d'%iproc, 'DERINFO.cluster/%d'%iproc, 'DERINFO.environ/%d'%iproc)

# create POSCARs
ofile1 = open( 'POSCAR.cluster', 'w' )
ofile2 = open( 'POSCAR.environ', 'w' )
ifile = open( 'ref/POSCAR', 'r' )
for i in range(5):
    line = ifile.readline()
    print(line,end="",file=ofile1)
    print(line,end="",file=ofile2)
line = ifile.readline()
atypes = line.split()
line = ifile.readline()
words = line.split()
n_atoms = [ int(words[i]) for i in range(len(atypes)) ]
na_environ = {}
for itype in range(len(atypes)):
    atype = atypes[itype]
    na_environ[atype] = n_atoms[itype]

for itype in range(len(atypes)):
    atype = atypes[itype]
    if atype in na_cluster:
        na_environ[atype] -= na_cluster[atype]

for itype in range(len(atypes)):
    atype = atypes[itype]
    if atype in na_cluster and na_cluster[atype] > 0:
        print('%4s'%atype,end="",file=ofile1)
    if atype in na_environ and na_environ[atype] > 0:
        print('%4s'%atype,end="",file=ofile2)
print('',file=ofile1)
print('',file=ofile2)

for itype in range(len(atypes)):
    atype = atypes[itype]
    if atype in na_cluster and na_cluster[atype] > 0:
        print('%4d'%na_cluster[atype],end="",file=ofile1)
    if atype in na_environ and na_environ[atype] > 0:
        print('%4d'%na_environ[atype],end="",file=ofile2)
print('',file=ofile1)
print('',file=ofile2)

line = ifile.readline()
print(line,end="",file=ofile1)
print(line,end="",file=ofile2)

na_tot = sum(n_atoms)
for i_atom in range(na_tot):
    line = ifile.readline()
    if i_atom+1 in cluster:
        print(line,end="",file=ofile1)
    else:
        print(line,end="",file=ofile2)

ofile1.close()
ofile2.close()
ifile.close()
