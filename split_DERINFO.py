#!/usr/bin/env python2.7
import sys
import commands
import os

n_atom = 2
cluster = [1]
na_cluster = {'Cl':1}
environ = range(1,n_atom+1)
for i in cluster:
    environ.remove(i)

def split_DERINFO( ifn, ofn_cluster, ofn_environ ):
    ifile = file( ifn, 'r' )
    i_atom_cluster = 0
    i_atom_environ = 0
    i_atom_prev = 0
    ofile1 = file( ofn_cluster, 'w' )
    ofile2 = file( ofn_environ, 'w' )
    outtag = 0
    for line in ifile:
        words = line.split()
        i_atom = int(words[0])
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
            print >> ofile1, '%4d%s'%(i_atom_cluster,line[4:]),
        else:
            print >> ofile2, '%4d%s'%(i_atom_environ,line[4:]),
    ifile.close()
    ofile1.close()
    ofile2.close()

os.system('mkdir DERINFO.cluster DERINFO.environ')
nproc = int( commands.getoutput('ls -l ref/DERINFO/* | wc | cut -c -10') )
for iproc in range(1, nproc+1):
    split_DERINFO('ref/DERINFO/%d'%iproc, 'DERINFO.cluster/%d'%iproc, 'DERINFO.environ/%d'%iproc)

# create POSCARs
ofile1 = file( 'POSCAR.cluster', 'w' )
ofile2 = file( 'POSCAR.environ', 'w' )
ifile = file( 'ref/POSCAR', 'r' )
for i in range(5):
    line = ifile.readline()
    print >> ofile1, line,
    print >> ofile2, line,
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
        print >>ofile1, '%4s'%atype,
    if atype in na_environ and na_environ[atype] > 0:
        print >>ofile2, '%4s'%atype,
print >>ofile1, ''
print >>ofile2, ''

for itype in range(len(atypes)):
    atype = atypes[itype]
    if atype in na_cluster and na_cluster[atype] > 0:
        print >>ofile1, '%4d'%na_cluster[atype],
    if atype in na_environ and na_environ[atype] > 0:
        print >>ofile2, '%4d'%na_environ[atype],
print >>ofile1, ''
print >>ofile2, ''

line = ifile.readline()
print >>ofile1, line,
print >>ofile2, line,

na_tot = sum(n_atoms)
for i_atom in range(na_tot):
    line = ifile.readline()
    if i_atom+1 in cluster:
        print >> ofile1, line,
    else:
        print >> ofile2, line,

ofile1.close()
ofile2.close()
ifile.close()
