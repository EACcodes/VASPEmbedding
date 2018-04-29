#!/usr/bin/python
import sys

ifile = file( sys.argv[1], 'r' )

iread = 1
ix = 0
iy = 0
iz = 0
for line in ifile:
    words = line.split()
    if len(words) == 0:
        iread = 1
        continue
    elif iread == 1:
        iread = 2
        [ nx, ny, nz ] = [ int(words[i]) for i in range(3) ]
        continue
    if iread == 2:
        for word in words:
            val = float(word)
            if ix == 50 and iy == 50:
                print '%5d%20.8e'%(iz, val)
            ix += 1
            if ix%nx == 0:
                ix = 0
                iy += 1
                if iy%ny == 0:
                    iy = 0
                    iz += 1
ifile.close()
