#!/usr/bin/env python2.7
import sys
import numpy as np

ifile = file( sys.argv[1], 'r' )

iread = 0
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
        data = np.zeros([nx, ny, nz])
        continue
    if iread == 2:
        for word in words:
            val = float(word)
            data[ix, iy, iz] = val
            # if ix == 50 and iy == 50:
            #     print '%5d%20.8e'%(iz, val)
            # if iy == 50 and iz == 50:
            #     print '%5d%20.8e'%(ix, val)
            # if ix == 50 and iz == 50:
            #     print '%5d%20.8e'%(iy, val)
            ix += 1
            if ix%nx == 0:
                ix = 0
                iy += 1
                if iy%ny == 0:
                    iy = 0
                    iz += 1
ifile.close()

for v in data[50, 50, :]:
    print v
