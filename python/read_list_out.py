#!/bin/python3

import numpy as np
import os
import sys



f = open( "list.out", 'rb' )
N = np.fromfile( f, dtype=np.int64, count=1 )
id_size = np.fromfile( f, dtype=np.int32, count=1 )

if id_size == 4:
    ptype = np.dtype([('ID','u4'), ('fof', 'u4'), ('sub', 'u4')])
else :
    ptype = np.dtype([('ID','u8'), ('fof', 'u4'), ('sub', 'u4')])


particles = np.fromfile( f, dtype=ptype )

print(particles[0:(10 if (N>10) else N)])

f.close()
