#!/bin/python3

import numpy as np
import os
import sys



f = open( "list.out", 'rb' )
N = np.fromfile( f, dtype=np.int64, count=1 )
id_size = np.fromfile( f, dtype=np.int32, count=1 )

if id_size == 4:
    ptype = np.dtype([('ID','u4'), ('fof', 'i4'), ('sub', 'i4')])
else :
    ptype = np.dtype([('ID','u8'), ('fof', 'i4'), ('sub', 'i4')])


particles = np.fromfile( f, dtype=ptype )

print(particles[0:(10 if (N>10) else N)])

f.close()
