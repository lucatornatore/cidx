#!/bin/python3

import numpy as np
import os
import sys


# -----------------------------------------------------
# id_size is either 4 or 8
id_size = 4



# -----------------------------------------------------
# let's say that the original id list is a python list
id_list = [102345, 67895423, 34567123, 1234567890, 2034987654]

# then, convert it into a numpy array
if id_size == 4:
    id_array = np.array(id_list, dtype='<u4')
else
    id_array = np.array(id_list, dtype='<u8')

header_t = np.dtype([('id_size','u4'), ('N', 'u8')])
header = np.ndarray(shape=(1,1), dtype=header_t)
header[0]=[(id_size, id_array.shape[0])]

# note: id_array.shape[0] is the size of the list

# ----------------------------------------------------
# now write the input id-list file to be processed by
# cidx

f = open( "list.in", 'wb' )

header.tofile(f)
np.tofile( f )

f.close()
