#!/usr/bin/env python3
import twixtools, h5py
import numpy as np
import sys
def get_noise_for(file):
    multi_twix = twixtools.read_twix(file)

    noise_data = []

    for meas in multi_twix:
        if type(meas) is not dict:
            continue
        mdb_list = [x for x in meas['mdb'] if x.is_flag_set('NOISEADJSCAN')]
        if len(mdb_list) == 0:
            continue
        nLin = 1 + max([mdb.cLin for mdb in mdb_list])
        nAcq = 1 + max([mdb.cAcq for mdb in mdb_list])
        nCha, nCol = mdb_list[0].data.shape
        out = np.zeros([nAcq, nLin, nCha, nCol], dtype=np.complex64)
        for mdb in mdb_list:
            out[mdb.cAcq, mdb.cLin] = mdb.data
        noise_data.append(out)
    return np.asarray(noise_data)
	
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("USAGE: get_noise.py INPUT.dat OUTPUT.h5")
        exit(1)
    data = get_noise_for(sys.argv[1])
    if len(data) == 0:
        print("Error: No noise data found.")
        exit(1)
    print(data.shape)
    try:
        h5f = h5py.File(sys.argv[2], 'w')
        h5f.create_dataset('noise', data=data)
    finally:
        h5f.close()