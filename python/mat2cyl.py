import numpy as np
import scipy.io
import sys

cyl = scipy.io.loadmat(str(sys.argv[1]))

LEN = len(cyl['Sta'])
arr = np.zeros((LEN, 14))

arr[:, 0] = cyl['Rad'].reshape(-1, LEN)
arr[:, 1] = cyl['Len'].reshape(-1, LEN)
arr[:, 2:5] = cyl['Sta']
arr[:, 5:8] = cyl['Axe']
arr[:, 8] = cyl['CPar'].reshape(-1, LEN)
arr[:, 9] = cyl['CExt'].reshape(-1, LEN)
# arr[:, 10] = cyl['CExt'].reshape(-1, LEN)
arr[:, 11] = cyl['BoC'][:, 1]
# arr[:, 12] = cyl['CExt'].reshape(-1, LEN)
# arr[:, 13] = cyl['CExt'].reshape(-1, LEN)

np.savetxt(sys.argv[1][:-4] + '.cyl', arr)