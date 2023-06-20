from cloudmrhub.cm2D import cm2DKellmanRSS, cm2DReconRSS
import raider.twixtools
import numpy as np
import matplotlib.pyplot as plt

s='../data/meas_MID00024_FID188178_Multislice.dat'
n='../data/meas_MID00027_FID188181_Multislice_no_RF.dat'

twix=twixtools.map_twix(s)
im_array = twix[0]['image']
im_array.flags['remove_os'] = True  # activate automatic os removal
im_array.flags['average']['Rep'] = True  # average all repetitions
im_array.flags['average']['Ave'] = True  # average all repetitions
S=np.transpose(im_array[0,0,0,0,0,0,0,0,0,0,0,0,0,:,:,:],[2,0,1])


twix=twixtools.map_twix(n)
im_array = twix[0]['image']
im_array.flags['remove_os'] = False  # activate automatic os removal
im_array.flags['average']['Rep'] = False  # average all repetitions
im_array.flags['average']['Ave'] = False # average all repetitions
N=np.transpose(im_array[0,0,0,0,0,0,0,0,0,0,0,0,0,:,:,:],[2,0,1])

# N=N[0:2,0:3,0:4]
# S=S[0:2,0:3,0:4]

import cm


L=cm2DReconRSS()
L.setSignalKSpace(S)
L.setNoiseKSpace(N)

NBW=L.getNoiseBandWidth()
print(NBW)
nc2=L.getNoiseCovariance()
plt.figure()
plt.imshow(np.abs(nc2))
plt.show()
p=L.getPrewhitenedSignal()
plt.figure()
plt.imshow(np.abs(p[:,:,0]))
plt.show()



im=L.getOutput()

import matplotlib.pyplot as plt

plt.imshow(im)
plt.show()

L=cm2DKellmanRSS()
L.setSignalKSpace(S)
L.setNoiseKSpace(N)
snr=L.getOutput()


plt.imshow(snr)
plt.show()

