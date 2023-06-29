import numpy as np
import cloudmrhub.cm2D as cm2D

import matplotlib.pyplot as plt
import  twixtools

s='cloudmrhub/data/meas_MID00024_FID188178_Multislice.dat'
n='cloudmrhub/data/meas_MID00027_FID188181_Multislice_no_RF.dat'

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
import devcm2D as cm2D
from pynico_eros_montin import pynico as pn
import devcm as cm
import numpy as np
import matplotlib.pyplot as plt


FA=1
PA=2
ACL=20
L=cm2D.cm2DReconSense()
L.AccelerationF=FA
L.AccelerationP=PA

US=cm.undersample2DDataSENSE(S,frequencyacceleration=FA,phaseacceleration=PA)

L.setSignalKSpace(US)
L.setNoiseKSpace(N)
L.setCoilSensitivityMatrixSource(S)
L.setCoilSensitivityMatrixCalculationMethod('inner')

plt.figure()
plt.imshow(np.abs(L.getOutput()))
plt.colorbar()
plt.title('recon')


# transform in Kellman!!!
L.__class__=cm2D.cm2DKellmanSense

plt.figure()
plt.imshow(np.abs(L.getOutput()))
plt.colorbar()
plt.title('SNR')

# transform in Kellman!!!
L.__class__=cm2D.cm2DGfactorSense

plt.figure()
plt.imshow(np.abs(L.getOutput()))
plt.colorbar()
plt.title('GFactor')
plt.show()
