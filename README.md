# Cloudmrhub


# Installation
```
pip install git+https://github.com/cloudmrhub-com/CMRCode.git

```

# Suggestions
```
#create an nevironment 
python3 -m venv CMT
source CMT/bin/activate
pip install git+https://github.com/cloudmrhub-com/CMRCode.git
```
# Example
Examples can be found in the [git test directoy](https://github.com/cloudmrhub-com/CMRCode/tree/main/cloudmrhub/test)


```
from cloudmrhub.cm2D import cm2DKellmanRSS, cm2DReconRSS, cm2DKellmanB1,cm2DReconB1
import cloudmrhub.cm as cm

import twixtools
import numpy as np
import matplotlib.pyplot as plt


import sys
s= sys.argv[1]
n=sys.argv[2]

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



L=cm2DReconRSS()
L.setSignalKSpace(S)
L.setNoiseKSpace(N)

NBW=L.getNoiseBandWidth()
print(NBW)
nc2=L.getNoiseCovariance()
plt.figure()
plt.imshow(np.abs(nc2))
plt.colorbar()
plt.title(f'Noise Covariance Matrix BW:{NBW:03f}')
p=L.getPrewhitenedSignal()
plt.figure()
plt.imshow(np.abs(p[:,:,0]))
plt.colorbar()
plt.title('Prewhitened Kspace First Coil')



im=L.getOutput()

import matplotlib.pyplot as plt

plt.imshow(im)
plt.colorbar()
plt.title('RSS Recon')

L=cm2DKellmanRSS()
L.setSignalKSpace(S)
L.setNoiseKSpace(N)
snr=L.getOutput()


plt.imshow(snr)
plt.colorbar()
plt.title('RSS SNR')



L=cm2DReconB1()
L.setSignalKSpace(S)
L.setNoiseKSpace(N)
L.setCoilSensitivityMatrixSource(S)
L.setCoilSensitivityMatrixCalculationMethod('inner')
for t in range(S.shape[-1]):
    NN=np.sqrt(S.shape[-1])
    plt.subplot(int(NN),int(NN),t+1)
    plt.imshow(np.abs(cm.calculate_simple_sense_sensitivitymaps(S,'ref')[:,:,t]))
    plt.title(f'Coil {t}')

plt.figure()
plt.imshow(L.getOutput())
plt.title('B1 Recon')

L=cm2DKellmanB1()
L.setSignalKSpace(S)
L.setNoiseKSpace(N)
L.setCoilSensitivityMatrixSource(S)
L.setCoilSensitivityMatrixCalculationMethod('inner')
plt.figure()
plt.imshow(L.getOutput())
plt.colorbar()
plt.title('SNR B1')




plt.show()
```
# Roadmap

[ ] SNR
    [ ] Kellman
        [x] RSS v0
        [x] B1  v1
        [ ] Sense v2

    [ ] Pseudo Multiple Replicas 
        [x] RSS v1
        [ ] B1  v1
        [ ] Sense v2
        [ ] Grappa v2
        [ ] Adapt v3

    [ ] Pseudo Multiple Replicas Wien
        [x] RSS v1
        [ ] B1 v1
        [ ] Sense v2
        [ ] Grappa v2
        [ ] Adapt v3
        [ ] CompressSenseBart v4

[ ] Recon
    [x] RSS v0
    [x] B1 v1
    [ ] Sense v2
    [ ] Grappa
    [ ] Adapt v3
    [ ] Compress Sense Bart v4

[ ] Coilsensitivity Maps:
    [x] Inner v1
    [ ] Adapt v3
    [ ] BodyCoil v3
    [ ] Espirit
    



- v1:
    - recon RSS,B1
    - SNR Kellman, RSS,B1
    - SNR Pseudo MR  RSS,B1
- v0
    - recon RSS

# Cite Us

- Montin E, Lattanzi R. Seeking a Widely Adoptable Practical Standard to Estimate Signal-to-Noise Ratio in Magnetic Resonance Imaging for Multiple-Coil Reconstructions. J Magn Reson Imaging. 2021 Dec;54(6):1952-1964. doi: 10.1002/jmri.27816. Epub 2021 Jul 4. PMID: 34219312; PMCID: PMC8633048.


[*Dr. Eros Montin, PhD*](http://me.biodimensional.com)

**46&2 just ahead of me!**
