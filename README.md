# Cloudmrhub


# Installation
```
pip install git+https://github.com/cloudmrhub-com/cloudmrhub.git

```

# Suggestions
```
#create an environment 
python3 -m venv CMT
source CMT/bin/activate
pip install git+https://github.com/cloudmrhub-com/cloudmrhub.git
```
# Example
Examples can be found in the [git test directoy](https://github.com/cloudmrhub-com/cloudmrhub/tree/main/cloudmrhub/test)

The codes reads a signal and a noise kspace from siemens scanner and computes reconstructions and snr for RSS and B1

```
from cloudmrhub.cm2D import cm2DKellmanRSS, cm2DReconRSS, cm2DKellmanB1,cm2DReconB1
import cloudmrhub.cm as cm

import twixtools
import numpy as np
import matplotlib.pyplot as plt


import sys
s= sys.argv[1]
n=sys.argv[2]

# read the files
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

# instantiate th classes
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

plt.imshow(np.abs(im))
plt.colorbar()
plt.title('RSS Recon')

L=cm2DKellmanRSS()
L.setSignalKSpace(S)
L.setNoiseKSpace(N)
snr=np.abs(L.getOutput())


plt.imshow(snr)
plt.colorbar()
plt.title('RSS SNR')



```



# Example SNR
```
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

L2=cm2D.cm2DSignalToNoiseRatioPseudoMultipleReplicas()
L=cm2D.cm2DReconRSS()
L2.numberOfReplicas=100
L.setNoiseKSpace(N)
L.setSignalKSpace(S)
L2.reconstructor=L

plt.subplot(221)
plt.imshow(np.abs(L2.getOutput()))
plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
plt.title('PMR')


R2=cm2D.cm2DSignalToNoiseRatioPseudoMultipleReplicasWein()
R=cm2D.cm2DReconRSS()
R2.numberOfReplicas=3
R2.boxSize=10
R.setNoiseKSpace(N)
R.setSignalKSpace(S)
R2.reconstructor=R



plt.subplot(222)
plt.imshow(np.abs(R2.getOutput()))
plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
plt.title('CR')


R=cm2D.cm2DKellmanRSS()
R.setNoiseKSpace(N)
R.setSignalKSpace(S)
plt.subplot(223)
plt.imshow(np.abs(R.getOutput()))
plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
plt.title('Kellman')

plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
plt.colorbar(plt.gcf().axes[0].images[0], ax=plt.gcf().axes)
# add a title for the figure
plt.suptitle("SNR", fontsize=16)


plt.show()
```
# Roadmap
    [x] SNR
        [x] Kellman
            [x] RSS
            [x] B1
            [x] Sense

        [x] Pseudo Multiple Replicas 
            [x] RSS
            [x] B1
            [x] Sense
            [x] Grappa

        [x] Pseudo Multiple Replicas Wien
            [x] RSS
            [x] B1
            [x] Sense
            [x] Grappa

        [ ] Multiple Replicas

    [x] Recon

        [x] RSS
        [x] B1
        [x] Sense
        [x] Grappa
        

    [ ] Coilsensitivity Maps:

        [x] Inner
        [x] BodyCoil
        [ ] Espirit
    


# Cite Us

- Montin E, Lattanzi R. Seeking a Widely Adoptable Practical Standard to Estimate Signal-to-Noise Ratio in Magnetic Resonance Imaging for Multiple-Coil Reconstructions. J Magn Reson Imaging. 2021 Dec;54(6):1952-1964. doi: 10.1002/jmri.27816. Epub 2021 Jul 4. PMID: 34219312; PMCID: PMC8633048.


[*Dr. Eros Montin, PhD*](http://me.biodimensional.com)

**46&2 just ahead of me!**
