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
import twixtools
import numpy as np
import matplotlib.pyplot as plt

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


L=cm2DKellmanRSS()
L.setSignalKSpace(S)
L.setNoiseKSpace(N)
snr=L.getOutput()


plt.imshow(snr)
plt.colorbar()
plt.title('RSS SNR')
plt.show()

```

# Roadmap

[ ] SNR
    [ ] Kellman
        [x] RSS v0
        [x] B1  v1
        [ ] Sense v2
    [ ] Pseudo Multiple Replicas 
        [ ] RSS v1
        [ ] B1  v1
        [ ] Sense v2
        [ ] Grappa v2
        [ ] Adapt v3

    [ ] Pseudo Multiple Replicas Wien
        [ ] RSS v1
        [ ] B1 v1
        [ ] Sense v2
        [ ] Grappa v2
        [ ] Adapt v3

[ ] Recon
    [x] RSS v0
    [x] B1 v1
    [ ] Sense v2
    [ ] Adapt v3

[ ] Coilsensitivity Maps:
    [x] Inner v1
    [ ] Adapt v3
    [ ] BodyCoil v2
    



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
