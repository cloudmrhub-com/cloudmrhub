import twixtools
import json
from pynico_eros_montin import pynico as pn
import numpy as np

class k2d:
    def __init__(self,k=None):
        self.k=np.array([])
        if k is not None:
            self.set(k)
    def getSize(self):
        return self.k.shape[0:2]
    def getNCoils(self):
        return self.k.shape[-1]
    def get(self):
        return self.k
    def set(self,k):
        self.k =k
    def isEmpty(self):
        return np.sum(self.k.shape)>0

    
        
        
class cmOutput:
    def __init__(self,message=None):
        self.Exporter = []
        self.Log = pn.Log(message)
        self.Type = "general"
        self.SubType = "zero"
        self.OUTPUTLOGFILENAME = pn.Pathale()
        self.OUTPUTFILENAME = pn.Pathale()

    def addToExporter(self, type, name, value):
        self.Exporter.append([type, name, float(value)])

    def setOutputLogFileName(self, L):
        self.OUTPUTLOGFILENAME = pn.Pathable(L)

    def getOutputLogFileName(self):
        return self.OUTPUTLOGFILENAME.getPosition()

    def setOutputFileName(self, L):
        self.OUTPUTFILENAME = pn.Pathable(L)

    def getOutputFileName(self):
        return self.OUTPUTFILENAME.getPosition()

    def setLog(self, L):
        self.Log = L

    def getLog(self):
        return self.Log

    def setTypeOutput(self, L):
        self.Type = L

    def getTypeOutput(self):
        return self.Type

    def appendLog(self, L):
        self.Log.append(L)

    def logIt(self, W, t=None):
        
        self.Log.append(W,t)

    def getResults(self):
        return self.exportResults()

    # def getResultsAsImages(self):
    #     return self.getImagesFromResuls(self.getResults())

    # def exportResults(self, fn):
    #     O = {}
    #     O["version"] = "20201223"
    #     O["author"] = "eros.montin@gmail.com"
    #     if self.Type == "":for t in range(len(self.Log)):
    #         print(self.Log[t])
    #         O["type"] = "DATA"
    #     else:
    #         O["type"] = self.Type
    #     if self.SubType == "":
    #         O["subtype"] = ""
    #     else:
    #         O["subtype"] = self.SubType
    #     for t in range(len(self.Exporter)):
    #         if self.Exporter[t][0] == "image2D":
    #             im = {}
    #             im["slice"] = self.image2DToJson(self.Exporter[t][2])
    #             im["imageName"] = self.Exporter[t][1]
    #             O["images"].append(im)
    #     if len(fn) > 0:
    #         with open(fn, "w") as f:
    #             json.dump(O, f)
    #     else:
    #         with open(self.getOutputFileName(), "w") as f:
    #             json.dump(O, f)

    def exportLog(self, fn=None):
        if fn==None:
            fn=self.getOutputFileName()
        
        self.Log.writeLogAs(fn)


    def whatHappened(self):
        print(self.Log.getWhatHappened())

    def errorMessage(self):
        self.logIt("ERROR", "error")

    def outputError(self, fn):
        self.errorMessage()
        self.exportLog(fn)



def prewhiteningSignal(signalrawdata, psi):
    L = np.linalg.cholesky(psi)
    L_inv = np.linalg.inv(L)
    nc = signalrawdata.shape[-1]
    nf = signalrawdata.shape[0]
    nph = signalrawdata.shape[1]
    pw_signalrawdata = np.transpose(signalrawdata, [2, 1, 0])  # [ncoil nphase nfreq]
    shape= pw_signalrawdata.shape
    output = pw_signalrawdata.reshape((shape[0], shape[1]*shape[2]))
    pw_signalrawdata = np.matmul(L_inv, output)
    pw_signalrawdata = np.reshape(pw_signalrawdata, [nc, nph, nf])
    pw_signalrawdata = np.transpose(pw_signalrawdata, [2, 1, 0])  # [nfreq nphase ncoil]
    return pw_signalrawdata

def calculate_covariance_matrix(noise, bandwidth=1):
    """Calculates the covariance matrix of noise.
       update 09/24/2020 Prof. Riccardo Lattanzi, Ph.D.

    Args:
        noise: The noise as a NumPy array.
        bandwidth: The bandwidth of the noise.
    Returns:
        The covariance matrix of the noise.
    """
    
    nf,nph,nchan = noise.shape
    Rn = np.zeros((nchan, nchan))
    noise_samples = np.swapaxes(noise.reshape((nph*nf,nchan)),0,1) # non complex conj
    Rn = ( noise_samples @ noise_samples.conj().T)/ (2.0 * float(nf)*float(nph))
    Rn = Rn / float(bandwidth)
    return Rn

import PIL

def getMRoptimumTestData(nc=None):
    """Create some fake data to test reconstructions
    Args:
        path: The path to the image file.
    Returns:
        Kspace singal data and nc.
    """
    image = PIL.Image.open('./eros.jpg')
    K = MRfft(np.array(image),[0,1])
    
    # Create a fake noise covariance matrix
    if nc is None:
        nc = np.eye(3) * np.random.rand(3) / 1000000

    KN=create_fake_noise_kspace_data(shape=K.shape,correlation_matrix=nc)
    return K, KN, nc
    


def get_pseudo_noise(msize, corr_noise_factor):
    """
    msize (freq,phase,coil)
    """
    N = np.sqrt(0.5) * (np.random.randn(*msize) + 1j * np.random.randn(*msize))
    nrow = msize[0]
    ncol = msize[1]
    nchan = msize[2]
    gaussian_whitenoise = np.reshape(N, (nrow * ncol, nchan))
    gaussian_whitenoise = np.matmul(corr_noise_factor,gaussian_whitenoise.T)
    gaussian_whitenoise = np.reshape(gaussian_whitenoise.T, (nrow, ncol, nchan))

    return gaussian_whitenoise

def create_fake_noise_kspace_data(shape, correlation_matrix=None):
    """Creates a fake noise matrix with some correlations between the channels.

    Args:
        shape: The shape of the matrix.
        correlation_matrix: The correlation matrix between the channels.

    Returns:
        The noise matrix.
    """
    if correlation_matrix is None:
        correlation_matrix=np.eye(len(shape))
    D,V = np.linalg.eigh(correlation_matrix, UPLO='U')
    corr_noise_factor = V * np.sqrt(np.diag(D)) * np.linalg.inv(V)  # 1 i
    return get_pseudo_noise(shape, corr_noise_factor)


def MRifft(k,dim):
    """
    Reconstruct individual coils' images and apply FFT scale factor
    iFFT scales data by 1/sqrt(N), so I need to scale back, as the noise covariance
    matrix is calculated in k-space (not Fourier transformed)

    Args:
        k: The k-space data

    Returns:
        o: The reconstructed images
    """
    tmp = np.fft.ifftshift(k)
    for idim in dim:
        tmp = np.fft.ifft(tmp, axis=idim)
    output = np.fft.ifftshift(tmp)
    return output


def MRfft(k,dim):
    """
    Reconstruct individual coils' images and apply FFT scale factor
    iFFT scales data by 1/sqrt(N), so I need to scale back, as the noise covariance
    matrix is calculated in k-space (not Fourier transformed)

    Args:
        k: The k-space data

    Returns:
        o: The reconstructed images
    """
    tmp = np.fft.ifftshift(k)
    for idim in dim:
        tmp = np.fft.fft(tmp, axis=idim)
    output = np.fft.ifftshift(tmp)
    return output


    
def writeResultsAsCmJSONOutput(results,filename,info=None):
    """_summary_

    Args:
        images (_type_): array of [type,name,imaginable,outputName]
        filename (_type_): oyutput filename
        info (_type_, optional): type, subType. Defaults to None.
    """
    OUTPUT={'version':'20220314',
    'language':'python',
    'type':'CMOutput',
    'subType':"default",
    'author':'eros.montin@gmail.com',
            }

        
    if info is not None:
        if info["type"] is not  None:
            OUTPUT["type"]=info['type']
        if info['subType'] is not None:
            OUTPUT["subType"]=info['subType']
    IMAGES=[]
    for r in results:
        if r['type']=='imaginable2D':
            if ((isinstance(r['imaginable'],im.Imaginable))):
                theim=r['imaginable']
            elif isinstance(r['imaginable'],sitk.Image):
                theim=im.Imaginable()
                theim.setImage(r['imaginable'])
            else:
                try:
                    if r['imaginable'].isImaginable():
                        theim=r['imaginable']
                    else:
                        return False
                except:
                    return False


            
            getSlice=theim.getAxialSlice
            orientation=2
            try:
                if r['orientation']==1:
                    getSlice=theim.getSagittalSlice
                    print('sagittal slices:)')
                elif r['orientation']==0:
                    getSlice=theim.getCoronalSlice
                    print('corona; slices:)')
                orientation=r['orientation']
                
            except:
                print('axial slices:)')
            
            S=theim.getImageSize()

            SL=[]
            for s in range(S[orientation]):
                o={"type":"double complex"}
                sl=getSlice(s)
                o['w'],o['h']=sl.GetSize()
                o['s']=sl.GetSpacing()
                o['o']=sl.GetOrigin()
                o['d']=sl.GetDirection()
                nda = sitk.GetArrayFromImage(sl) #YX
                s=nda.flatten('F').astype(complex)
                o['Vi']=s.imag.tolist()
                o['Vr']=s.real.tolist()
                SL.append(o)
            R={
                "slice":SL, # this should be slices... i know
                "imageName":r["name"]
            }
            
            try:
                R["imageOutputName"]=r["outputName"]
            except:
                R["imageOutputName"]=r["name"]
            
            IMAGES.append(R)
    
    OUTPUT["images"]=IMAGES

    if filename is None:
        return OUTPUT
    else:

        # create json object from dictionary
        js = json.dumps(OUTPUT)

        # open file for writing, "w" 
        f = open(filename,"w")

        # write json object to file
        f.write(js)

        # close file
        f.close()
        return True
    
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
	
import scipy
def savemat(fn,var):
    scipy.io.savemat(fn, {"var":var})


import matplotlib.pyplot as plt

def mrir_noise_bandwidth(noise):
    """Calculates the noise bandwidth.

    Args:
        noise: The noise data.
        display: Whether to display the noise power spectrum.

    Returns:
        The noise bandwidth.
    """

    dims = noise.shape
    power_spectrum = np.abs(np.fft.fft(noise, axis=0))**2
    power_spectrum_chan = np.mean(power_spectrum, axis=1)
    N=np.concatenate([power_spectrum_chan[:int(dims[0] / 4)],power_spectrum_chan[int(3*dims[0] / 4)::]])
    power_spectrum_norm = np.mean(N, axis=0)
    norm_power_spectrum_chan = power_spectrum_chan / np.repeat(power_spectrum_norm[None, :], dims[0], axis=0)
    noise_bandwidth_chan = np.sum(norm_power_spectrum_chan, axis=0) / dims[0]
    noise_bandwidth = np.mean(noise_bandwidth_chan)

    return noise_bandwidth