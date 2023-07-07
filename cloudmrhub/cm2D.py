import numpy as np
try:
    import cm #dev
except:
    import cloudmrhub.cm as cm #runtime
    
import matplotlib.pyplot as plt
import scipy
from types import MethodType

class cm2DRecon(cm.cmOutput):
    """
    Python implementation of the cm2DRecon MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2023
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.


    :Attributes:
        SignalKSpace (cm.k2d): The signal k-space data.
        
        NoiseKSpace (cm.k2d): The noise k-space data.
                
        NoiseCovariance (np.ndarray): The noise covariance matrix.
        
        InverseNoiseCovariance (np.ndarray): The Inverse Noise Covariance matrix.
        
        SignalPrewhitened (np.ndarray): The prewhitened signal k-space data.
        
        HasSensitivity (bool): Whether the reconstruction contains sensitivity information.
        
        HasAcceleration (bool): Whether the reconstruction has accelerations.
   
    """
    def __init__(self):
        """
        Initializes the cm2DRecon object.
        """
        self.SignalKSpace = cm.k2d()        
        self.NoiseKSpace = cm.k2d()
        self.NoiseCovariance = np.array([])
        self.InverseNoiseCovariance = np.array([])
        self.SignalPrewhitened = cm.k2d()
        self.HasSensitivity = False
        self.HasAcceleration = False
        self.NoiseBandWidth = None

    def setSignalKSpace(self, signalKSpace):
        """
        Sets the signal k-space data.
        
        :param signalkspace: The signal k-space data
        :type: np.ndarray
        """
        
        self.SignalKSpace.set(signalKSpace)
        self.SignalPrewhitened.reset()

    def getSignalKSpace(self):
        """
        Gets the signal k-space data.
        
        Returns:
            _type_: nd.array(f,p,c)
        """
        return self.SignalKSpace.get()
    def getSignalKSpaceSize(self):
        """
        Gets the signal k-space size.
        
        Returns:
            NoiseKspace: nd.array(f,p,c)
        """
        return self.SignalKSpace.getSize()
    def getsignalNCoils(self):
        """
        Gets the number of coils in the signal k-space data.
        
        Returns:
            int: ncoils
        """
        return self.SignalKSpace.getNCoils()
    def setNoiseKSpace(self, noiseKSpace):
        """ Sets the noise k-space data.

        :param noiseKspace: The noise k-space data
        :type: np.ndarray(f,p,c)
        """
        self.NoiseKSpace.set(noiseKSpace)
    def getNoiseKSpace(self):
        """
        Gets the noise k-space data.
        
        Returns:
            NoiseKspace: nd.array(f,p,c)
        """
        return self.NoiseKSpace.get()
    def getNoiseKSpaceSize(self):
        """
        Gets the noise k-space size.
        
        Returns:
            NoiseKspace: nd.array(f,p,c)
        """
        return self.NoiseKSpace.getSize()
    def getNoiseNCoils(self):
        """
        Gets the number of coils in the noise k-space data.
        
        Returns:
            int: ncoils
        """
        return self.NoiseKSpace.getNCoils()
    
    def getNoiseBandWidth(self):
            if(self.NoiseBandWidth is None):
                noise_bandwidth = cm.mrir_noise_bandwidth(self.getNoiseKSpace());
                self.NoiseBandWidth=noise_bandwidth
            try:
                return  self.NoiseBandWidth
            except:
                print("---problem in the noisebadwidth----")
                self.appendLog("---problem in the noisebadwidth----","warning")
                return  1.0
            
    def getInverseNoiseCovariancePrewhitened(self):
        """when prewhitened the correlation is an eye

        """
        return np.eye(self.getsignalNCoils())
    def setNoiseCovariance(self, noiseCovariance):
        """
        Sets the noise covariance matrix.

        :param noiseCovariance: The noise covariance matrix.
        :type: mp.ndarray(c,c)

        """
        self.NoiseCovariance = noiseCovariance
        self.InverseNoiseCovariance = np.linalg.inv(noiseCovariance)
    def getNoiseCovariance(self):
        """
        Return the covariance matrix
        
        Returns:
            np.ndarray(c,c): Covariance Matrix
        """
        if not self.NoiseCovariance.any():
            self.NoiseCovariance = cm.calculate_covariance_matrix(self.getNoiseKSpace(),self.getNoiseBandWidth())

        return self.NoiseCovariance
    def getInverseNoiseCovariance(self):
        """
        Gets the inverse noise covariance matrix.

        Returns:
            np.ndarray(c,c): inverse of the covariance matrix
        """
        return self.InverseNoiseCovariance

    def getPrewhitenedSignal(self):
        """
        Gets the prewhitened signal.

        Returns:
            np.ndarray(f,p,c): prewhitened signal
        """
        if self.SignalPrewhitened.isEmpty():
            self.SignalPrewhitened.set(cm.prewhiteningSignal(self.getSignalKSpace(), self.getNoiseCovariance()))
        return self.SignalPrewhitened.get()
    def setPrewhitenedSignal(self, prewhitenedSignal):
        self.SignalPrewhitened.set(prewhitenedSignal)
    def plotImageAfterTest(self,IM,tit):
        # Create the figure and subplots
        fig, axarr = plt.subplots(2, 1)
        axarr[0].imshow(IM)
        axarr[0].set_title(tit)
        axarr[1].imshow(abs(IM))
        axarr[1].set_title(tit + ' abs')

        ha = plt.axes([0, 0, 1, 1], frameon=False, visible=False,
              xlim=[0, 1], ylim=[0, 1], aspect='equal')
        # Add text to the axes object
        text = ha.text(0.5, 0.98, f'{IM.shape[0]}x{IM.shape[1]}',
                    transform=ha.transAxes, horizontalalignment='center',
                    verticalalignment='top')
        plt.show()
    def test(self):
        TEST = self.testrecon()
        self.plotImageAfterTest(TEST.getOutput(), "recon")
        return TEST
    def get2DKSIFFT(self,K=None):
        if K is None:
            K=self.getPrewhitenedSignal()
        SC=np.sqrt(np.prod(np.array(K.shape[0:2])))
        return cm.MRifft(K,[0,1])*SC
    def getOutput(self):
        pass
    def resetAfterSignal(self):
        self.SignalPrewhitened.reset()
    def resetAfterNoise(self):
        if not isinstance(self.NoiseBandWidth,str):
            self.NoiseBandWidth=None
            print('reset the NBW')
        self.NoiseCovariance = np.array([])
        self.InverseNoiseCovariance = np.array([])
        self.SignalPrewhitened.reset()

        

class cm2DReconRSS(cm2DRecon):
    """
    Python implementation of the cm2DReconRSS MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2023
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    
    def __init__(self):
        """
        Initializes the RSS reconstruction.
        
        """
        super().__init__()
        self.HasAcceleration= False
        self.HasSensitivity=False
    def getOutput(self):
        img_matrix = self.get2DKSIFFT()
        im = np.sqrt(np.sum(np.abs(img_matrix)**2,axis=-1))
        return im

    @staticmethod
    def testrecon():
        TEST = cm2DReconRSS()
        [K, N, nc] = cm.getMRoptimumTestData()
        TEST.setSignalKSpace(K)
        TEST.setNoiseCovariance(nc)
        return TEST

class cm2DReconSS(cm2DReconRSS):
    """
    Python implementation of the cm2DReconSS MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2023
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    
    def getOutput(self):
        img_matrix = self.get2DKSIFFT()
        im = np.sum(img_matrix**2,axis=-1)
        return im

    @staticmethod
    def testrecon():
        TEST = cm2DReconSS()
        [K, N, nc] = cm.getMRoptimumTestData()
        TEST.setSignalKSpace(K)
        TEST.setNoiseCovariance(nc)
        return TEST
class cm2DReconRSSunAbs(cm2DReconSS):
    """
    Python implementation of the cm2DReconRSSunAbs MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2023
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    def getOutput(self):
        im = np.sqrt(super().getOutput())
        return im

    @staticmethod
    def testrecon():
        TEST = cm2DReconRSSunAbs()
        [K, N, nc] = cm.getMRoptimumTestData()
        TEST.setSignalKSpace(K)
        TEST.setNoiseCovariance(nc)
        return TEST


class cm2DKellmanRSS(cm2DReconRSS):
    def __init__(self):
        super().__init__()
    def getOutput(self):
        nf,nph =self.getSignalKSpaceSize()
        img_matrix = self.get2DKSIFFT()
        snr = np.zeros((nf,nph))
        for irow in range(nf):
            for icol in range(nph):
                        B=np.expand_dims(img_matrix[irow,icol],axis=-1)
                        A=B.conj().T
                        snr[irow,icol] = np.abs(np.sqrt(2*(A @ B)))
        return snr
    


class cm2DReconWithSensitivity(cm2DRecon):
    """
    Python implementation of the cm2DReconWithSensitivity MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2003
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    
    def __init__(self):
        """
        Initializes the reconstruction.
        
        """
        super().__init__()
        self.HasAcceleration= False
        self.HasSensitivity=True
        self.CoilSensitivityMatrixCalculationMethod='inner'
        self.CoilSensitivityMatrix=cm.k2d()
        self.CoilSensitivityMatrixSource=cm.k2d()
        self.CoilSensitivityMatrixSourcePrewhitened=cm.k2d()
        self.MaskCoilSensitivityMatrix='ref'

    
    def setCoilSensitivityMatrix(self, S):
        self.CoilSensitivityMatrix.set(S)

    def resetCoilSensitivityMatrix(self):
        self.CoilSensitivityMatrix.reset()

    def __selectCoilSensitivityMapMethod__(self):
        s=None
        if ((self.getCoilSensitivityMatrixCalculationMethod() == 'simplesense') or (self.getCoilSensitivityMatrixCalculationMethod() =='inner')):
                s=self.getCoilSensitivityMatrixSimpleSense
        return s

    def calculateCoilSensitivityMatrix(self):
        if self.CoilSensitivityMatrix.isEmpty():
            CSM=self.__selectCoilSensitivityMapMethod__()
            s=CSM()
        else:
            s=self.getCoilSensitivityMatrix()
        return s

    def getCoilSensitivityMatrix(self):
        if self.CoilSensitivityMatrix.isEmpty():
            coilsens_set = self.calculateCoilSensitivityMatrix()
            self.setCoilSensitivityMatrix(coilsens_set)
        return self.CoilSensitivityMatrix.get()

    def setCoilSensitivityMatrixSourcePrewhitened(self, x):
        self.CoilSensitivityMatrixSourcePrewhitened.set(x)

    def getCoilSensitivityMatrixSource(self):
        return self.CoilSensitivityMatrixSource.get()
    
    def getCoilSensitivityMatrixSourcePrewhitened(self):
        if self.CoilSensitivityMatrixSourcePrewhitened.isEmpty():
            S = self.getCoilSensitivityMatrixSource()
            Rn = self.getNoiseCovariance()
            pw_S = cm.prewhiteningSignal(S, Rn)
            self.setCoilSensitivityMatrixSourcePrewhitened(pw_S)
        return pw_S

    def setCoilSensitivityMatrixSource(self, IM):
        self.CoilSensitivityMatrixSource.set(IM)

    def setMaskCoilSensitivityMatrix(self, x):
        self.MaskCoilSensitivityMatrix = x

    def getMaskCoilSensitivityMatrix(self):
        return self.MaskCoilSensitivityMatrix

    def setCoilSensitivityMatrixCalculationMethod(self, x):
        self.CoilSensitivityMatrixCalculationMethod = x

    def getCoilSensitivityMatrixCalculationMethod(self):
        return self.CoilSensitivityMatrixCalculationMethod

    def getCoilSensitivityMatrixSimpleSense(self):
        # MASK 
        if self.CoilSensitivityMatrix.isEmpty():       
            if self.CoilSensitivityMatrixSource.isEmpty():
                s=self.getSignalKSpace()
            else:
                s=self.getCoilSensitivityMatrixSource()

            self.setCoilSensitivityMatrix(cm.prewhiteningSignal(cm.calculate_simple_sense_sensitivitymaps(s,self.MaskCoilSensitivityMatrix), self.getNoiseCovariance() ))
        return self.CoilSensitivityMatrix.get()
    
 
    def resetAfterSignal(self):
        super().resetAfterSignal()
        self.CoilSensitivityMatrix.reset()        

        
    def resetAfterNoise(self):
        super().resetAfterNoise()
        self.CoilSensitivityMatrix.reset()

    




class cm2DReconB1(cm2DReconWithSensitivity):
    """
    Python implementation of the cm2DReconB1 MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2003
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    
    def __init__(self):
        """
        Initializes the RSS reconstruction.
        
        """
        super().__init__()
        self.HasAcceleration= False
        self.HasSensitivity=True
        
    def getOutput(self):
        img_matrix = self.get2DKSIFFT()
        pw_sensmap=self.getCoilSensitivityMatrix()
        invRn=self.getInverseNoiseCovariancePrewhitened()
        nf,nph =self.getSignalKSpaceSize()
        im = np.zeros((nf,nph),dtype=pw_sensmap.dtype)
        for irow in range(nf):
            for icol in range(nph):
                        s_matrix=pw_sensmap[irow,icol,:]
                        if s_matrix.sum() !=0:
                            im[irow,icol] = s_matrix.conj().T @ invRn @ img_matrix[irow,icol,:]
        return im
    

    @staticmethod
    def testrecon():
        TEST = cm2DReconB1()
        [K, N, nc] = cm.getMRoptimumTestData()
        TEST.setSignalKSpace(K)
        TEST.setNoiseCovariance(nc)
        return TEST

class cm2DKellmanB1(cm2DReconB1):
    """
    Python implementation of the cm2DReconB1 MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2003
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    
    def getOutput(self):
        img_matrix = self.get2DKSIFFT()
        pw_sensmap=self.getCoilSensitivityMatrix()
        invRn=self.getInverseNoiseCovariancePrewhitened()
        nf,nph =self.getSignalKSpaceSize()
        im = np.zeros((nf,nph))
        SR=np.sqrt(2.0)
        for irow in range(nf):
            for icol in range(nph):
                        s_matrix=pw_sensmap[irow,icol,:]
                        if s_matrix.sum() !=0:
                            S=s_matrix
                            ST=s_matrix.conj().T
                            I=img_matrix[irow,icol,:]
                            num=np.dot(SR,np.abs(ST @ invRn @ I))
                            den=np.sqrt(np.abs(ST @ invRn @ S))
                            im[irow,icol] = np.divide(num,den)
        return im

class cm2DReconWithSensitivityAutocalibrated(cm2DReconWithSensitivity):
    """
    Python implementation of the m2DReconWithSensitivityAutocalibrated MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2003
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    
    def __init__(self):
        """
        Initializes the reconstruction.
        
        """
        super().__init__()
        self.HasAcceleration= True
        self.HasSensitivity=True
        self.AutocalibrationF =None
        self.AutocalibrationP =None
        self.AccelerationF =0
        self.AccelerationP=0

    def setAutocalibrationLines(self,ACL):
        self.AutocalibrationF=ACL[0]
        self.AutocalibrationP=ACL[1]

    def getCoilSensitivityMatrixSimpleSenseACL(self):
        # MASK 
        if self.CoilSensitivityMatrix.isEmpty():       
            if self.CoilSensitivityMatrixSource.isEmpty():
                s=self.getSignalKSpace()
            else:
                s=self.getCoilSensitivityMatrixSource()
            self.setCoilSensitivityMatrix(cm.prewhiteningSignal(cm.calculate_simple_sense_sensitivitymaps_acl(s,self.AutocalibrationF,self.AutocalibrationP,self.MaskCoilSensitivityMatrix), self.getNoiseCovariance() ))
        return self.CoilSensitivityMatrix.get()


    def __selectCoilSensitivityMapMethod__(self):
        s=None
        s=super().__selectCoilSensitivityMapMethod__()
        if s is None:
            if ((self.getCoilSensitivityMatrixCalculationMethod() == 'simplesenseACL') or (self.getCoilSensitivityMatrixCalculationMethod() =='innerACL')):
                s=self.getCoilSensitivityMatrixSimpleSenseACL
        return s



class cm2DReconSense(cm2DReconWithSensitivityAutocalibrated):
    """_summary_
    class to reconstruct as SENSE the kspace.
    AutocalibrationsP are set to 0
    AutocalibrationF are set after setting the songalkspace
    to KS.shape[0]

    """
    def __init__(self):
        super().__init__()
        delattr(self,'AutocalibrationF')
        delattr(self,'AutocalibrationP')

    def setAutocalibrationLines(self,ACL):
        [nf,_]=self.getSignalKSpaceSize()
        self.AutocalibrationF=nf
        self.AutocalibrationP=ACL[1]

    def undersamplerZeroPadded(self,k):
        return cm.undersample2DDataSENSE(k,self.AccelerationF,self.AccelerationP)

    def setSignalKSpace(self, signalKSpace):
        super().setSignalKSpace(signalKSpace)
        self.AutocalibrationF=signalKSpace.shape[0]
        
        
        
    def getOutput(self):
        R1=self.AccelerationF
        R2=self.AccelerationP
        Rtot = R1 * R2
        nc=self.getsignalNCoils()
        nf,nph = self.getSignalKSpaceSize()
        invRn=self.getInverseNoiseCovariancePrewhitened()
        pw_signalrawdata=self.getPrewhitenedSignal()     
        k_temp = self.undersamplerZeroPadded(pw_signalrawdata)
        img_matrix = self.get2DKSIFFT(k_temp) * np.sqrt(Rtot)
        imfold=cm.shrinktoaliasedmatrix_2d(img_matrix,R1,R2)
        
        pw_sensmap=self.getCoilSensitivityMatrix()
        MRimage = np.zeros((nf, nph),dtype=complex)
        for irow in range(nf // R1):
            for icol in range(nph//R2):
                r1=np.arange(irow,nf,nf//R1)
                r2=np.arange(icol,nph,nph//R2)
                current_R1 = len(r1)
                current_R2 = len(r2)
                current_Rtot = current_R1 * current_R2
                s=np.zeros((current_R1,current_R2,nc),dtype=complex)
                for i,_x in enumerate(r1):
                    for j,_y in enumerate(r2):
                        s[i,j,:] = pw_sensmap[_x,_y, :]
                s = s.reshape((current_Rtot, nc),order='F')
                s = np.transpose(s,[1,0])
                s[np.isnan(s)] = 0 + 1j*0
                u = np.linalg.pinv(s.conj().T @ invRn @ s) @ (s.conj().T @ invRn)
                u[np.isnan(u)] = 0 + 1j*0
                U=np.reshape(u @ imfold[irow,icol],(current_R1,current_R2),order='F')
                for i,_x in enumerate(r1):
                    for j,_y in enumerate(r2):
                        MRimage[_x,_y] = U[i,j]

        if (cm.needs_regridding(self.getSignalKSpace(),R1,R2)):
            MRimage=cm.resizeIM2D(MRimage,self.getSignalKSpaceSize());
        return MRimage
                

class cm2DKellmanSense(cm2DReconSense):
    def getOutput(self):
        R1=self.AccelerationF
        R2=self.AccelerationP
        Rtot = R1 * R2
        nc=self.getsignalNCoils()
        nf,nph = self.getSignalKSpaceSize()
        invRn=self.getInverseNoiseCovariancePrewhitened()
        pw_signalrawdata=self.getPrewhitenedSignal()     
        k_temp = self.undersamplerZeroPadded(pw_signalrawdata)
        img_matrix = self.get2DKSIFFT(k_temp) * np.sqrt(Rtot)
        imfold=cm.shrinktoaliasedmatrix_2d(img_matrix,R1,R2)
        
        pw_sensmap=self.getCoilSensitivityMatrix()
        Rn = np.eye(nc); #it's prewhitened
        MRimage = np.zeros((nf, nph),dtype=complex)
        _c=np.sqrt(2.0)
        for irow in range(nf // R1):
            for icol in range(nph//R2):
                r1=np.arange(irow,nf,nf//R1)
                r2=np.arange(icol,nph,nph//R2)
                current_R1 = len(r1)
                current_R2 = len(r2)
                current_Rtot = current_R1 * current_R2
                s=np.zeros((current_R1,current_R2,nc),dtype=complex)
                for i,_x in enumerate(r1):
                    for j,_y in enumerate(r2):
                        s[i,j,:] = pw_sensmap[_x,_y, :]
                s = s.reshape((current_Rtot, nc),order='F')
                s = np.transpose(s,[1,0])
                s[np.isnan(s)] = 0 + 1j*0
                u = np.linalg.pinv(s.conj().T @ invRn @ s) @ (s.conj().T @ invRn)
                u[np.isnan(u)] = 0 + 1j*0
                U=np.reshape(u @ imfold[irow,icol]/np.diag(np.sqrt(u @ Rn @ u.conj().T)),(current_R1,current_R2),order='F')
                for i,_x in enumerate(r1):
                    for j,_y in enumerate(r2):
                        MRimage[_x,_y] = _c*U[i,j]

        if (cm.needs_regridding(self.getSignalKSpace(),R1,R2)):
            MRimage=cm.resizeIM2D(MRimage,self.getSignalKSpaceSize())
        return MRimage


class cm2DGfactorSense(cm2DKellmanSense):
    def getOutput(self):
        R1=self.AccelerationF
        R2=self.AccelerationP
        Rtot = R1 * R2
        nc=self.getsignalNCoils()
        nf,nph = self.getSignalKSpaceSize()
        invRn=self.getInverseNoiseCovariancePrewhitened()
        pw_signalrawdata=self.getPrewhitenedSignal()     
        k_temp = self.undersamplerZeroPadded(pw_signalrawdata)
        pw_sensmap=self.getCoilSensitivityMatrix()
        MRimage = np.zeros((nf, nph),dtype=complex)
        _c=np.sqrt(2.0)
        for irow in range(nf // R1):
            for icol in range(nph//R2):
                r1=np.arange(irow,nf,nf//R1)
                r2=np.arange(icol,nph,nph//R2)
                current_R1 = len(r1)
                current_R2 = len(r2)
                current_Rtot = current_R1 * current_R2
                s=np.zeros((current_R1,current_R2,nc),dtype=complex)
                for i,_x in enumerate(r1):
                    for j,_y in enumerate(r2):
                        s[i,j,:] = pw_sensmap[_x,_y, :]
                s = s.reshape((current_Rtot, nc),order='F')
                s = np.transpose(s,[1,0])
                s[np.isnan(s)] = 0 + 1j*0
                u1=s.conj().T @ invRn @ s
                u = np.diag(np.linalg.pinv(u1))*np.diag(u1)
                u[np.isnan(u)] = 0 + 1j*0
                U=np.reshape(u,(current_R1,current_R2),order='F')
                for i,_x in enumerate(r1):
                    for j,_y in enumerate(r2):
                        MRimage[_x,_y] = _c*U[i,j]

        if (cm.needs_regridding(self.getSignalKSpace(),R1,R2)):
            MRimage=cm.resizeIM2D(MRimage,self.getSignalKSpaceSize())
        return MRimage
    
class cm2DReconmSense(cm2DReconSense):
    
    def __init__(self):
        super().__init__()
        self.AutocalibrationP=0


class cm2DGfactormSense(cm2DGfactorSense):
    def __init__(self):
        super().__init__()
        self.AutocalibrationP=0


class cm2DKellmanmSense(cm2DKellmanSense):
    def __init__(self):
        super().__init__()
        self.AutocalibrationP=0


class cm2DReconGrappa(cm2DReconWithSensitivityAutocalibrated):
    """Grappa REconstruction
    Args:
        cm2DReconWithSensitivityAutocalibrated (_type_): _description_
    """
    def __init__(self):
        super().__init__()
        self.HasSensitivity=False
        self.GrappaKernel=[3,2]
        self.PrewhitenedSignalKspaceACL=cm.k2d()
        # delattr(self,'AccelerationF')
        self.reconstructor=cm2DReconRSS()

    
    def getPrewhitenedSignalACL(self):
        """
        Gets the prewhitened signalACL.

        Returns:
            np.ndarray(f,p,c): prewhitened signal
        """
        if self.PrewhitenedSignalKspaceACL.isEmpty():
             
            if self.SignalPrewhitened.isEmpty():
                S=cm.prewhiteningSignal(self.getSignalKSpace(), self.getNoiseCovariance())
                self.SignalPrewhitened.set(S)
            else:
                S=self.getPrewhitenedSignal()
            ACL = cm.getAutocalibrationsLines2DKSpace(S,self.AutocalibrationF,self.AutocalibrationP)
            self.PrewhitenedSignalKspaceACL.set(ACL)
            return ACL
        else:
            return self.PrewhitenedSignalKspaceACL.get()



    def setGrappaKernel(self,GK):
        self.GrappaKernel=GK
        
    
    def getR(this):
        SS = this.getSignalKSpaceSize()
        R = this.AccelerationP
        np = SS[1]

        if np % R != 0:
            tempR = R
            while np % tempR != 0:
                tempR -= 1
            R = tempR
            ss = f'Acceleration R reduced to {R} ' \
                f', so the number of lines can be exactly divided by R'
            print(ss)
        return R


    def getOutput(self):
        grappa_kernel=self.GrappaKernel
        data_acs=self.getPrewhitenedSignalACL()
        pw_signalrawdata=self.getPrewhitenedSignal(); 
        K=cm.getGRAPPAKspace(pw_signalrawdata,data_acs,grappa_kernel)
        # R=self.recon()
        R=self.reconstructor
        R.setPrewhitenedSignal(K)
        return R.getOutput()



##MR PMR
class cm2DSignalToNoiseRatio(cm.cmOutput):
    def __init__(self, message=None):
        super().__init__(message)
        """
        the class expects a 3D matrix composed by a tile of 2D numpy images
        """
        self.appendLog("SNR Calculation started instantiated", "start")
        self.SNR = None
        self.Type = "MR"
        self.SubType = ""

def resetASPMR(self):
    self.SignalPrewhitened.reset()
def resetANPMR(self):
    if not isinstance(self.NoiseBandWidth,str):
        self.NoiseBandWidth=None
        print('reset the NBW')
    self.NoiseCovariance = np.array([])
    self.InverseNoiseCovariance = np.array([])
    self.SignalPrewhitened.reset()

        
    

class cm2DSignalToNoiseRatioMultipleReplicas(cm2DSignalToNoiseRatio):
    def __init__(self,x=None,message=None):
        """
        the class expects a 3D matrix composed by a tile of 2D images
        """
        super().__init__(message)
        self.imageArray=cm.i3d()
        self.Mean=None
        self.STD=None
        self.Max=None
        self.Type = "MR"
        self.reconstructor=None
        self.referenceImage=cm.i2d()

        
        if x is not None:
            self.add2DStackOfImages(x)
    def reset(self):
        self.Max=None
        self.Mean=None
        self.STD=None
        self.SNR=None
    def getReferenceImage(self):
        return self.referenceImage.get()
    def setReferenceImage(self,x):
        self.referenceImage.set(x)
    def add2DImage(self, x):
        self.reset()
        """
        add a 2D image to the class
        """
        if len(x.shape)==2:
            x=np.expand_dims(x,axis=-1)
        if self.imageArray.isEmpty():
            self.setImageArray(x)
        else:
            self.setImageArray(np.concatenate((self.getImageArray(), x), axis=-1))

    def add2DStackOfImages(self, x):
        """
        add a stack of 2D images to the class
        """
        self.add2DImage(x)
        # for t in range(x.shape[2]):
            # self.add2DImage(x[:, :, t])
    def add2DKspace(self,signal):        
        self.reconstructor.setSignalKSpace(signal)
        self.add2DImage(self.reconstructor.getOutput())
    
    def setReconstructor(self,recon):
        """We dont' want to recalculate the Sensitivity map at every replica

        Args:
            recon (cm2dRecon): REconstructor class
        """
        self.reconstructor=recon
        self.reconstructor.resetAfterSignal=MethodType(resetASPMR,self.reconstructor)
        self.reconstructor.resetAfterNoise=MethodType(resetANPMR,self.reconstructor)

    def getReconstructor(self):
        return self.reconstructor
       
       
    def resetAfterNoise(self,keepsensitivity=False):
        super().resetAfterNoise()
        if not keepsensitivity:
            self.setCoilSensitivityMatrix.reset()


    def getImageArray(self):
        """
        return the image array
        """
        return self.imageArray.get()
    def setImageArray(self,x):
        """
        return the image array
        """
        self.imageArray.set(x)

    def getOutput(self):
        """
        calculate the mean and standard deviation of the image array
        """
        if self.SNR is None:
            self.SNR = np.divide(self.getImageArrayMean(), self.getImageArraySTD())
        return self.SNR
    def getImageArrayMean(self):
        """
        return the mean of the image array
        """
        if self.Mean is None:
            self.Mean=np.nanmean(np.abs(self.getImageArray()), axis=-1)
        return self.Mean
    def getImageArraySTD(self):
        """
        return the standard deviation of the image array
        """
        if self.STD is None:
            self.STD=np.nanstd(np.abs(self.getImageArray()), axis=-1,ddof=1) #matlab 0
        return self.STD
    def getImageArrayMax(self):
        """
        return the standard deviation of the image array
        """
        if self.Max is None:
            self.Max=np.nanmax(np.abs(self.getImageArray()), axis=-1) #matlab 0
        return self.Max
 
    def plotImageArray(self, p=0.5):
        """
        plot the image array
        """
        im = self.getImageArray()
        for t in range(im.shape[-1]):
            plt.subplot(121)
            plt.imshow(im[:, :, t])
            plt.colorbar()
            plt.title("Replicas number: " + str(t+1))
            if t>0:
                plt.subplot(122)
                plt.imshow(im[:, :, t]-im[:, :, t-1])
                plt.title(f'differerence {t+1} - {t}')
                plt.colorbar()
            plt.pause(interval=1)
            plt.show()

class cm2DSignalToNoiseRatioPseudoMultipleReplicas(cm2DSignalToNoiseRatioMultipleReplicas):
    def __init__(self, x=None, message=None):
        super().__init__(x, message)
        self.numberOfReplicas=20
        self.D=None
    def createPseudoReplica(self,S,corr_noise_factor):
        sh=self.reconstructor.getSignalKSpaceSize()
        N= cm.get_pseudo_noise(msize=[*sh, self.reconstructor.getsignalNCoils()],corr_noise_factor=corr_noise_factor)
        self.add2DKspace(S+N)
    def getSNRDenumerator(self):
        if self.D is None:
            D=np.nanstd(np.abs(self.getImageArray())+np.max(np.abs(self.getImageArrayMax())), axis=-1,ddof=1)
            D[D<=np.finfo(np.float64).eps]=1
            self.D=D
        return self.D
    def reset(self):
        super().reset()
        self.D=None
    
    def getOutput(self):
        # set the reference image
        # self.setReferenceImage(self.reconstructor.getOutput())
        corr_noise_factor=cm.get_correlation_factor(correlation_matrix=self.reconstructor.getNoiseCovariance())
        S=self.reconstructor.getSignalKSpace()

        for a in range(self.numberOfReplicas):
            #add in the queue
            self.createPseudoReplica(S,corr_noise_factor)
        if self.referenceImage.isEmpty():
            self.setReferenceImage(self.reconstructor.getOutput())
        D=self.getSNRDenumerator()
        SNR=np.divide(self.getReferenceImage(),D)
        return SNR

class cm2DSignalToNoiseRatioPseudoMultipleReplicasWen(cm2DSignalToNoiseRatioPseudoMultipleReplicas):
    def __init__(self, x=None, message=None):
        super().__init__(x, message)
        self.Type='CR'
        self.boxSize=2

    def getSNRDenumerator(self):
        if self.referenceImage.isEmpty():
            self.setReferenceImage(self.reconstructor.getOutput())
        r=self.getReferenceImage()-self.getImageArrayMean()
        return cm.get_wien_noise_image(r,self.boxSize)
