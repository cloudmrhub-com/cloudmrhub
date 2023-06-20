
import numpy as np
import cm2D
import cm


class cm2DReconRSS(cm2DReconW.cm2DRecon):
    """
    Python implementation of the cm2DReconRSS MATLAB class
   
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

if __name__ == "__main__":
    L=cm2DReconRSS()
    L.test()