
import numpy as np
import cm2D
import cm


class cm2DReconWithSensitivity(cm2D.cm2DRecon):
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
        Initializes the RSS reconstruction.
        
        """
        super().__init__()
        self.HasAcceleration= False
        self.HasSensitivity=True
    

    def calculateCoilSensitivityMatrixFullySampled(self):
        if not isempty(self.CoilSensitivityMatrixSourcePrewhitened):
            this.logIt(['a Source Coil sensitivity map has been correctly set so i can calculate the senstivity map'], 'ok')
        switch(lower(self.CoilSensitivityMatrixCalculationMethod)):
            case {'espirit'}:
                this.logIt(['sensitivity map calculated as espirit'], '?')
                coilsens_set = this.getCoilSensitivityMatrixForFullySampledKSpaceEspirit()
            case {'simplesense', 'internal reference', 'inner'}:
                this.logIt(['sensitivity map calculated as simplesense'], '?')
                coilsens_set = this.getCoilSensitivityMatrixForFullySampledKSpaceSimpleSense()
            case 'adaptive':
                this.logIt(['sensitivity map calculated as aspative'], '?')
                coilsens_set = self.getCoilSensitivityMatrixForFullySampledKSpaceAdapt2D()
            case 'bodycoil':
                this.logIt(['sensitivity map calculated as aspative'], '?')
                coilsens_set = self.getCoilSensitivityMatrixForFullySampledKSpaceBodyCoil()
        return coilsens_set

    def setCoilSensitivityMatrix(self, S):
        self.CoilSensitivityMatrix = S

    def resetCoilSensitivityMatrix(self):
        self.setSensitivityMatrix([])

    def getCoilSensitivityMatrix(self):
        if isempty(self.CoilSensitivityMatrix):
            coilsens_set = self.calculateCoilSensitivityMatrix()
            if self.getCoilSensitivityMatrixSourceSmooth():
                for aa in range(size(coilsens_set, 3)):
                    coilsens_set[:, :, aa] = medfilt2(real(coilsens_set[:, :, aa]), (3, 3), 'symmetric') + 1j * medfilt2(
                        imag(coilsens_set[:, :, aa]), (3, 3), 'symmetric')
            self.CoilSensitivityMatrix = coilsens_set
        return self.CoilSensitivityMatrix

    def setCoilSensitivityMatrixSourcePrewhitened(self, x):
        self.CoilSensitivityMatrixSourcePrewhitened = x
        if isempty(self.getCoilSensitivityMatrixSourceNCoils()):
            self.setCoilSensitivityMatrixSourceNCoils(size(x, 3))

    def getCoilSensitivityMatrixSourcePrewhitened(self):
        pw_S = self.getCoilSensitivityMatrixSourcePrewhitened()
        if isempty(pw_S):
            S = self.getCoilSensitivityMatrixSource()
            Rn = self.getNoiseCovariance()
            if not isempty(Rn):
                pw_S = self.prewhiteningSignal(S, Rn)
                self.setCoilSensitivityMatrixSourcePrewhitened(pw_S)
        return pw_S

    def setCoilSensitivityMatrixSource(self, IM):
        self.CoilSensitivityMatrixSource = IM
        self.setCoilSensitivityMatrixSourceNCoils(size(IM, 3))

    def setMaskCoilSensitivityMatrix(self, x):
        self.MaskCoilSensitivityMatrix = x

    def getMaskCoilSensitivityMatrix(self):
        return self.MaskCoilSensitivityMatrix

    def setCoilSensitivityMatrixCalculationMethod(self, x):
        self.CoilSensitivityMatrixCalculationMethod = x
