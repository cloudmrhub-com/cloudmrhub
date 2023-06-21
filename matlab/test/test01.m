% we are testing the python implementation of B1 recon and it's SNR

load d.mat
% N =N(1:2,1:3,1:4);
% S =S(1:2,1:3,1:4);

%B1reconstruction
TITLE='B1 '  
RE= cm2DReconB1();
RE.setCoilSensitivityMatrixSource(S)
RE.setCoilSensitivityMatrixCalculationMethod('simplesense');
RE.setSignalKSpace(S)
RE.setNoiseKSpace(N);


imagesc(abs(RE.getOutput()))