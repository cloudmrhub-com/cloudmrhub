% we are testing the python implementation of B1 recon and it's SR

load d.mat
N =N(1:8,1:6,1:2);
S =S(1:8,1:6,1:2);




R1=1; %acceleration on frequency direction
R2=2; % acceleration of phase direction
ACL=NaN; %number of autocalibration lines for sense reconstruction
ACLG=[ACL ACL]; %number of autocalibration lines for grappa reconstruction (frequency, phase)
GK=[7 8]; %grappa kernel
SOSR=1; %sumofsquare 

% 
% %GRAPPA
% TITLE='GRAPPA '  
% RE= cm2DReconGRAPPA();
% %mimik a grappa acquisition from a fully sample
% k_temp=RE.mimicmGrappaDataFromFullysampledZeroPadded(S,R1,R2,ACL,ACL);
% RE.setSignalKSpace(k_temp)
% RE.setNoiseKSpace(N)
% RE.setAcceleration(R2)
% RE.setAutocalibration(ACLG)
% RE.setGrappaKernel(GK);
% RE.setSosRecon(SOSR);
% figure();
% imshow(abs(RE.getOutput()),[]);
% colorbar();title(['Reconstruction ' TITLE])



%SENSE
TITLE='SENSE ' ;
RE=cm2DGFactorSENSE()
ACCK=RE.mimicmSenseDataFromFullysampledZeroPadded(S,R1,R2,ACL);
RE.setSignalKSpace(ACCK)
RE.setNoiseKSpace(N)
RE.setCoilSensitivityMatrixSource(S)
RE.setCoilSensitivityMatrixCalculationMethod('simplesense');
%don't need to set the sensitivity source since it will be the signal for simple
RE.setAccelerationFrequency(R1)
RE.setAccelerationPhase(R2);
RE.setAutocalibration(ACL);
OUT.sense=RE.getOutput();
figure();
imagesc(abs(RE.getOutput()));
colorbar();title(['Reconstruction ' TITLE])