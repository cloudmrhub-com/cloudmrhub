% we are tteesting the python implementation of B1 recon and it's SR

load d.mat
x=[1:18]+(96/4);
y=[1:16]+(96/4);
N =N(x,y,:);
S =S(x,y,:);




R1=1; %acceleration on frequency direction
R2=2; % acceleration of phase direction
ACL=16; %number of autocalibration lines for sense reconstruction
ACLG=[ACL, ACL]; %number of autocalibration lines for grappa reconstruction (frequency, phase)
GK=[7 8]; %grappa kernel
SOSR=1; %sumofsquare 

RE2= cm2DReconRSS();
RE2.setSignalKSpace(S)
RE2.setNoiseKSpace(N)


%GRAPPA
TITLE='GRAPPA ';
RE= cm2DReconGRAPPA();
%mimik a grappa acquisition from a fully sample
k_temp=RE.mimicmGrappaDataFromFullysampledZeroPadded(S,R1,R2,ACLG(1),ACLG(2));
RE.setSignalKSpace(k_temp)
RE.setNoiseKSpace(N)
RE.setAcceleration(R2)
RE.setAutocalibration(ACLG)
RE.setGrappaKernel(GK);
RE.setSosRecon(SOSR);
figure();
subplot(121);
imagesc(abs(RE.getOutput()));
colorbar();title(['Reconstruction ' TITLE])

subplot(122);
imagesc(abs(RE2.getOutput()));
colorbar();title(['Reconstruction RSS'])


