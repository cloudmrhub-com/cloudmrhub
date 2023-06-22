% we are testing the python implementation of B1 recon and it's SNR

load d.mat
% N =N(1:2,1:3,1:4);
% S =S(1:2,1:3,1:4);

%B1reconstruction
%let's  reconstruct as Sum of square
%instantiate the reconstruct class
RE= cm2DReconRSS();
%instantiate the SNR calculation
F=cm2DSignalToNoiseRatioPseudoMultipleReplicas();
F.setNumberOfPseudoReplicas(20)
F.setReconstructor(RE)
F.setSignalKSpace(S)
F.setNoiseKSpace(N);
figure();
imagesc(abs(F.getOutput()),[]);
colorbar();
title('SNR RSS 20 replicas')
