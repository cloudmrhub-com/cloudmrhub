%let's reconstruct the phantom data in our release:
% - read the signal data
% - read the noise data
% - slect or loop over the slices:
% - customize the reconstruction parameter
% reconstruct the image

    
%the file that we are going to use    
NOISEFILE='../Data/AC1SLICE_1P5MM_NOISE.dat';
SIGNALFILE='../Data/AC1SLICE_1P5MM.dat';
 
%instantiate the readingFilter
N=readRawDataFilter('noise'); %we need to inform the class that it's going tor ead a noise data so it' won't remove the fewquency oversample
N.setFilename(NOISEFILE); %set the filename
nReader=N.getOutput(); %the ouptut is a rawdata reader class (cm2DRawDataReader)


%instantiate the readingFilter
S=readRawDataFilter('signal'); %we need to inform the class that it's going tor ead a signal data so it will reomve the frequency overesampling
S.setFilename(SIGNALFILE); %set the filename
sReader=S.getOutput(); %the ouptut is a rawdata reader class (cm2DRawDataReader)


% we use a7D kspace reppresentation
%    {'1: Average'}    {'2: contrast'}    {'3: repetition'}    {'4: Frequency Encode'}    {'5: Phase Encode'}    {'6: Slice'}    {'7: Coils'}
sK=sReader.getRawDataImageKSpace();

% get the dimensions name of the kspace
display(sReader.getKSpaceDimensionsName())
% get the image dimensions
display(sReader.getImageDimensions())
% get the number of slices
display(sReader.getNumberImageSlices())
% get the number of coils
display( sReader.getNumberCoils())
 % get the number of repetitions
display( sReader.getNumberRepetitions())
 % get the number of averages
display( sReader.getNumberAverages())


%select a slice and reconstruct
average=1;
contrast=1;
repetition=1;
slice=1;

%get slice data 
KS=sReader.getRawDataImageKSpaceSlice(average,contrast,repetition,slice);
KN=nReader.getRawDataImageKSpaceSlice(average,contrast,repetition,slice);


%let's explore the noise covariance matrix and noise 
RE=cm2DRecon();
RE.setNoiseKSpace(KN);
figure
subplot(121);
nc=RE.getNoiseCovariance();
imshow(nc,[]);colorbar();title('Noise Covariance')
subplot(122);
imshow(RE.getNoiseCoefficient,[]);colorbar();title('Noise Coefficient')
clear RE

%let's  reconstruct as Sum of square
%instantiate the reconstruct class
RE= cm2DReconRSS();
RE.setSignalKSpace(KS)
RE.setNoiseKSpace(KN);
figure();
imshow(abs(RE.getOutput()),[]);
colorbar();
title('Reconstruction RSS')


% it is poossible even to set a noise covartiance matrix instead of the
% noise kspace
RE2= cm2DReconRSS();
RE2.setSignalKSpace(KS)
RE2.setNoiseCovariance(nc)
figure();
imshow(abs(RE.getOutput()),[]);
colorbar();
title('Reconstruction RSS using noise covariance instead of the noise kspace')

%B1reconstruction
TITLE='B1 '  
RE= cm2DReconB1();
RE.setSignalKSpace(KS)
RE.setNoiseKSpace(KN)
RE.setCoilSensitivityMatrixSource(KS)
RE.setCoilSensitivityMatrixCalculationMethod('simplesense');
figure();
imshow(abs(RE.getOutput()),[]);colorbar();
colorbar();title(['Reconstruction ' TITLE])


R1=4; %acceleration on frequency direction
R2=2; % acceleration of phase direction
ACL=32; %number of autocalibration lines for sense reconstruction
ACLG=[ACL ACL]; %number of autocalibration lines for grappa reconstruction (frequency, phase)
GK=[7 8]; %grappa kernel
SOSR=1; %sumofsquare 


%GRAPPA
TITLE='GRAPPA '  
RE= cm2DReconGRAPPA();
%mimik a grappa acquisition from a fully sample
k_temp=RE.mimicmGrappaDataFromFullysampledZeroPadded(KS,R1,R2,ACL,ACL);
RE.setSignalKSpace(k_temp)
RE.setNoiseKSpace(KN)
RE.setAcceleration(R2)
RE.setAutocalibration(ACLG)
RE.setGrappaKernel(GK);
RE.setSosRecon(SOSR);
figure();
imshow(abs(RE.getOutput()),[]);
colorbar();title(['Reconstruction ' TITLE])



%SENSE
TITLE='SENSE ' ;
RE=cm2DReconSENSE();
ACCK=RE.mimicmSenseDataFromFullysampledZeroPadded(KS,R1,R2,ACL);
RE.setSignalKSpace(ACCK)
RE.setNoiseKSpace(KN)
RE.setCoilSensitivityMatrixCalculationMethod('simplesense');
%don't need to set the sensitivity source since it will be the signal for simple
RE.setAccelerationFrequency(R1)
RE.setAccelerationPhase(R2);
RE.setAutocalibration(ACL);
OUT.sense=RE.getOutput();
figure();
imshow(abs(RE.getOutput()),[]);
colorbar();title(['Reconstruction ' TITLE])
