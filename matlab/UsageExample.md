# Array Combining
The Classes for kellman array combining are 3 [cit]:
- cm2DKellmanRSS
- cm2DKellmanB1
- cm2DKellmanSENSE

Similar to the recons classes you can use the classes in different way:
- setting data (Signal) KSpace and a noise covariance matrix
- setting prewhithened data
- setting Signal Kspace and Noise Kspace

SNR will be available by calling CLASS.getOutput() 

## Example cm2DKellmanB1

```
%input data
S= single slice 3D matrix of size #frequencyencoding, #phaseencoding, #coils
N= noise 3D matrix of size #frequencyencoding, #phaseencoding, #coils
NC= noise covariance matrix sie #coils x #coils
PWS = prewhithened single slice 3D matrix of size #frequencyencoding, #phaseencoding, #coils



	AC=cm2DReconB1();
AC.setNoiseCovariance(NC)
AC.setSignalKSpace(S);
AC.setCoilSensitivityMatrixSource(S);
AC.setCoilSensitivityMatrixCalculationMethod('simplesense');
image=AC.getOutput();


AC.setNoiseCovariance(NC)
AC.setSignalKSpace(S);
AC.setCoilSensitivityMatrixSource(S);
AC.setCoilSensitivityMatrixCalculationMethod('simplesense');


%SNR
AC= cm2DKellmanB1();


%alterenative you can use
%AC=cm2DReconB1();
%AC.setSignalKSpace(S);
%AC.setNoiseKSpace(N)

%alterenative 3 you can use
%AC=cm2DReconB1([],[],PW);
%but the sensitivity matrix need to be passed
%AC.AC.setCoilSensitivityMatrixSourcePrewhitened(PW);


```

## Multiple Replicas

Calculate the SNR as multiple replicas



```
%inpput data
% thereplicas is a stck of 2D siingle sliced images

MR=cm2DSignalToNoiseRatioMultipleReplicas(thereplicas);
MR.getOutput()

```


## Pseudo Muslitple Replicas

```
%input data
S= single slice 3D matrix of size #frequencyencoding, #phaseencoding, #coils
N= noise 3D matrix of size #frequencyencoding, #phaseencoding, #coils

%create a reconstructor
RE= cm2DKellmanB1(S,NC);
RE.setCoilSensitivityMatrixSource(S);
RE.setCoilSensitivityMatrixCalculationMethod('simplesense');

PMR=cm2DSignalToNoiseRatioPseudoMultipleReplicas(RE);
PMR.setNumberOfPseudoReplicas(100);
SNR=PMR.getOutput();

PMR=cm2DSignalToNoiseRatioPseudoMultipleReplicasWien(RE);
PMR.setNumberOfPseudoReplicas(100);
SNR=PMR.getOutput();

```
