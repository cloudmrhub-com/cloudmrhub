# Mr Optimum matlab implementation
[![sssN|Solid](http://cloudmrhub.com/img/cloud.png)](https://https://cloudmrhub.com) 

## Matlab sudgested coding style:
v1
- classes:
    - methods first charachter of each word composing the function name shoud be capital, only the first one is lower case: 
        - getTheFile
        - completeThatArticleYouStartedIn2009
    - attributes first charachter of each word composing the function name shoud be capital
        - SubType
        - ThatThing
- functions: all lowercases
    -    tranformthatthinginthisother
- constants: all capital
    - MYPI
    - CLOUDMRWEBSITE
## Classes
*= mandatory argument
- ### cmOutput
    - base Class here you find all the common method needed from our applications
        -  Attributes(Access=private)
            - Exporter % a cell array of name and values to be exported
            - Log
            - Type % Output type
            - SubType % Output subtype (if needed)
            - OUTPUTLOGFILENAME
            - OUTPUTFILENAME
        - Methods
            - setOutputLogFileName(*filename)
            - getOutputLogFileName
            - setOutputFileName(*filename)
            - getOutputFileName
            - setLog(log)
            - getLog
            - appendLog(*log)
            - logIt(*text,*type)
            - addToExporter(*type,*name,*value)
            - getResults
            - exportResults(filename)
            - exportLog(filename)
            - whatHappened()
            - errorMessage %standrd ERROR message for the log
            - add2DImagetoExport(*im,*name)
            - getResultsAsImages()
        - Methods (static)
            - image2DToJson(*im2d)
            - get2DKSIFFT(*2dkspace)
            - bartMy2DKSpace(*2dkspace,cartesianFLAG)
            - debartMy2DKSpace(*2dkspace,cartesianFLAG)
            - rescale01(*array1D)
            - getJsonResultFromJsonFile(*jsonresultfile)
            - getImagsesFromJsonResultFile(*jsonresultfile)
            - getImagesFromResuls(*jsondecodedresult)
            - write2DCartesianKspacedatainISMRMRDv1(*fullysampled2DKspce,*filename)
            - getMRoptimumTestData() %routine to test the class
            - plotImageAfterTest(*image,*plottitle)
            - prewhiteningSignal(*signalrawdata,*psi)
            - rescale01(*array)
            - shrinktoaliasedmatrix_2d(*KSPACE,R1,R2)
            - test()

## Recon Classes
### cm2DRecon<cmOutput 
 base recon class here you find all the common method needed from our applications
-  Attributes(Access=private)
    - SignalKSpace %freq,phase,coils
    - SignalNCoils %freq,phase,coils
    - SignalSize %freq,phase
    - NoiseCovariance %ncoils.ncoils
    - InverseNoiseCovariance %ncoils.ncoils
    - SignalPrewhitened %freq,phase,coils    
- Methods
    - getHasAcceleration
    - getHasSensitivity
    - getInverseNoiseCovariance
    - getInverseNoiseCovariancePrewhithened
    - getNoiseCovariance
    - getNoiseCovariancePrewhithened
    - getPrewhitenedSignal
    - getSignalKSpace
    - getSignalNCoils
    - getSignalSize
    - setInverseNoiseCovariance
    - setNoiseCovariance
    - setPrewhitenedSignal
    - setSignalKSpace
    - setSignalNCoils
    - setSignalSize
    - test                                  

### cm2DReconRSS<cm2DRecon Root sum of squares class 
-  Methods
    -  getOutput() reconstruct RSS
### cm2DReconWithSensitivity<cm2DRecon base class for reconstruciton with fully smapled kspace (B1 weighted)
    - Attibutes (Access=private)
        CoilSensitivityMatrix %Smap 
        CoilSensitivityMatrixSourcePrewhitened %the matrix for the sensitivity calculation to be set prewhitened
        CoilSensitivityMatrixSource %the matrix for the sensitivity calculation to be set
        CoilSensitivityMatrixCalculationMethod %'simplesense' ,'adaptive','espirit'
        CoilSensitivityMatrixSourceNCoils %number of coils set with the source
        CoilSensitivityMatrixSourceSmooth=false
        WithSensitivity=true;        
    end
-  Methods
    - getCoilSensitivityMatrix                   
    - getCoilSensitivityMatrixCalculationMethod  
    - getCoilSensitivityMatrixSource             
    - getCoilSensitivityMatrixSourceNCoils       
    - getCoilSensitivityMatrixSourcePrewhitened  
    - getCoilSensitivityMatrixSourceSmooth       
    - getWithSensitivity                         
    - resetCoilSensitivityMatrix                 
    - setCoilSensitivityMatrix(*sensitivitymatrix)                   
    - setCoilSensitivityMatrixCalculationMethod(*prewhithenedsensitivitymatrix)  
    - setCoilSensitivityMatrixSource(*sourcefilekspace)  
    - setCoilSensitivityMatrixSourceNCoils(*numbercoil)       
    - setCoilSensitivityMatrixSourcePrewhitened(*prewhitenedsource2dkspace)  
    - setCoilSensitivityMatrixSourceSmooth       
    - setWithSensitivity(*boolean)       
    - testSensitivityMatrixvalidity         
### cm2DReconB1<cm2DReconWithSensitivity
reconstruct images as b1
    -  Methods
        -  getOutput() reconstruct RSS
### cm2DReconWithSensitivityAutocalibrated<cm2DReconWithSensitivity
base class for reconstrcution of images that need senitivity map calculated from accelerated kspace
    -  Attributes (Access=private)
        - AutocalibrationF
        - AutocalibrationP
        - AccelerationF
        - AccelerationP
    - Methods
        - getAccelerationFrequency                
        - getAccelerationPhase                    
        - getAutocalibrationFrequency             
        - getAutocalibrationPhase                 
        - setAccelerationFrequency                
        - setAccelerationPhase                    
        - setAutocalibrationFrequency             
        - setAutocalibrationPhase  



### cm2DReconGRAPPA<cm2DReconWithSensitivityAutocalibrated 
    - Methods (access=private)
        - getAccelerationFrequency                  
        - getAutocalibration                        
        - getGrappaKernel                           
        - getOutput                                 
        - getR                                      
        - getSignalPrewhitenedAutocalibrationsArea  
        - getSosRecon                               
        - setAccelerationFrequency                  
        - setAutocalibration                        
        - setGrappaKernel                           
        - setSosRecon                               

    - Static methods:
        - getGRAPPAKspace       
## Kellman SNR unit reconstrutions
###  cm2DKellmanRSS<cm2DReconRSS
reconstruct in snr units as rss fully sampled kspace
- Methods
    - getOutput()
### cm2DKellmanB1<cm2DReconB1

- Methods
-   getOutput()


### cm2DKellmanB1<cm2DReconB1

- Methods
-   getOutput()


### cm2DKellmanSENSE<cm2DReconSENSE
- Methods
-   getOutput()


## G factor ricotruction
### cm2DGFactorSENSE<cm2DReconSENSE 
- Methods
-   getOutput()


# Requirements
- matlab 2019x
- 4GB RAM

# Tests
against Riccardo's code
- [x] Recon Grappa 
- [x] Recon SENSE
	- [x] set signal, nc and sensitiviy
	- [x] set signal, nc and source senritivity
	- [x] set sognal, nc andsource sensitivity as simplesense
- [x] Recon B1
- [x] Recon RSS 










[*Dr. Eros Montin, PhD*](http://me.biodimensional.com)
**46&2 just ahead of me!**
