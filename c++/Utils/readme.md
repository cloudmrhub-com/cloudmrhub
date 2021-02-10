will have all the calsses and functionsmainly for cm

consider the Kspace as 3D matrix, that can be used both for a bidimensional acquisition and a threeD.
most part of the filter use the cm::Acquiisiton flag to automatically use the 2D or the 3D filter. 
Check it before going to 3D.
Right now the version is only for multislice 2D;


#NYi -> Cmclasses
ISMRMR/JSNO
- ISMRMRDTomage
- covarianceRiccardo
- differenceImageSNR
- sliceWriteJson 3 times

-writemage



# utils ->ImageUtils
- KSpace2DToVNLVector
- VNLVctorTOSTDVector
- axiallySliceThisImage
- cast2DVNLMatrixfromto
- compyREferenceInfooonImage
- GetPixelNumber
- noiseBand
- readImage
- shrinkImage
- vectorImageToVNLMatrix
- writeImage
- ForwardFFT
- InverseFFT
- ConstPointerToPointer
- addImageAsDimensionToVectorImage
- PixelCountOld
- VectorImageElementAsArray
- VectorImageElementAsImage
- fftshift
- ImageToArray
- multiplyImageTimesScalar
- orImages
- castImage
- composeVectorImage
- stdVectorToVnl
- subtractImage
