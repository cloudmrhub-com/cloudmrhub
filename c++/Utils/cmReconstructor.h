//https://itk.org/Wiki/ITK/Examples/Developer/ImageFilter



#ifndef __cmReconstructor_h
#define __cmReconstructor_h

#include "itkImageToImageFilter.h"
#include "cm.h"

namespace cm
{

/** \class Reconstructor
 * \brief Base class for phased array image.
 *	input is a signal kspace and noise is a 2x vnl_matrix size of number of channles X number of pixels.
 *	Noise covariance matrx may be corrected  (SetNoiseBandWidthCorrection(true))by a value (SetNoiseBandWidth(XX)) if not set SetNoiseBandWidth will be calculated.
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 *
 */


template< class TImage,class TOImage>
class Reconstructor:public itk::ImageToImageFilter< TImage, TOImage >
{
public:
  /* Standard class typedefs. (similar to using..)*/
  typedef Reconstructor                                                 Self;
  typedef itk::ImageToImageFilter< TImage, TOImage > Superclass;
  typedef itk::SmartPointer< Self >                                Pointer;
  typedef itk::SmartPointer< const Self >                                ConstPointer;

  /* Method for creation through the object factory. */
  itkNewMacro(Self);

  /* Run-time type information (and related methods). */
  itkTypeMacro(Reconstructor, ImageToImageFilter);

/*  the main image type*/
   typedef TImage                   VectorImageType;
   typedef typename TImage::Pointer VectorImageTypePointer;
   typedef typename TImage::PixelType VectorImagePixelType;
   typedef typename TImage::InternalPixelType VectorImageInternalPixelType;

   typedef TOImage                   ScalarImageType;
   typedef typename TOImage::Pointer ScalarImagePointerType;
   typedef typename TOImage::PixelType ScalarImagePixelType;
   typedef typename TOImage::InternalPixelType scalarImageInternalPixelType;

   typedef  vnl_matrix<VectorImageInternalPixelType> ChannelArrayType;


      void SetNoiseKSpace(VectorImageTypePointer Noise);
      itkGetMacro(NoiseKSpace,VectorImageTypePointer);
      //itkSetMacro(NoiseKSpace,VectorImageTypePointer); This method is implemented, every time a noise kspace is set
      // the Noise matrix is updated
   //   this is the noise kspace reshaped as 2D matrix of number of channels x number of pixels
      itkGetMacro(Noise,ChannelArrayType);
      itkSetMacro(Noise,ChannelArrayType)


//  NBW correction
itkGetMacro(NoiseBandWidthCorrection,bool);
itkSetMacro(NoiseBandWidthCorrection,bool);

//  NBW
itkGetMacro(NoiseBandWidth,float);
itkSetMacro(NoiseBandWidth,float);



//itkGetMacro(NoiseCovarianceMatrix,ChannelArrayType);
itkSetMacro(NoiseCovarianceMatrix,ChannelArrayType);
cm::ChannelArrayType GetNoiseCovarianceMatrix();


itkGetMacro(NoiseCoefficientMatrix,ChannelArrayType);
itkSetMacro(NoiseCoefficientMatrix,ChannelArrayType);

//Kspacedimension
 	itkGetMacro(KSpaceDimension,cm::KSpaceAcquisitionDimension);
 	itkSetMacro(KSpaceDimension,cm::KSpaceAcquisitionDimension);


//IFFT of the input image useful for threads
 	//itkGetMacro(InputIFFT,VectorImageTypePointer);
 	itkSetMacro(InputIFFT,VectorImageTypePointer);
 	VectorImageTypePointer GetInputIFFT();

ChannelArrayType CalculateNoiseCovarianceMatrix();
VectorImageInternalPixelType CalculateNoiseBW();
ChannelArrayType CalculateNoiseCoefficientMatrix();

//public can query but can't set it:)
itkGetMacro(HasSensitivity,bool);

itkGetMacro(HasAcceleration,bool);

virtual void SetSensitivityMapCalculationMode(cm::SensitivityMapCalculation x){};

virtual void SetAcceleration(cm::VectorImageType::IndexType a){};
virtual void SetAutocalibrationLines(cm::VectorImageType::IndexType a){};

protected:
  Reconstructor(){
	  
	  this->m_NoiseBandWidth=std::nanf("1");
	  
	  this->m_InputIFFT= nullptr;
	  
	  this->m_NoiseBandWidthCorrection=true;

	  this->m_NoiseCoefficientMatrix= std::nanf("1");;

	  this->m_NoiseCovarianceMatrix= std::nanf("1");;

  }

  ~Reconstructor(){}

  //public can query but can't set it:)
  itkSetMacro(HasSensitivity,bool);

  itkSetMacro(HasAcceleration,bool);


VectorImageTypePointer m_InputIFFT;

private:
  Reconstructor(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

  /** noise is a 2D matrix of Nchannel X Nvoxels*/
  ChannelArrayType  m_Noise;

  bool m_NoiseBandWidthCorrection;
  //instantiated as nan if is nan is calculated
  float m_NoiseBandWidth;
  
  ChannelArrayType m_NoiseCovarianceMatrix;

  ChannelArrayType m_NoiseCoefficientMatrix;
  //  KSpace of noise
    VectorImageTypePointer m_NoiseKSpace;

  //2D or 3D sequence
  cm::KSpaceAcquisitionDimension m_KSpaceDimension;

  bool m_HasAcceleration;
  bool m_HasSensitivity;


    
};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmReconstructor.hxx"
#endif


#endif // __cmReconstructor_h
