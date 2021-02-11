//https://itk.org/Wiki/ITK/Examples/Developer/ImageFilter



#ifndef __cmPhasedArraySensitivityMapFilter_h
#define __cmPhasedArraySensitivityMapFilter_h

#include "itkImageToImageFilter.h"
#include <typeinfo>
#include "itkImageRegionIterator.h"
#include "cm.h"
namespace cm
{

/** \class PhasedArraySensitivityMapFilter
 * \brief This class represents the abstarct class of phased array sensitivity map.
 * in order to speed up the run time exectution, m_Ifft is a pointer to the input ifft so that it can be calculated just once before the threaded part.  cm::KSpaceAcquisitionDimension m_KSpaceDimension can be bidimensional or    
 * threedimensional as describe din cm.h. rescaleOutputVectorImage() is a method that rescale the coil sensitivity from 0 to 1 by dividing each moponent maps by its max. where max is intended as the max of abs.
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 *
 *
 */


template< class VectorImageType>
class PhasedArraySensitivityMapFilter:public itk::ImageToImageFilter<VectorImageType, VectorImageType >
{
public:
  /* Standard class typedefs. (similar to using..)*/
  typedef PhasedArraySensitivityMapFilter                                                 Self;
  typedef itk::ImageToImageFilter< VectorImageType, VectorImageType > Superclass;
  typedef itk::SmartPointer< Self >                                Pointer;
  typedef itk::SmartPointer< const Self >                                ConstPointer;

  /* Method for creation through the object factory. */
  itkNewMacro(Self);

  /* Run-time type information (and related methods). */
  itkTypeMacro(PhasedArraySensitivityMapFilter, ImageToImageFilter);

/*  kljk the main image type*/
  
   typedef typename VectorImageType::Pointer VectorImageTypePointer;
   typedef typename VectorImageType::PixelType VectorImagePixelType;
   typedef typename VectorImageType::InternalPixelType VectorImageInternalPixelType;
   typedef typename VectorImageType::RegionType          VectorImageRegionType;
    typedef typename itk::ImageRegionIterator<VectorImageType>        VectorImageIteratorType;

	
	void rescaleOutputVectorImage();

		itkGetMacro(Ifft,VectorImageTypePointer);
		itkSetMacro(Ifft,VectorImageTypePointer);

		itkGetMacro(KSpaceDimension,	cm::KSpaceAcquisitionDimension);
		itkSetMacro(KSpaceDimension,	cm::KSpaceAcquisitionDimension);


protected:

	PhasedArraySensitivityMapFilter(){};
  ~PhasedArraySensitivityMapFilter(){};



private:
  PhasedArraySensitivityMapFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

  
	VectorImageTypePointer m_Ifft;


	cm::KSpaceAcquisitionDimension m_KSpaceDimension;

};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmPhasedArraySensitivityMapFilter.hxx"
#endif


#endif // __cmPhasedArraySensitivityMapFilter_h
