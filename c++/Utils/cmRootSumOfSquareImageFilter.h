//https://itk.org/Wiki/ITK/Examples/Developer/ImageFilter



#ifndef __cmRootSumOfSquareImageFilter_h
#define __cmRootSumOfSquareImageFilter_h

#include "itkImageToImageFilter.h"
#include <typeinfo>
#include "cm.h"

namespace cm
{

/** \class RootSumOfSquareImageFilter
 * \This Filter reconstructs image from 2D fully sampled  k-space phased array raw data as Root Sum of Squares.
 * input image is a multicoil vector image of complex number the output is a scalar matrix of complex.
 * the acquisition can be 2d or 3d and it will change the way the ifft is computed.
 * in a 2d acquisition we would do a 2D ifft. the 3D is yet to develop
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 */


template< class TImage,class TOImage>
class RootSumOfSquareImageFilter:public itk::ImageToImageFilter< TImage, TOImage >
{
public:
  /* Standard class typedefs. (similar to using..)*/
  typedef RootSumOfSquareImageFilter                                                 Self;
  typedef itk::ImageToImageFilter< TImage, TImage > Superclass;
  typedef itk::SmartPointer< Self >                                Pointer;
  typedef itk:: SmartPointer< const Self >                                ConstPointer;

  /* Method for creation through the object factory. */
  itkNewMacro(Self);

  /* Run-time type information (and related methods). */
  itkTypeMacro(RootSumOfSquareImageFilter, ImageToImageFilter);


   
   /*  the main image type*/
   typedef TImage                   VectorImageType;
   typedef typename TImage::Pointer VectorImageTypePointer;
   typedef typename TImage::PixelType VectorImagePixelType;
   typedef typename TImage::InternalPixelType VectorImageInternalPixelType;

   typedef TOImage                   ScalarImageType;
   typedef typename TOImage::Pointer ScalarImagePointerType;
   typedef typename TOImage::PixelType ScalarImagePixelType;
   typedef typename TOImage::InternalPixelType scalarImageInternalPixelType;


   	   //Kspacedimension
    	itkGetMacro(KSpaceDimension,cm::KSpaceAcquisitionDimension);
    	itkSetMacro(KSpaceDimension,cm::KSpaceAcquisitionDimension);



  
protected:
  RootSumOfSquareImageFilter(){
	  std::cout<<"RSOS with scalarImageInternalPixelType type: " << typeid(scalarImageInternalPixelType).name()<<std::endl;

  }
  ~RootSumOfSquareImageFilter(){}

  /* Does the real work. */
  virtual void GenerateData();
  

private:
  RootSumOfSquareImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
  cm::KSpaceAcquisitionDimension m_KSpaceDimension;


};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmRootSumOfSquareImageFilter.hxx"
#endif


#endif // __RootSumOfSquareImageFilter_h
