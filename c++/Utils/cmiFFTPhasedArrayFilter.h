//https://itk.org/Wiki/ITK/Examples/Developer/ImageFilter



#ifndef __cmiFFTPhasedArrayFilter_h
#define __cmiFFTPhasedArrayFilter_h

#include "itkImageToImageFilter.h"

#include "cm.h"

namespace cm
{

/** \class iFFTPhasedArrayFilter
 * \This Filter do an ifft from fully sampled  k-space phased array raw data .
 * input image is a multicoil vector image of complex number the output is a vectorimage matrix of complex.
 * the acquisition can be 2d or 3d and it will change the way the ifft is computed.
 * in a 2d acquisition we would do a 2D ifft. the 3D is yet to be developed
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 */


template< class VectorImageType>
class iFFTPhasedArrayFilter:public itk::ImageToImageFilter< VectorImageType, VectorImageType >
{
public:
  /* Standard class typedefs. (similar to using..)*/
  typedef iFFTPhasedArrayFilter                                                 Self;
  typedef itk::ImageToImageFilter< VectorImageType, VectorImageType > Superclass;
  typedef itk::SmartPointer< Self >                                Pointer;
  typedef itk:: SmartPointer< const Self >                                ConstPointer;

  /* Method for creation through the object factory. */
  itkNewMacro(Self);

  /* Run-time type information (and related methods). */
  itkTypeMacro(iFFTPhasedArrayFilter, ImageToImageFilter);


   
   /*  the main image type*/
 
   typedef typename VectorImageType::Pointer VectorImageTypePointer;
   typedef typename VectorImageType::PixelType VectorImagePixelType;
   typedef typename VectorImageType::InternalPixelType VectorImageInternalPixelType;


   	   //Kspacedimension
    	itkGetMacro(KSpaceDimension,cm::KSpaceAcquisitionDimension);
    	itkSetMacro(KSpaceDimension,cm::KSpaceAcquisitionDimension);


    	void ifft3D();
   	void ifft2D();
  
protected:
  iFFTPhasedArrayFilter(){}
  ~iFFTPhasedArrayFilter(){}

  /* Does the real work. */
  virtual void GenerateData();
  

private:
  iFFTPhasedArrayFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
  cm::KSpaceAcquisitionDimension m_KSpaceDimension;


};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmiFFTPhasedArrayFilter.hxx"
#endif


#endif // __cmiFFTPhasedArrayFilter_h
