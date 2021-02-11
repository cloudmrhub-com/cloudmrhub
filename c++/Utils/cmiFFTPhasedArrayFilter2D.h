//https://itk.org/Wiki/ITK/Examples/Developer/ImageFilter



#ifndef __iFFTPhasedArrayFilter2D_h
#define __iFFTPhasedArrayFilter2D_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "cm.h"

namespace cm
{

/** \class iFFTPhasedArrayFilter2D
 * \brief This class computes the 2DIFFT on a phased array fully sampled k-space data. Need test for 3D
 *
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 *
 */


template< class TImage>
class iFFTPhasedArrayFilter2D:public itk::ImageToImageFilter< TImage, TImage >
{
public:
//	ITK_DISALLOW_COPY_AND_ASSIGN(iFFTPhasedArrayFilter2D);
	/* Standard class typedefs. (similar to using..)*/
	typedef iFFTPhasedArrayFilter2D                                                 Self;
	typedef itk::ImageToImageFilter< TImage, TImage > Superclass;
	typedef itk::SmartPointer< Self >    Pointer;
	typedef itk::SmartPointer< const Self >          ConstPointer;

	/* Method for creation through the object factory. */
	itkNewMacro(Self);

	/* Run-time type information (and related methods). */
	itkTypeMacro(iFFTPhasedArrayFilter2D, ImageToImageFilter);

	typedef TImage                   InputImageType;
	typedef typename TImage::Pointer InputImageTypePointer;
	typedef typename TImage::PixelType InputImagePixelType;
	typedef typename TImage::InternalPixelType InputImageInnerPixelType;



	/** Do you want to normalize by the square root of the number of pixels in the k-space? Default NO*/
	enum NormalizeIFFTType {
		YES = 1,
		NO = 0
	};

	itkSetMacro(Normalize, NormalizeIFFTType);
	itkGetConstMacro(Normalize, NormalizeIFFTType);


	   //Kspacedimension
 	itkGetMacro(KSpaceDimension,cm::KSpaceAcquisitionDimension);
 	itkSetMacro(KSpaceDimension,cm::KSpaceAcquisitionDimension);

protected:
	iFFTPhasedArrayFilter2D(){
		m_Normalize=NormalizeIFFTType::YES;
	}
	~iFFTPhasedArrayFilter2D(){}

	/* Does the real work. */
	virtual void GenerateData();


private:
//	iFFTPhasedArrayFilter2D(const Self &); //purposely not implemented
//	void operator=(const Self &);  //purposely not implemented

	NormalizeIFFTType m_Normalize;
	cm::KSpaceAcquisitionDimension m_KSpaceDimension;


};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmiFFTPhasedArrayFilter2D.hxx"
#endif


#endif // __cmiFFTPhasedArrayFilter2D_h
