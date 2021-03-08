//https://itk.org/Wiki/ITK/Examples/Developer/ImageFilter



#ifndef __cmSNRMultipleReplicasFilter_h
#define __cmSNRMultipleReplicasFilter_h

#include "itkImageToImageFilter.h"
#include "itkNumericTraits.h"

#include "itkImageRegionIteratorWithIndex.h"


namespace cm
{

/** \class SNRMultipleReplicasFilter
 * \brief This class calculates SNR for Multiple replicas.
 * replicas are stored in a std::vectorobj
 *
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 *
 */


template< class ScalarImageType>
class SNRMultipleReplicasFilter:public itk::ImageToImageFilter< ScalarImageType,ScalarImageType >
{
public:
	/* Standard class typedefs. (similar to using..)*/
	typedef SNRMultipleReplicasFilter                                                 Self;
	typedef itk::ImageToImageFilter< ScalarImageType,ScalarImageType > Superclass;
	typedef itk::SmartPointer< Self >                                Pointer;
	typedef itk::SmartPointer< const Self >                                ConstPointer;

	/* Method for creation through the object factory. */
	itkNewMacro(Self);

	/* Run-time type information (and related methods). */
	itkTypeMacro(cm2DSNRArrayCombiningMethod, ImageToImageFilter);




	typedef typename ScalarImageType::Pointer ScalarImageTypePointer;
	typedef typename ScalarImageType::PixelType ScalarImagePixelType;
	typedef typename ScalarImageType::InternalPixelType ScalarImageInternalPixelType;
	typedef typename ScalarImageType::IndexType ScalarImageIndexType;

	typedef typename std::vector<ScalarImageTypePointer> 	VectorofScalarImageType;
	typedef typename std::vector<ScalarImageInternalPixelType> 	VectorofScalarInternalImagePixelType;
   typedef typename itk::ImageRegionIteratorWithIndex<ScalarImageType>        ScalarImageIteratorType;
	typedef typename ScalarImageType::RegionType          ScalarImageRegionType;
	//  I need the first image
	itkGetMacro(stack,VectorofScalarImageType);
	itkSetMacro(stack,VectorofScalarImageType);





	void	pushReconstructedImage(ScalarImageTypePointer im);

	VectorofScalarInternalImagePixelType	getReplicasonPoint(ScalarImageIndexType p);
	std::vector<float>	getReplicasonPointReal(ScalarImageIndexType p);

protected:
	SNRMultipleReplicasFilter(){}
	~SNRMultipleReplicasFilter(){}


	  /* Does the real work. */
//	  virtual void GenerateData();

	virtual void ThreadedGenerateData(const ScalarImageRegionType& outputRegionForThread,
	      itk::ThreadIdType threadId);



private:
	SNRMultipleReplicasFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented

	VectorofScalarImageType m_stack;


};
} //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmSNRMultipleReplicasFilter.hxx"
#endif


#endif // __cmSNRMultipleReplicasFilter_h
