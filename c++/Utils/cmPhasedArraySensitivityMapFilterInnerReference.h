//https://itk.org/Wiki/ITK/Examples/Developer/ImageFilter



#ifndef __cmPhasedArraySensitivityMapFilterInnerReference_h
#define __cmPhasedArraySensitivityMapFilterInnerReference_h

#include "cmPhasedArraySensitivityMapFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImage.h"
#include "cm.h"
#include "cmiFFTPhasedArrayFilter.h"

namespace cm
{

/** \class PhasedArraySensitivityMapFilterInnerReference
 * \brief This class calculates the sensitivity maps as described by Riccardo Lattanzi code.
 *input is the signal Kspace data ifftof the phased array is sivided by the repetition of the sos;

 *
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 *
 */


template< class VectorImageType>
class PhasedArraySensitivityMapFilterInnerReference:public PhasedArraySensitivityMapFilter< VectorImageType >
{
public:
	/* Standard class typedefs. (similar to using..)*/
	typedef PhasedArraySensitivityMapFilterInnerReference                                                 Self;
	typedef PhasedArraySensitivityMapFilter<VectorImageType > Superclass;
	typedef itk::SmartPointer< Self >                                Pointer;
	typedef itk::SmartPointer< const Self >                                ConstPointer;

	/* Method for creation through the object factory. */
	itkNewMacro(Self);

	/* Run-time type information (and related methods). */
	itkTypeMacro(PhasedArraySensitivityMapFilterInnerReference, PhasedArraySensitivityMapFilter);

	/*the main image type*/
	typedef typename VectorImageType::Pointer VectorImageTypePointer;
	typedef typename VectorImageType::PixelType VectorImagePixelType;
	typedef typename VectorImageType::InternalPixelType VectorImageInternalPixelType;
	typedef typename VectorImageType::RegionType          VectorImageRegionType;
	typedef typename itk::ImageRegionIterator<VectorImageType>        VectorImageIteratorType;

	typedef typename itk::Image<typename VectorImageType::InternalPixelType,3> ScalarImageType;
	typedef typename ScalarImageType::Pointer ScalarImageTypePointer;
	typedef typename ScalarImageType::PixelType ScalarImagePixelType;


	itkGetMacro(Sos,ScalarImageTypePointer);
	itkSetMacro(Sos,ScalarImageTypePointer);





protected:
	PhasedArraySensitivityMapFilterInnerReference(){

	}
	~PhasedArraySensitivityMapFilterInnerReference(){}

	/* Does the real work. */
	// virtual void GenerateData();
	virtual void ThreadedGenerateData(const VectorImageRegionType& outputRegionForThread,
			itk::ThreadIdType threadId);

	//calculate the ifftt and the sos
	void BeforeThreadedGenerateData( );


private:
	PhasedArraySensitivityMapFilterInnerReference(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented


	ScalarImageTypePointer m_Sos;





};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmPhasedArraySensitivityMapFilterInnerReference.hxx"
#endif


#endif // __cmPhasedArraySensitivityMapFilterInnerReference_h
