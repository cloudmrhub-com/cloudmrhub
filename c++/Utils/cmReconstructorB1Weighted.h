#ifndef __cmReconstructorB1Weighted_h
#define __cmReconstructorB1Weighted_h

#pragma once
#include "cm.h"
#include "cmReconstructorWithSensitivity.h"



namespace cm
{

/** \class ReconstructorB1Weighted
 * \brief This class calculates SNR on a 2D k-space phased array image using B1.
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 */


template< class VectorImageType,class ScalarImageType>
class ReconstructorB1Weighted:public ReconstructorWithSensitivity< VectorImageType, ScalarImageType >
{
public:
	/* Standard class typedefs. (similar to using..)*/
	typedef ReconstructorB1Weighted                                                 Self;
	typedef ReconstructorWithSensitivity< VectorImageType, ScalarImageType > Superclass;
	typedef itk::SmartPointer< Self >                                Pointer;
	typedef itk::SmartPointer< const Self >                                ConstPointer;

	/* Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(ReconstructorB1Weighted, ReconstructorWithSensitivity);


	/*  the main image type*/

	typedef typename VectorImageType::Pointer VectorImageTypePointer;
	typedef typename VectorImageType::PixelType VectorImagePixelType;
	typedef typename VectorImageType::InternalPixelType VectorImageInternalPixelType;
	typedef typename VectorImageType::RegionType          VectorImageRegionType;
	typedef typename itk::ImageRegionIterator<VectorImageType>        VectorImageIteratorType;


	typedef typename ScalarImageType::Pointer ScalarImageTypePointer;
	typedef typename ScalarImageType::PixelType ScalarImagePixelType;
	typedef typename ScalarImageType::InternalPixelType scalarImageInternalPixelType;

	typedef  vnl_matrix<VectorImageInternalPixelType> ChannelArrayType;






protected:
	ReconstructorB1Weighted(){
		  this->SetHasAcceleration(false);
		  this->SetHasSensitivity(true);
	}
	~ReconstructorB1Weighted(){}



	// virtual void GenerateData();
	virtual void ThreadedGenerateData(const VectorImageRegionType& outputRegionForThread,
			itk::ThreadIdType threadId);


	//calculate the ifftt and the sos
	void BeforeThreadedGenerateData( ){
		//calculate the sensitivity map
		this->GetSensitivityMap();
		this->GetInputIFFT();
		this->CalculateInverseCovariance();
	};

private:
	ReconstructorB1Weighted(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented



};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmReconstructorB1Weighted.hxx"
#endif


#endif // __ReconstructorB1Weighted_h
