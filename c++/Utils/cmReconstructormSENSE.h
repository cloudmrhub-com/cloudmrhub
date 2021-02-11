#ifndef __cmReconstructormSENSE_h
#define __cmReconstructormSENSE_h

#pragma once
#include "cm.h"
#include "cmReconstructorWithSensitivityAccelerated.h"

#include "itkShrinkImageFilter.h"

namespace cm
{

/** \class ReconstructormSENSE
 * \brief This class calculates SNR on a 2D k-space phased array image using SENSE.
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 */


template< class VectorImageType,class ScalarImageType>
class ReconstructormSENSE:public ReconstructorWithSensitivityAccelerated< VectorImageType, ScalarImageType >
{
public:
	/* Standard class typedefs. (similar to using..)*/
	typedef ReconstructormSENSE                                                 Self;
	typedef ReconstructorWithSensitivityAccelerated< VectorImageType, ScalarImageType > Superclass;
	typedef itk::SmartPointer< Self >                                Pointer;
	typedef itk::SmartPointer< const Self >                                ConstPointer;

	/* Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(ReconstructormSENSE,   ReconstructorWithSensitivityAccelerated);



	typedef typename VectorImageType::Pointer VectorImageTypePointer;
	typedef typename VectorImageType::PixelType VectorImagePixelType;
	typedef typename VectorImageType::InternalPixelType VectorImageInternalPixelType;
	typedef typename VectorImageType::RegionType          VectorImageRegionType;
	typedef typename itk::ImageRegionIterator<VectorImageType>        VectorImageIteratorType;


	typedef typename ScalarImageType::Pointer ScalarImageTypePointer;
	typedef typename ScalarImageType::PixelType ScalarImagePixelType;
	typedef typename ScalarImageType::InternalPixelType scalarImageInternalPixelType;
   
// OutputImageTypePointer GetGfactor();

protected:
	ReconstructormSENSE(){
this->SetHasAcceleration(true);
this->SetHasSensitivity(true);
	}
	~ReconstructormSENSE(){}


//   	itkSetMacro(ImFolded,VectorImageTypePointer);
//   	itkGetMacro(ImFolded,VectorImageTypePointer);
	 virtual void GenerateData();
//	virtual void ThreadedGenerateData(const VectorImageRegionType& outputRegionForThread,	itk::ThreadIdType threadId);


	//calculate the ifftt and the sos
	void BeforeThreadedGenerateData( ){
//		calculate the sensitivity map
//		this->GetSensitivityMap(); //TODO reewrtite it
//		this->GetInputIFFT();//TODO reewrtite it
//
//		typename VectorImageType::SizeType acceleration=this->GetAcceleration();
//		const int Kdimension =this->GetInput()->GetNumberOfComponentsPerPixel();
//		cm::VectorImagePixelElelmentType c(Kdimension);
//		cm::PixelType RTOT=(acceleration[0]*acceleration[1],0.0);
//		c.Fill(RTOT);
//		this->SetInputIFFT(multiplyImageTimesScalar(this->GetInputIFFT(),RTOT));
//
//		//FOLDE IT
//    	  using ShrinkImageFilterType = itk::ShrinkImageFilter<VectorImageType ,VectorImageType>;
//    	  typename ShrinkImageFilterType::Pointer shrinkFilter = ShrinkImageFilterType::New();
//    	    shrinkFilter->SetInput(this->GetInputIFFT());
//    	    shrinkFilter->SetShrinkFactor(0, acceleration[0]);
//    	    shrinkFilter->SetShrinkFactor(1, acceleration[1]);
//    	    shrinkFilter->Update();
//		this->SetInputIFFT(shrinkFilter=>GetOutput());
//		this->CalculateInverseCovariance();


	};


  /* this filter need to expand the image */
//  virtual void GenerateOutputInformation() ITK_OVERRIDE;

  /* this filter need to expand the image */
//  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;
  
  
private:
	ReconstructormSENSE(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
	
	
//	void senseRecon2D(const VectorImageRegionType& outputRegionForThread, itk::ThreadIdType threadId);
	void senseRecon2D();
//	VectorImageType m_ImFolded;
};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmReconstructormSENSE.hxx"
#endif


#endif // __ReconstructormSENSE_h
