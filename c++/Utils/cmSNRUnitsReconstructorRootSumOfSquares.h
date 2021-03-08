#ifndef __cmSNRUnitsReconstructorRootSumOfSquares_h
#define __cmSNRUnitsReconstructorRootSumOfSquares_h

#pragma once
#include "cm.h"
#include "cmReconstructorRootSumOfSquares.h"



namespace cm
{

/** \class SNRUnitsReconstructorRootSumOfSquares
 * \brief This class calculates SNR on a 2D k-space phased array image using Root sum of squares (RSS).
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 */


template< class TImage,class TOImage>
class SNRUnitsReconstructorRootSumOfSquares:public ReconstructorRootSumOfSquares< TImage, TOImage >
{
public:
  /* Standard class typedefs. (similar to using..)*/
  typedef SNRUnitsReconstructorRootSumOfSquares                                                 Self;
  typedef ReconstructorRootSumOfSquares< TImage, TOImage > Superclass;
  typedef itk::SmartPointer< Self >                                Pointer;
  typedef itk::SmartPointer< const Self >                                ConstPointer;

  /* Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SNRUnitsReconstructorRootSumOfSquares, ReconstructorRootSumOfSquares);


  /*  the main image type*/
     typedef TImage                   VectorImageType;
     typedef typename TImage::Pointer VectorImageTypePointer;   
     typedef typename TImage::PixelType VectorImagePixelType;
     typedef typename TImage::InternalPixelType VectorImageInternalPixelType;   
 	 typedef typename VectorImageType::RegionType          VectorImageRegionType;
 	 typedef typename itk::ImageRegionIterator<VectorImageType>        VectorImageIteratorType;

     typedef TOImage                   ScalarImageType;
     typedef typename TOImage::Pointer ScalarImagePointerType;
     typedef typename TOImage::PixelType ScalarImagePixelType;
     typedef typename TOImage::InternalPixelType scalarImageInternalPixelType;

     typedef  vnl_matrix<VectorImageInternalPixelType> ChannelArrayType;







protected:
     SNRUnitsReconstructorRootSumOfSquares() : Superclass() {}
  ~SNRUnitsReconstructorRootSumOfSquares(){}


  // virtual void GenerateData();
  	virtual void ThreadedGenerateData(const VectorImageRegionType& outputRegionForThread,
  			itk::ThreadIdType threadId);


  	//calculate the ifftt and the sos
  	void BeforeThreadedGenerateData( ){
  		//calculate the sensitivity map
//  		this->GetSensitivityMap();
  		this->GetInputIFFT();
  	};

private:
  SNRUnitsReconstructorRootSumOfSquares(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented


};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmSNRUnitsReconstructorRootSumOfSquares.hxx"
#endif


#endif // __SNRUnitsReconstructorRootSumOfSquares_h
