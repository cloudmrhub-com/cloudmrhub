#ifndef __cmReconstructorRootSumOfSquares_h
#define __cmReconstructorRootSumOfSquares_h

#pragma once
#include "cm.h"
#include "cmReconstructor.h"



namespace cm
{

/** \class ReconstructorRootSumOfSquares
 * \brief This class calculates SNR on a 2D k-space phased array image using Root sum of squares (RSS).
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 */


template< class TImage,class TOImage>
class ReconstructorRootSumOfSquares:public Reconstructor< TImage, TOImage >
{
public:
  /* Standard class typedefs. (similar to using..)*/
  typedef ReconstructorRootSumOfSquares                                                 Self;
  typedef Reconstructor< TImage, TOImage > Superclass;
  typedef itk::SmartPointer< Self >                                Pointer;
  typedef itk::SmartPointer< const Self >                                ConstPointer;

  /* Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ReconstructorRootSumOfSquares, ImageToImageFilter);


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






protected:
  ReconstructorRootSumOfSquares(){
	  this->SetHasAcceleration(false);
	  this->SetHasSensitivity(false);

  }
  ~ReconstructorRootSumOfSquares(){}




private:
  ReconstructorRootSumOfSquares(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
 // bool m_UseCovarianceMatrix;
  /* Does the real work. */
  virtual void GenerateData();

};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmReconstructorRootSumOfSquares.hxx"
#endif


#endif // __ReconstructorRootSumOfSquares_h
