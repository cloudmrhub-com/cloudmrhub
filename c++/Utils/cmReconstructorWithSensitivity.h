#ifndef __cmReconstructorWithSensitivity_h
#define __cmReconstructorWithSensitivity_h

#pragma once
#include "cm.h"
#include "cmReconstructor.h"

#include "vnl/algo/vnl_matrix_inverse.h"

namespace cm
{

/** \class ReconstructorWithSensitivity
 * \brief This class is the parent class for the SNR or image reconstruction which need Sensitivity map.
 * m_SensitivityMapCalculationMode default is inner reference but you can select also BC
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 */


template< class VectorImageType,class ScalarImageType>
class ReconstructorWithSensitivity:public Reconstructor< VectorImageType, ScalarImageType >
{
public:
  /* Standard class typedefs. (similar to using..)*/
  typedef ReconstructorWithSensitivity                                                 Self;
  typedef Reconstructor< VectorImageType, ScalarImageType > Superclass;
  typedef itk::SmartPointer< Self >                                Pointer;
  typedef itk::SmartPointer< const Self >                                ConstPointer;

  /* Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ReconstructorWithSensitivity, Reconstructor);


  /*  the main image type*/

     typedef typename VectorImageType::Pointer VectorImageTypePointer;
     typedef typename VectorImageType::PixelType VectorImagePixelType;
     typedef typename VectorImageType::InternalPixelType VectorImageInternalPixelType;


     typedef typename ScalarImageType::Pointer ScalarImagePointerType;
     typedef typename ScalarImageType::PixelType ScalarImagePixelType;
     typedef typename ScalarImageType::InternalPixelType scalarImageInternalPixelType;

     typedef  vnl_matrix<VectorImageInternalPixelType> ChannelArrayType;



     VectorImageTypePointer GetSensitivityMap();


     void ResetSensivityMap(){
    	 this->m_SensitivityMap= nullptr;
     };


     	itkSetMacro(SensitivityMap,VectorImageTypePointer);
		//the get is on the hxx:)

     	itkSetMacro(SensitivityMapSource,VectorImageTypePointer);
     	itkGetMacro(SensitivityMapSource,VectorImageTypePointer);

     	itkSetMacro(BodyCoilSource,VectorImageTypePointer);
     	itkGetMacro(BodyCoilSource,VectorImageTypePointer);

     	itkGetMacro(SensitivityMapCalculationMode,	cm::SensitivityMapCalculation);
     	itkSetMacro(SensitivityMapCalculationMode,	cm::SensitivityMapCalculation);




protected:
  ReconstructorWithSensitivity(){
	  this->m_SensitivityMapCalculationMode=cm::SensitivityMapCalculation::INNER;

	  this->m_SensitivityMap= nullptr;

	  this->SetHasSensitivity(true);


  }
  ~ReconstructorWithSensitivity(){}

  



private:
  ReconstructorWithSensitivity(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

  VectorImageTypePointer m_SensitivityMap;

  VectorImageTypePointer m_SensitivityMapSource;

  VectorImageTypePointer m_BodyCoilSource;

  cm::SensitivityMapCalculation m_SensitivityMapCalculationMode;

  ChannelArrayType m_InverseNoiseCovariance;

};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmReconstructorWithSensitivity.hxx"
#endif


#endif // __ReconstructorWithSensitivity_h
