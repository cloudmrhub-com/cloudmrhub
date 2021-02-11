#ifndef __cmReconstructorWithSensitivityAccelerated_h
#define __cmReconstructorWithSensitivityAccelerated_h

#pragma once
#include "cm.h"
#include "cmReconstructorWithSensitivity.h"

#include "vnl/algo/vnl_matrix_inverse.h"

namespace cm
{

/** \class ReconstructorWithSensitivityAccelerated
 * \brief This class is the parent class for the SNR or image reconstruction which need Sensitivity map.
 * m_SensitivityMapCalculationMode default is inner reference but you can select also BC
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 */


template< class VectorImageType,class ScalarImageType>
class ReconstructorWithSensitivityAccelerated:public ReconstructorWithSensitivity< VectorImageType, ScalarImageType >
{
public:
  /* Standard class typedefs. (similar to using..)*/
  typedef ReconstructorWithSensitivityAccelerated                                                 Self;
  typedef ReconstructorWithSensitivity< VectorImageType, ScalarImageType > Superclass;
  typedef itk::SmartPointer< Self >                                Pointer;
  typedef itk::SmartPointer< const Self >                                ConstPointer;

  /* Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ReconstructorWithSensitivityAccelerated, ReconstructorWithSensitivity);





     	itkSetMacro(AutocalibrationLines,cm::VectorImageType::IndexType);
     	itkGetMacro(AutocalibrationLines,cm::VectorImageType::IndexType);

     	itkSetMacro(Acceleration,cm::VectorImageType::IndexType);
     	itkGetMacro(Acceleration,cm::VectorImageType::IndexType);

     	itkSetMacro(PreMimicAcceleratedSize,cm::VectorImageType::SizeType);
     	itkGetMacro(PreMimicAcceleratedSize,cm::VectorImageType::SizeType);

     	itkSetMacro(MimicAcceleratedAcquisition,bool);
     	itkGetMacro(MimicAcceleratedAcquisition,bool);


protected:
  ReconstructorWithSensitivityAccelerated(){
	  this->SetSensitivityMapCalculationMode(cm::SensitivityMapCalculation::INNER);
	  this->SetSensitivityMap(nullptr);
	  this->SetHasSensitivity(true);


  }
  ~ReconstructorWithSensitivityAccelerated(){}

//  void GetInputIFFT(){};



private:
  ReconstructorWithSensitivityAccelerated(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

cm::VectorImageType::IndexType m_AutocalibrationLines;
cm::VectorImageType::IndexType m_Acceleration;
bool m_MimicAcceleratedAcquisition;
cm::ScalarImageType::SizeType m_PreMimicAcceleratedSize;

};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmReconstructorWithSensitivityAccelerated.hxx"
#endif


#endif // __ReconstructorWithSensitivityAccelerated_h
