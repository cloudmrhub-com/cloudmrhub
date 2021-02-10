#ifndef __cmReconstructorWithSensitivityAccelerated_hxx
#define __cmReconstructorWithSensitivityAccelerated_hxx


#include "cmReconstructorWithSensitivityAccelerated.h"
#include "itkImageAlgorithm.h"

#include "cmPhasedArraySensitivityMapFilterInnerReference.h"
#include "cmPhasedArraySensitivityMapFilterBodyCoil.h"
#include "cmPrewhiteningPhasedArray.h"

#include "cm.h"


namespace cm
{

//template< class VectorImageType,class ScalarImageType>
//typename VectorImageType::Pointer ReconstructorWithSensitivityAccelerated< VectorImageType,ScalarImageType>
//::GetSensitivityMap()
// {
//
//	std::cout<<"Start Sensitivity Matrix Reonstruction\n\n";
//
//	typedef cm::PhasedArraySensitivityMapFilterInnerReference<VectorImageType> sensitivityFilterType;
//	typename sensitivityFilterType::Pointer sens=sensitivityFilterType::New();
//	typename sensitivityFilterType::Pointer sens1=sensitivityFilterType::New();
//
//	typedef cm::PhasedArraySensitivityMapFilterBodyCoil<VectorImageType> sensitivityFilterType0;
//	typename sensitivityFilterType0::Pointer sens0=sensitivityFilterType0::New();
//
//	if ( this->m_SensitivityMap == nullptr ){
//
//		switch(this->GetSensitivityMapCalculationMode()){
//		case cm::SensitivityMapCalculation::INNER:
//
//					sens->SetKSpaceDimension(this->GetKSpaceDimension());
//					sens->SetInput(this->GetInput());
//					sens->Update();
//					this->SetSensitivityMap( sens->GetOutput());
//			break;
//		case cm::SensitivityMapCalculation::BODYCOIL:
//
//
//					sens0->SetKSpaceDimension(this->GetKSpaceDimension());
//					sens0->SetInput(this->GetSensitivityMapSource());
//					sens0->SetBodyCoil(this->GetBodyCoilSource());
//
//					sens0->Update();
//
//					this->SetSensitivityMap( sens0->GetOutput());
//			break;
//
//		default:
//			std::cout<<"NO selection has been done for the coil sensitivity matrix computation so i'll use INNER\n\n\n";
//
//			sens1->SetKSpaceDimension(this->GetKSpaceDimension());
//			sens1->SetInput(this->GetSensitivityMapSource());
//			sens1->Update();
//			this->SetSensitivityMap( sens1->GetOutput());
//
//		}
//
//
//		//prewhitning
//		using PrewhiteningFilter=PrewhiteningPhasedArray<VectorImageType,VectorImageType>;
//		typename PrewhiteningFilter::Pointer prewhite =PrewhiteningFilter::New();
//		prewhite->SetInput(this->m_SensitivityMap);
//		prewhite->SetNoiseCovarianceMatrix(this->GetNoiseCovarianceMatrix());
//		prewhite->Update();
//		this->SetSensitivityMap(prewhite->GetOutput());
//
//
//
//	}
//
//
//		return this->m_SensitivityMap;
//
//
// }

}// end namespace


#endif






