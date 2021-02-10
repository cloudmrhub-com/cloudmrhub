#ifndef __cmReconstructorWithSensitivity_hxx
#define __cmReconstructorWithSensitivity_hxx


#include "cmReconstructorWithSensitivity.h"
#include "itkImageAlgorithm.h"

#include "cmPhasedArraySensitivityMapFilterInnerReference.h"
#include "cmPhasedArraySensitivityMapFilterBodyCoil.h"
#include "cmPrewhiteningPhasedArray.h"

#include "cm.h"


namespace cm
{

template< class VectorImageType,class ScalarImageType>
typename VectorImageType::Pointer ReconstructorWithSensitivity< VectorImageType,ScalarImageType>
::GetSensitivityMap()
 {

	std::cout<<"Start Sensitivity Matrix Reonstruction\n\n";

	typedef cm::PhasedArraySensitivityMapFilterInnerReference<VectorImageType> sensitivityFilterType;
	typename sensitivityFilterType::Pointer sens=sensitivityFilterType::New();
	typename sensitivityFilterType::Pointer sens1=sensitivityFilterType::New();

	typedef cm::PhasedArraySensitivityMapFilterBodyCoil<VectorImageType> sensitivityFilterType0;
	typename sensitivityFilterType0::Pointer sens0=sensitivityFilterType0::New();

	if ( this->m_SensitivityMap == nullptr ){

		switch(this->GetSensitivityMapCalculationMode()){
		case cm::SensitivityMapCalculation::INNER:

					sens->SetKSpaceDimension(this->GetKSpaceDimension());
					sens->SetInput(this->GetInput());
					sens->Update();
					this->SetSensitivityMap( sens->GetOutput());
			break;
		case cm::SensitivityMapCalculation::BODYCOIL:


					sens0->SetKSpaceDimension(this->GetKSpaceDimension());
					sens0->SetInput(this->GetSensitivityMapSource());
					sens0->SetBodyCoil(this->GetBodyCoilSource());

					sens0->Update();

					this->SetSensitivityMap( sens0->GetOutput());
			break;

		default:
			std::cout<<"NO selection has been done for the coil sensitivity matrix computation so i'll use INNER\n\n\n";

			sens1->SetKSpaceDimension(this->GetKSpaceDimension());
			sens1->SetInput(this->GetSensitivityMapSource());
			sens1->Update();
			this->SetSensitivityMap( sens1->GetOutput());

		}


		//prewhitning
		using PrewhiteningFilter=PrewhiteningPhasedArray<VectorImageType,VectorImageType>;
		typename PrewhiteningFilter::Pointer prewhite =PrewhiteningFilter::New();
		prewhite->SetInput(this->m_SensitivityMap);
		prewhite->SetNoiseCovarianceMatrix(this->GetNoiseCovarianceMatrix());
		prewhite->Update();
		this->SetSensitivityMap(prewhite->GetOutput());



	}
	
	
		return this->m_SensitivityMap;
	

 }

}// end namespace


#endif




//
//
///* Allocate the output */
//
//typedef typename TOImage::PixelType    ScalarImagePixelType;
//typedef typename TImage::PixelType    VectorImagePixelType;
//
///* The first step is to take the signal and the noise data*/
//typename TImage::ConstPointer s = this->GetInput();
//typename TImage::Pointer noise = this->GetNoise();
//
//
///* allocate the Output file (scalar)*/
//typename TOImage::Pointer output = this->GetOutput();
//output->SetBufferedRegion(output->GetRequestedRegion());
//output->Allocate();
//
//
//
//typename TImage::Pointer signal=ConstPointerToPointer<TImage>(s);
//
//vnl_matrix<typename TImage::InternalPixelType>C;
///* Then we calculate the covariance matrix (yes is more a correlation one but the mean of the noise is supposed to be zero)*/
//C=this->calculateCovarianceMatrix();
//int nChan=noise->GetNumberOfComponentsPerPixel();
//
///* then we build the tensor image of the*/
//std::vector<typename TOImage::Pointer> tmp;
//typename TImage::Pointer CIFFT;
//typename TOImage::Pointer tmpIm;
//
//int NP= PixelCount<TImage>(signal);
//
////	typename TImage::IndexType TTT;
////	TTT[0]==128;
////	TTT[1]==64;
//
//for (int t=0;t<nChan; t++){
//	tmpIm=VectorImageElementAsImage<TImage, TOImage>(t,signal);
//	tmpIm=InverseFFT<TOImage>(tmpIm);
//	tmpIm=multiplyImageTimesScalar<TOImage>(tmpIm,sqrt((ScalarImagePixelType)NP) );
//	tmp.push_back(tmpIm);
//};
//
///* compose the */
//CIFFT=composeVectorImage<TOImage,TImage>(tmp);
//
//
//
//
//
//itk::ImageRegionIterator<TOImage> snr(output,output->GetLargestPossibleRegion());
//itk::ImageRegionIterator<TImage> it(CIFFT,CIFFT->GetLargestPossibleRegion());
//
//it.GoToBegin();
//snr.GoToBegin();
//
//
//VectorImagePixelType v;
//VectorImagePixelType v2;
////	ScalarImagePixelType sMag;
////	ScalarImagePixelType nPow;
//
//
//
//
//
//vnl_vector<ScalarImagePixelType>S(nChan);
//vnl_matrix<ScalarImagePixelType> X(nChan,1);
//vnl_matrix<ScalarImagePixelType> sMag;
//vnl_matrix<ScalarImagePixelType> nPow;
//
//
//
//
//
//
//ScalarImagePixelType N= sqrt(2);
//
//
//
//
//typename vnl_matrix<ScalarImagePixelType>::iterator II;
//
//
//
//if (this->GetUseCovarianceMatrix())
//{
//	while( !snr.IsAtEnd() )
//	{
//		v=it.Get();
//
//		for (int t=0;t<v.GetSize();t++)
//		{
//
//			X(t,0)=v.GetElement(t);
//
//		}
//
//
////			std::cout<<"X is";
////			vcl_cout<<X;
////			std::cout<<std::endl;
//
//
//		sMag=(X.conjugate_transpose()*X);
//
////			sMag=(X.transpose()*X);
//
//
//
//		for (II= sMag.begin(); II != sMag.end(); II++){
//			*II=abs(*II)*N;
//		};
//
//
//
//		nPow= X.transpose()*C*X;
//
//
////			vcl_cout<<nPow;
////			std::cout<<std::endl;
//		for (II= nPow.begin(); II != nPow.end(); II++){
//			*II=sqrt(abs(*II));
//		};
//
////			std::cout<<"after sqrt abs";
////			vcl_cout<<nPow;
////			std::cout<<std::endl;
//
//	//std::cout<<"SNR"<<sMag(0,0)/nPow(0,0)<<std::endl;
//		snr.Set((ScalarImagePixelType)(sMag(0,0)/nPow(0,0)));
//		++it;
//		++snr;
//	}
//
//}else{
//
//	try
//				{
//				sMag=vnl_matrix_inverse<ScalarImagePixelType>(C);
//				}catch (const std::exception& e) { std::cout<<"NO"; }
//
//
//
//	while( !snr.IsAtEnd() )
//	{
//
//		v=it.Get();
//
//		for (int t=0;t<v.GetSize();t++)
//		{
//
//			X(t,0)=v.GetElement(t);
//
//		}
//
//
//		//nPow=X.transpose()*sMag*X;
//
//		nPow=X.conjugate_transpose()*sMag*X;
//
//		snr.Set(N*(ScalarImagePixelType)sqrt(nPow(0,0)));
//
//		++it;
//		++snr;
//	}
//
//
//}
//
//
//
//
//
//}



