#ifndef __cmSNRUnitsReconstructorB1Weighted_hxx
#define __cmSNRUnitsReconstructorB1Weighted_hxx


#include "cmSNRUnitsReconstructorB1Weighted.h"
#include "itkImageAlgorithm.h"
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_matrix_inverse.h"


#include "itkVector.h"

namespace cm
{

template< class VectorImageType,class ScalarImageType>
void SNRUnitsReconstructorB1Weighted< VectorImageType,ScalarImageType>
::ThreadedGenerateData(const VectorImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
 {



	VectorImageTypePointer SENSITIVITYMAP=this->GetSensitivityMap();
	VectorImageTypePointer KSPACEDATAIFFT=this->GetInputIFFT();
	ScalarImageTypePointer OUT=this->GetOutput();
	ChannelArrayType iCM=this->GetInverseNoiseCovariance();
		const int Kdimension =KSPACEDATAIFFT->GetNumberOfComponentsPerPixel();
		const int Sdimension =SENSITIVITYMAP->GetNumberOfComponentsPerPixel();




				itk::ImageRegionIterator<VectorImageType> is(SENSITIVITYMAP,outputRegionForThread);
				itk::ImageRegionIterator<VectorImageType> iif(KSPACEDATAIFFT,outputRegionForThread);
				itk::ImageRegionIterator<ScalarImageType> io(OUT,outputRegionForThread);

//				/* prepare */
				is.GoToBegin();
				iif.GoToBegin();
				io.GoToBegin();



				 VectorImagePixelType S;

				 ChannelArrayType SV(Sdimension,1);

				 ChannelArrayType D(1,1);


				 VectorImagePixelType I;
				 ChannelArrayType IV(Kdimension,1);

				 ChannelArrayType OV(1,1);

				 ChannelArrayType REG(1,1);
								 REG(0,0)= (ScalarImagePixelType) std::sqrt(2.0);

				/* BANG!!*/
								while( !iif.IsAtEnd() )
								{
									S=is.Get();
									I=iif.Get();

									/* reset counter and calculate the sos*/

									for(auto y=0;y<Sdimension;y++)
									{
										SV(y,0)=S.GetElement(y);
										IV(y,0)=I.GetElement(y);

									}


									D=SV.conjugate_transpose()*iCM*SV;
									D=D.apply(sqrt);
									OV= SV.conjugate_transpose()*iCM*IV;
								
								
									io.Set((ScalarImagePixelType)   REG(0,0) * std::abs(OV(0,0)) / std::abs(D(0,0)));
									++is;
									++iif;
									++io;
								}



 }



}// end namespace


#endif




//
//
///* Allocate the output */
//
//typedef typename ScalarImageType::PixelType    ScalarImagePixelType;
//typedef typename VectorImageType::PixelType    VectorImagePixelType;
//
///* The first step is to take the signal and the noise data*/
//typename VectorImageType::ConstPointer s = this->GetInput();
//typename VectorImageType::Pointer noise = this->GetNoise();
//
//
///* allocate the Output file (scalar)*/
//typename ScalarImageType::Pointer output = this->GetOutput();
//output->SetBufferedRegion(output->GetRequestedRegion());
//output->Allocate();
//
//
//
//typename VectorImageType::Pointer signal=ConstPointerToPointer<VectorImageType>(s);
//
//vnl_matrix<typename VectorImageType::InternalPixelType>C;
///* Then we calculate the covariance matrix (yes is more a correlation one but the mean of the noise is supposed to be zero)*/
//C=this->calculateCovarianceMatrix();
//int nChan=noise->GetNumberOfComponentsPerPixel();
//
///* then we build the tensor image of the*/
//std::vector<typename ScalarImageType::Pointer> tmp;
//typename VectorImageType::Pointer CIFFT;
//typename ScalarImageType::Pointer tmpIm;
//
//int NP= PixelCount<VectorImageType>(signal);
//
////	typename VectorImageType::IndexType TTT;
////	TTT[0]==128;
////	TTT[1]==64;
//
//for (int t=0;t<nChan; t++){
//	tmpIm=VectorImageElementAsImage<VectorImageType, ScalarImageType>(t,signal);
//	tmpIm=InverseFFT<ScalarImageType>(tmpIm);
//	tmpIm=multiplyImageTimesScalar<ScalarImageType>(tmpIm,sqrt((ScalarImagePixelType)NP) );
//	tmp.push_back(tmpIm);
//};
//
///* compose the */
//CIFFT=composeVectorImage<ScalarImageType,VectorImageType>(tmp);
//
//
//
//
//
//itk::ImageRegionIterator<ScalarImageType> snr(output,output->GetLargestPossibleRegion());
//itk::ImageRegionIterator<VectorImageType> it(CIFFT,CIFFT->GetLargestPossibleRegion());
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



