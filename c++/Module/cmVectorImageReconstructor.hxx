#ifndef __cmVectorImageReconstructor_hxx
#define __cmVectorImageReconstructor_hxx


#pragma once
#include "cmVectorImageReconstructor.h"
#include "itkImageAlgorithm.h"



#pragma once
#include "utils/utils.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <omp.h>

#pragma once
#include <vnl/vnl_complex_traits.h>

#pragma once
#include "cmiFFTPhasedArrayFilter.h"

#include "itkImageRegionConstIterator.h"

#include <typeinfo>
namespace cm
{

template <class TImage, class TOImage>
 void VectorImageReconstructor<TImage,TOImage>::SetNoiseKSpace(VectorImageTypePointer NoiseKSpace)
{
	this->m_NoiseKSpace=NoiseKSpace;
//	set the noise
typename TImage::SizeType size=NoiseKSpace->GetLargestPossibleRegion().GetSize();

	std::cout << typeid(*NoiseKSpace).name() << '\n';

	typename VectorImageType::Pointer T=NoiseKSpace.GetPointer();

	std::cout << typeid(T).name() << '\n';

	using theIteratorType = itk::ImageRegionConstIteratorWithIndex< VectorImageType >;

	theIteratorType constIterator( T, T->GetLargestPossibleRegion() );

	std::cout<< "NC"<<NoiseKSpace->GetNumberOfComponentsPerPixel()<<std::endl;

    constIterator.GoToBegin();

	long int a=0;
	while(!constIterator.IsAtEnd()){
//		O[a]=constIterator.Get();
		std::cout<<a<< "Index: " << constIterator.GetIndex() <<": " << (typename TImage::PixelType) constIterator.Get()<<std::endl;
		++constIterator;
		a++;

	}



//std::cout<<size<<std::endl;
//std::cout<<"qui2";
	unsigned long long NP=size[0]*size[1];
	short int NC=NoiseKSpace->GetNumberOfComponentsPerPixel();
	ChannelArrayType Noise(NC,getPixelsnumber<TImage>(NoiseKSpace));
	ChannelArrayType sliceNoise;

	using SliceType = itk::VectorImage<typename TImage::InternalPixelType, 2>;
	typename SliceType::Pointer axial;
//	std::cout<<"qui1";
//#pragma omp parallel for shared(NoiseKSpace,Noise)
	for (int slice=0; slice<size[2];slice++){
		//sliced
		axial=axiallySliceThisImage<TImage,SliceType>(NoiseKSpace,slice);
//		vcl_cerr<<" SLICE size"<<axial->GetLargestPossibleRegion().GetSize()<<std::endl;
		sliceNoise=vectorImageToVNLMatrix<SliceType>(axial);
//		vcl_cerr<<" Vector size"<<sliceNoise.size()<<std::endl;
		Noise.update(sliceNoise,0,NP*slice);
	}



//	update the Noise Matrix
	this->m_Noise=Noise;

}

template <class TImage, class TOImage>
typename TImage::InternalPixelType VectorImageReconstructor<TImage,TOImage>::CalculateNoiseBW()
 {
	typename TImage::InternalPixelType NBW=0.0;


	ChannelArrayType noise = this->GetNoise();
	vcl_cerr<<noise;
	std::cout<<noise.rows()<<", "<<noise.cols()<<std::endl;
	std::vector<typename TImage::InternalPixelType> NV;
	for(auto t=0;t<noise.rows();t++)
	{
		NV=VNLVectorToSTDvector<typename TImage::InternalPixelType>(noise.get_row(t));
		NBW+=noiseBand<typename TImage::InternalPixelType>(NV);
	}
	return NBW/ (typename TImage::InternalPixelType)noise.rows();
 }





//template <class TImage,class TOImage>
//vnl_matrix<typename TImage::InternalPixelType> VectorImageReconstructor< TImage,TOImage>
//::GetNoiseCoefficientMatrix()
// {
//
//	vnl_matrix<typename TImage::InternalPixelType> NC;
//	NC=this->GetNoiseCovarianceMatrix();
//	typename vnl_matrix<typename TImage::InternalPixelType>::iterator IT;
//
//	for(IT=NC.begin();IT<NC.end();IT++){
//
//		*IT=*IT/std::sqrt(*IT* std::conj(*IT));
//	}
//
//	return NC;
// }
//
//template <class TImage,class TOImage>
//vnl_matrix<typename TImage::InternalPixelType> VectorImageReconstructor< TImage,TOImage>
//::GetNoiseCovarianceMatrix()
// {
//
//	std::cout<<"start calc noise cov"<<std::endl;
//
//
//
//	/* This Method calculate the Covariance matrix as implemented in Riccardo Lattanzi Code*/
//	ChannelArrayType noise = this->GetNoise();
//
//
//	vnl_matrix<typename TImage::InternalPixelType> COVARIANCEMATRIX(noise.rows(),noise.rows());
//
//
//
//	VectorImageInternalPixelType NB;
//	//std::cout<<"isnan"<<std::isnan(this->GetNoiseBandWidth())<<std::endl;
//	if(this->GetNoiseBandWidthCorrection()){
//		if (std::isnan(this->GetNoiseBandWidth())){
//			NB=this->CalculateNoiseBW();
//		}else{
//			NB=this->GetNoiseBandWidth();
//		}
//
//	}else{
//
//		NB=1.0;
//	}
//
//	std::cout<<"noise bandwith of covariance matrix ="<<NB<<std::endl;
//
//	vnl_vector<typename TImage::InternalPixelType> TMP;
//	vnl_vector<typename TImage::InternalPixelType> TMP2;
//
//	typename TImage::InternalPixelType EL=0.0;
//
//	typename vnl_vector<typename TImage::InternalPixelType>::iterator IT;
//	typename vnl_vector<typename TImage::InternalPixelType>::iterator IT2;
//
//// pay attention to the shared variables:)
////#pragma omp parallel for private(TMP,TMP2)
//	for(int x=0; x<noise.rows();x++){
//		for(int y=0;y<noise.rows();y++){
//
//			TMP=noise.get_row(x);
//			TMP2=noise.get_row(y);
//
//
//			IT2=TMP2.begin();
//			EL=0.0;
////#pragma omp parallel for private(TMP,TMP2)
//			for(IT=TMP.begin();IT<TMP.end();IT++){
//
//				EL+=(*IT * std::conj(*IT2))/(typename TImage::InternalPixelType)(TMP2.size());
//				IT2++;
//
//			};
//
//
//			COVARIANCEMATRIX(x,y)=EL/NB;
//
//
//		}
//	}
//
//
//
//	std::cout<<"end calc cov"<<std::endl;
//
//
//
//	return COVARIANCEMATRIX;
//
//
//
//
//
// }


//template <class TImage,class TOImage>
//typename TImage::Pointer VectorImageReconstructor< TImage,TOImage>
//::GetInputIFFT()
// {
//	 std::cout<<"requested the inverse ifft of the signal"<<std::endl;
//	if ( this->m_InputIFFT == nullptr ){
//
//	 std::cout<<"And i'm calculating"<<std::endl;
//
//					VectorImageTypePointer IF=this->GetSignal();
//
//				typedef cm::iFFTPhasedArrayFilter<VectorImageType> IFFTFilter;
//				typename IFFTFilter::Pointer ifFilter= IFFTFilter::New();//=this->GetIfft();
//
//				ifFilter->SetKSpaceDimension(this->GetKSpaceDimension());
//				ifFilter->SetInput(IF);
//				ifFilter->Update();
//
//				this->SetInputIFFT(ifFilter->GetOutput());
//
//
//	}
//	 std::cout<<"ifft served"<<std::endl;
//
//		return this->m_InputIFFT;
//
//
// }





}// end namespace


#endif
