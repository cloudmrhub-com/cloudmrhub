#ifndef __cmReconstructor_hxx
#define __cmReconstructor_hxx


#pragma once
#include "cmReconstructor.h"
#include "itkImageAlgorithm.h"




#include "utils.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <omp.h>

#include <vnl/vnl_complex_traits.h>

#include "cmiFFTPhasedArrayFilter.h"

#include "ImageUtils.h"

#include "cmPrewhiteningPhasedArray.h"

namespace cm
{


template< class TImage,class TOImage>
ChannelArrayType Reconstructor< TImage,TOImage>
::GetNoiseCovarianceMatrix(){

	if (this->m_NoiseCovarianceMatrix.size()==1){

		this->SetNoiseCovarianceMatrix(this->CalculateNoiseCovarianceMatrix());

	}else{
		return this->m_NoiseCovarianceMatrix;
	}

	return this->m_NoiseCovarianceMatrix;


};

template< class TImage,class TOImage>
void Reconstructor< TImage,TOImage>
::SetNoiseKSpace(VectorImageTypePointer NoiseKSpace)
{
	this->m_NoiseKSpace=NoiseKSpace;
//	set the noise
typename TImage::SizeType size=NoiseKSpace->GetLargestPossibleRegion().GetSize();


std::cout<<size<<"KSPACE SIZE";

	unsigned long long NP=size[0]*size[1];
	short int NC=NoiseKSpace->GetNumberOfComponentsPerPixel();

	ChannelArrayType Noise(NC,NP * size[2],0.0);

	ChannelArrayType sliceNoise;

	using SliceType = itk::VectorImage<typename TImage::InternalPixelType, 2>;
	typename SliceType::Pointer axial;

#pragma omp parallel for shared(Noise)
	for (int slice=0; slice<size[2];slice++){
		axial=axiallySliceThisImage<TImage,SliceType>(NoiseKSpace,slice);

		sliceNoise=vectorImageToVNLMatrix<SliceType>(axial);
		Noise.update(sliceNoise,0,NP*slice);
	}


//	update the Noise Matrix
	this->m_Noise=Noise;


	ScalarImagePixelType BW=this->CalculateNoiseBW();
	this->SetNoiseBandWidth( BW.real());


	this->SetNoiseCovarianceMatrix(this->CalculateNoiseCovarianceMatrix());

}

template< class TImage,class TOImage>
typename TImage::InternalPixelType Reconstructor< TImage,TOImage>
::CalculateNoiseBW()
 {
	typename TImage::InternalPixelType NBW=0.0;

	ChannelArrayType noise = this->GetNoise();
	std::vector<typename TImage::InternalPixelType> NV;
	for(auto t=0;t<noise.rows();t++)
	{


		NV=VNLVectorToSTDvector<typename TImage::InternalPixelType>(noise.get_row(t));
		NBW+=noiseBand<typename TImage::InternalPixelType>(NV);
	}
	return NBW/ (typename TImage::InternalPixelType)noise.rows();
 }





template< class TImage,class TOImage>
vnl_matrix<typename TImage::InternalPixelType> Reconstructor< TImage,TOImage>
::CalculateNoiseCoefficientMatrix()
 {

	vnl_matrix<typename TImage::InternalPixelType> NC;
	NC=this->GetNoiseCovarianceMatrix();
	typename vnl_matrix<typename TImage::InternalPixelType>::iterator IT;

	for(IT=NC.begin();IT<NC.end();IT++){

		*IT=*IT/std::sqrt(*IT* std::conj(*IT));
	}

	return NC;
 }

//
template< class TImage,class TOImage>
vnl_matrix<typename TImage::InternalPixelType> Reconstructor< TImage,TOImage>
::CalculateNoiseCovarianceMatrix()
 {



	/* This Method calculate the Covariance matrix as implemented in Riccardo Lattanzi Code*/
    //update 09/24/2020 Riccardo Lattanzi
	ChannelArrayType noise = this->GetNoise();
	vnl_matrix<typename TImage::InternalPixelType> COVARIANCEMATRIX(noise.rows(),noise.rows());



	VectorImageInternalPixelType NB;
	//std::cout<<"isnan"<<std::isnan(this->GetNoiseBandWidth())<<std::endl;
	if(this->GetNoiseBandWidthCorrection()){
		if (std::isnan(this->GetNoiseBandWidth())){
			NB=this->CalculateNoiseBW();
		}else{
			NB=this->GetNoiseBandWidth();
		}

	}else{

		NB=1.0;
	}



	vnl_matrix<typename TImage::InternalPixelType> TMP(1,noise.cols());
	vnl_matrix<typename TImage::InternalPixelType> TMP2(1,noise.cols());

	vnl_matrix<typename TImage::InternalPixelType> EL(1,1);
	typename TImage::InternalPixelType D= 2 * noise.cols();
	for(int x=0; x<noise.rows();x++){
		for(int y=0;y<noise.rows();y++){
			TMP.set_row(0,noise.get_row(x));
			TMP2.set_row(0,noise.get_row(y));
			EL=(TMP * TMP2.conjugate_transpose())/D;
			COVARIANCEMATRIX(x,y)=EL[0][0]/NB;
		}
	}






	return COVARIANCEMATRIX;





 }


template< class TImage,class TOImage>
typename TImage::Pointer Reconstructor< TImage,TOImage>
::GetInputIFFT()
 {
	 std::cout<<"requested the inverse ifft of the signal"<<std::endl;
	if ( this->m_InputIFFT == nullptr ){


	 	 	 	 VectorImageTypePointer IF=this->GetInput();

	 	 	 	std::cout<<"Prewhitening of the signal KSpace"<<std::endl;
	 	 	 	using PrewhiteningFilter=PrewhiteningPhasedArray<VectorImageType,VectorImageType>;
	 	 	 	typename PrewhiteningFilter::Pointer prewhite =PrewhiteningFilter::New();
	 	 	 	prewhite->SetInput(IF);
	 	 	 	prewhite->SetNoiseCovarianceMatrix(this->GetNoiseCovarianceMatrix());
	 	 	 	prewhite->Update();




				typedef cm::iFFTPhasedArrayFilter<VectorImageType> IFFTFilter;
				typename IFFTFilter::Pointer ifFilter= IFFTFilter::New();//=this->GetIfft();

				ifFilter->SetKSpaceDimension(this->GetKSpaceDimension());
				ifFilter->SetInput(prewhite->GetOutput());
				ifFilter->Update();
				this->SetInputIFFT(ifFilter->GetOutput());


	}
	 std::cout<<"ifft served"<<std::endl;

		return this->m_InputIFFT;


 }





}// end namespace


#endif
