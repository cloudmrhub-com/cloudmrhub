//#pragme once
//#include <vnl/vnl_complex_traits.h>
//
//template <class ValueType>

#pragma once
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include <iostream>
#include <vector>
#include <functional>
#include <future>
#include <functional>
#include <type_traits>
#include "itkCastImageFilter.h"


//OR
#include "itkOrImageFilter.h"

//MINMAX
#include "itkMinimumMaximumImageCalculator.h"
//SUBTRACT,SNR
#include "itkSubtractImageFilter.h"

// multiply
#include "itkMultiplyImageFilter.h"



#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToPhaseImageFilter.h"

#include "itkImage.h"

#include "itkVector.h"
#include "itkVariableLengthVector.h"

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>



#include "itkComplexToComplexFFTImageFilter.h"
#include "itkFFTShiftImageFilter.h"

#pragma once
#include "itkShrinkImageFilter.h"
template <class  VectorImageType>
typename VectorImageType::Pointer shrinkImage(typename VectorImageType::Pointer image, std::vector<int>& ShrinkFactor)
{

	/* Let's shrink the image the new size will be size./accelerations*/
	typedef itk::ShrinkImageFilter <VectorImageType,VectorImageType> 	ShrinkImageFilterType;
	typename ShrinkImageFilterType::Pointer shrinkFilter = ShrinkImageFilterType::New();
	shrinkFilter->SetInput(image);

	for (auto it= 0; it <ShrinkFactor.size(); it++){
		shrinkFilter->SetShrinkFactor(it, ShrinkFactor.at(it)); // shrink the first dimension by a factor of 2 (i.e. 100 gets changed to 50)
	};


	shrinkFilter->Update();

	return shrinkFilter->GetOutput();
}

template <class  SourceImageType,class  ReferenceImageType>
void copyReferenceInfoonImage(typename SourceImageType::Pointer image,typename ReferenceImageType::Pointer ref)
{

	image->SetOrigin(ref->GetOrigin());
	image->SetDirection(ref->GetDirection());
	image->SetSpacing(ref->GetSpacing());
}

#include "itkImageRegionIterator.h"

template <class TImage>
long int getPixelsnumber( typename TImage::Pointer image)
{
	typename TImage::SizeType P=image->GetLargestPossibleRegion().GetSize();

	long int m=1;
	for (auto s=0; s<P.Dimension;s++)
	{
		m *= P[s];
	}


	return m;
}



//template <class TImage>
//vnl_vector<typename TImage::InternalPixelType> imageToVNLVector( typename TImage::Pointer image)
//{
//	vnl_vector<typename TImage::InternalPixelType> O(getPixelsnumber<TImage>(image));
//
//	itk::ImageRegionIterator<TImage> imageIterator(image,image->GetLargestPossibleRegion());
//
//	long int a=0;
//	while(!imageIterator.IsAtEnd()){
//		O[a]=imageIterator.Get();
//		++imageIterator;
//		a++;
//
//	}
//
//
//	return O;
//}
#include "vnl/vnl_c_vector.h"

template <class TImage>
vnl_vector<typename TImage::InternalPixelType> KSpace2DToVNLVector( typename TImage::Pointer image)
{
	vnl_vector<typename TImage::InternalPixelType> O(getPixelsnumber<TImage>(image));

	typename TImage::SizeType size=image->GetLargestPossibleRegion().GetSize();

	long long a=0;

	typename TImage::IndexType PT;

	for (long long y=0; y<size[1];y++){ //for each phase encoding
		for (long long x=0; x<size[0];x++){

			PT[0]=x;
			PT[1]=y;
			PT[2]=0;
			O[a]=image->GetPixel(PT);
			a++;
		}

	};
//	std::cout<<"###############################################################3\n";
//	vcl_cerr<<"O"<<O<<vcl_endl;
//	std::cout<<"###############################################################3\n";


	return O;
}



#pragma once
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "vnl/vnl_matrix.h"
#include "itkImage.h"
template <class TImage>
vnl_matrix<typename TImage::InternalPixelType> vectorImageToVNLMatrix( typename TImage::Pointer image)
{
	short int NC=image->GetNumberOfComponentsPerPixel();
	vnl_matrix<typename TImage::InternalPixelType> O(NC,getPixelsnumber<TImage>(image));


	typedef itk::Image<typename TImage::InternalPixelType, TImage::ImageDimension> ScalarImageType;


	typedef itk::VectorIndexSelectionCastImageFilter< TImage,ScalarImageType> FilterType;

	for (auto nc=0;nc<NC;nc++){
		typename FilterType::Pointer filter = FilterType::New();
		filter->SetInput(image);
		filter->SetIndex(nc);
		filter->Update();


		O.set_row(nc,KSpace2DToVNLVector<ScalarImageType>(filter->GetOutput()));
	}

	return O;
}


#pragma once
#include "vnl/algo/vnl_fft_1d.h"
#include "vcl_complex.h"
template <class Pxl>
Pxl noiseBand(std::vector<Pxl> & X){
	/** from Riccardo lattanzi mrir_noise_bandwidth */

	//https://itk.org/Doxygen/html/Examples_2Numerics_2FourierDescriptors1_8cxx-example.html
	//https://github.com/InsightSoftwareConsortium/ITK/blob/master/Examples/Numerics/FourierDescriptors1.cxx





	/** basically powerspectrum of thee signal abs(fft(noise)).^2
	 * power_spectrum_norm = mean(power_spectrum([1:dim/4 3*dim/4:end]));
	 *norm_power_spectrum = power_spectrum./repmat(power_spectrum_norm,1, dim);
	 *noise_bandwidth = mean(norm_power_spectrum);
	 * */

	const auto powerOfTwo   = (unsigned int)std::ceil(
			std::log( (double)(X.size())) /
			std::log( (double)(2.0)) );

	const unsigned int spectrumSize = 1 << powerOfTwo;


	vnl_fft_1d<float> fft( spectrumSize);


	std::vector<Pxl>signal( spectrumSize );





	for (auto it=0; it < X.size(); it++){
		signal[it]=X.at(it);
	};

	/** zero padding*/
	for(unsigned int pad=X.size(); pad<spectrumSize; pad++)
	{
		signal[pad] = 0.0;
	}

	/** actually copute transform*/
	fft.fwd_transform(signal);

	typename std::vector< Pxl>::iterator OO;



	Pxl Su=(Pxl) 0.0;
	Pxl TT;

	/**calculate power_spectrum_norm  (Su) and accumulate the spectrum*/
	unsigned int LL=spectrumSize/4;
	unsigned int UL=3*spectrumSize/4;

	std::vector<Pxl>spectrum;
	for (OO= signal.begin(); OO != signal.end(); OO++){
		//		spectrum.at([signal.begin(),OO])=abs(*OO)*abs(*OO);

		//		std::cout<<std::distance(signal.begin(),OO);

		TT=std::pow(abs(*OO),2);

		if(std::distance(signal.begin(), OO)<LL  || std::distance(signal.begin(), OO)> UL )
		{
			Su+= TT;
			//			std::cout <<"yes";
		}

		spectrum.push_back(TT);
		//		std::cout<<std::endl;
	};

	//std::cout <<Su<<std::endl;
	Su/=(spectrumSize/2);
	//std::cout<< Su<<std::endl;

	/** normalize and comuputemeans*/
	Pxl BW=(Pxl)0.0;
	for (OO= spectrum.begin(); OO != spectrum.end(); OO++){
		BW+=(*OO/Su)/(Pxl)spectrumSize;
	};








	//std::cout<< "BANDWIDTH IS:"<<BW<< std::endl;

	return (Pxl) BW;

};



template <class Pxl>
std::vector<Pxl> VNLVectorToSTDvector(vnl_vector<Pxl>  X){

	std::vector<Pxl>V(X.size());
	for (auto it=0; it < X.size(); it++){
		V.at(it)=X[it];
	};
	return V;
}



#pragma once
#include "itkImageFileReader.h"

template <class TImage>
typename TImage::Pointer readImage(std::string name)
{
	typedef itk::ImageFileReader<TImage >  RType;
	typename RType::Pointer r = RType::New();
	r->SetFileName(name);
	r->Update();
	std::cout << name << " read" << std::endl;
	return r->GetOutput();
}


#pragma once
#include "itkImageFileWriter.h"
template <class TImage>
void writeImage(typename TImage::Pointer image, std::string name)
{
	typedef itk::ImageFileWriter<TImage >  WType;
	typename WType::Pointer w = WType::New();
	w->SetInput(image);
	w->SetFileName(name);
	w->Update();


	std::cout << name <<" written" << std::endl;
}



template <class from, class to>
vnl_matrix<to> cast2DVNLMatrixfromto(vnl_matrix<from>& f){
	vnl_matrix<to>Out(f.rows(),f.cols());
	for (auto a=0;a<f.rows();a++){

		for (auto b=0;b<f.cols();b++){
			Out(a,b)=(to)f(a,b);
		}
	}
	return Out;

};


#pragma once
#include "itkExtractImageFilter.h"
template <class From, class To>
typename To::Pointer axiallySliceThisImage(typename From::Pointer image,int slice){

	typename From::RegionType inputRegion =image->GetLargestPossibleRegion();
	typename From::SizeType size = inputRegion.GetSize();
	size[2] = 0; //collpase the third dimension

	typename From::IndexType 		 start = inputRegion.GetIndex();
	start[2] = slice;

	typename From::RegionType desiredRegion;

	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);


	using ExtractFilterType = itk::ExtractImageFilter< From,To >;
	typename ExtractFilterType::Pointer filter = ExtractFilterType::New();
	filter->SetInput(image);
	filter->InPlaceOn();
	filter->SetDirectionCollapseToSubmatrix();
	filter->SetExtractionRegion(desiredRegion);
	filter->Update();
	return filter->GetOutput();
}


template< class TImage>
typename TImage::Pointer ForwardFFT(typename TImage::Pointer image)
{

	typedef itk::ComplexToComplexFFTImageFilter<TImage> invFFTFilterType;
	typename invFFTFilterType::Pointer filter = invFFTFilterType::New();
	typename invFFTFilterType::TransformDirectionType T=invFFTFilterType::TransformDirectionType::FORWARD;
	filter->SetTransformDirection(T);
	filter->SetInput(image);


	typedef itk::FFTShiftImageFilter< typename invFFTFilterType::OutputImageType, typename invFFTFilterType::OutputImageType > FFTShiftFilterType;

	typename  FFTShiftFilterType::Pointer fftShiftFilter = FFTShiftFilterType::New();
	fftShiftFilter->SetInput( filter->GetOutput() );
	fftShiftFilter->Update();

	return fftShiftFilter->GetOutput();


}

template< class TImage>
typename TImage::Pointer InverseFFT(typename TImage::Pointer image)
{

	typedef itk::ComplexToComplexFFTImageFilter<TImage> invFFTFilterType;
	typename invFFTFilterType::Pointer filter = invFFTFilterType::New();
	typename invFFTFilterType::TransformDirectionType T=invFFTFilterType::TransformDirectionType::INVERSE;
	filter->SetTransformDirection(T);
	filter->SetInput(image);


	typedef itk::FFTShiftImageFilter< typename invFFTFilterType::OutputImageType, typename invFFTFilterType::OutputImageType > FFTShiftFilterType;

	typename  FFTShiftFilterType::Pointer fftShiftFilter = FFTShiftFilterType::New();
	fftShiftFilter->SetInput( filter->GetOutput() );
	fftShiftFilter->Update();



	return fftShiftFilter->GetOutput();


}

//template <class TImage>
//long int PixelCountOld( typename TImage::Pointer image)
//{
//	std::vector<typename TImage::PixelType>p=imageToArray<TImage>(image);
//
//	return p.size();
//}

#include "itkImageDuplicator.h"
template< class TImage>
typename TImage::Pointer ConstPointerToPointer(typename TImage::ConstPointer image)
{

	  typedef itk::ImageDuplicator< TImage > DuplicatorType;
		  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
		  duplicator->SetInputImage(image);
		  duplicator->Update();



	return duplicator->GetOutput();


}

#include "itkVectorImage.h"

#include "itkVectorIndexSelectionCastImageFilter.h"
/**this filter give an image of type TOut from a versor image */
template< class TVImage, class TSImage>
typename TSImage::Pointer VectorImageElementAsImage(unsigned int x, typename TVImage::Pointer image)
{



	typedef itk::VectorIndexSelectionCastImageFilter< TVImage,TSImage> FilterType;
	typename FilterType::Pointer filter = FilterType::New();
	filter->SetInput(image);
	filter->SetIndex(x);
	filter->Update();
	return filter->GetOutput();

}


template <class TImage>
std::vector<typename TImage::PixelType> imageToArray( typename TImage::Pointer image)
{
	std::vector<typename TImage::PixelType> O;

	itk::ImageRegionIterator<TImage> imageIterator(image,image->GetLargestPossibleRegion());

	while(!imageIterator.IsAtEnd()){
		O.push_back(imageIterator.Get());
		++imageIterator;

	}


	return O;
}



template<class VP, class TImage>
// VectorImageElementAsArray(int x, typename TImage::Pointer image)
std::vector<VP> VectorImageElementAsArray(int x, typename TImage::Pointer image)
{


	/** Create the new filetype*/
	typedef typename itk::Image<VP, TImage::ImageDimension >    ScalarImageType;
	/** get the image*/
	typename ScalarImageType::Pointer im = ScalarImageType::New();

	im=VectorImageElementAsImage<TImage,ScalarImageType>(x,image);


	std::vector<VP> o=imageToArray<ScalarImageType>(im);

	return o;

}

template <class TS,class TV>
typename TV::Pointer addImageAsDimensionToVectorImage(typename TS::Pointer s, typename TV::Pointer v)
{

	int NCV=v->GetNumberOfComponentsPerPixel();

	typedef itk::ComposeImageFilter<TS> FilterType;
	typename FilterType::Pointer Filter = FilterType::New();

	/** we want to decompose the image and add the last image in the queue*/
	for(int nim=0;nim <NCV; nim++){



		Filter->SetInput( nim ,VectorImageElementAsImage<TV,TS>(nim, v));
	};

	Filter->SetInput(NCV,s);

	Filter->Update();
	//obtain the mask in the image FOV
	return Filter->GetOutput();



}




template< class TImage>
typename TImage::Pointer fftShift(typename TImage::Pointer image)
{
	typedef itk::FFTShiftImageFilter< TImage,TImage > FFTShiftFilterType;

	typename  FFTShiftFilterType::Pointer fftShiftFilter = FFTShiftFilterType::New();
	fftShiftFilter->SetInput( image);
	fftShiftFilter->Update();

	return fftShiftFilter->GetOutput();
}







































template <class Tin,class Tout>
typename Tout::Pointer castImage(typename Tin::Pointer im)
{

	typedef itk::CastImageFilter< Tin, Tout > CastFilterType;
	typename CastFilterType::Pointer castFilter = CastFilterType::New();

	castFilter->SetInput(im);
	castFilter->Update();
	//obtain the mask in the image FOV
	return castFilter->GetOutput();

	std::cout << "image casted" << std::endl;


}






template <class TImage, class TRoi>
typename TRoi::Pointer roiOnImage(typename TImage::Pointer image, typename TRoi::Pointer roi)
{
	typedef itk::IdentityTransform< double, TImage::ImageDimension >  TransformType;
	typename TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();

	typedef itk::NearestNeighborInterpolateImageFunction< TRoi, double > InterpolatorType;
	typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

	typedef itk::ResampleImageFilter< TRoi, TRoi > ResampleFilterType;
	typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();

	resampleFilter->SetInput(roi);
	resampleFilter->SetTransform(transform);
	resampleFilter->SetInterpolator(interpolator);
	resampleFilter->SetOutputParametersFromImage(image);
	resampleFilter->Update();


	//writeImage<TRoi>(resampleFilter->GetOutput(), "p.mha");
	std::cout << "Roi Resampled" << std::endl;

	//obtain the mask in the image FOV
	return resampleFilter->GetOutput();


}


template<class TImage, class TImage_ >
typename TImage::Pointer orImages(typename TImage::Pointer r1, typename TImage_::Pointer r2) {

	//report image n the same size non so se serve but:)
	typename TImage::Pointer ROI2on1 = TImage::New();
	ROI2on1 = roiOnImage<TImage, TImage_>(r1, r2);

	typedef itk::OrImageFilter <TImage> OrImageFilterType;
	typename OrImageFilterType::Pointer orFilter = OrImageFilterType::New();

	orFilter->SetInput(0, r1);
	orFilter->SetInput(1, ROI2on1);
	orFilter->Update();
	std::cout << " OR--" << std::endl;

	return orFilter->GetOutput();
}


template<class TImage>
typename TImage::Pointer subtractImage(typename TImage::Pointer r1, typename TImage::Pointer r2) {


	typedef itk::SubtractImageFilter<TImage,TImage> SIFtype;
	typename SIFtype::Pointer filter = SIFtype::New();

	filter->SetInput1(r1);
	filter->SetInput2(r2);
	filter->Update();
	std::cout << " minus--" << std::endl;

	return filter->GetOutput();
}

template<class TImage>
typename TImage::Pointer multiplyImageTimesScalar(typename TImage::Pointer im, typename TImage::PixelType c) {


	typedef itk::MultiplyImageFilter<TImage,TImage,TImage> Ftype;
	typename Ftype::Pointer filter = Ftype::New();

	filter->SetInput(im);
	filter->SetConstant(c);
	filter->Update();

	return filter->GetOutput();
}



#include "itkComposeImageFilter.h"
template <class TS,class TV>
typename TV::Pointer composeVectorImage(std::vector<typename TS::Pointer> & vect)
{

	typedef itk::ComposeImageFilter<TS> FilterType;
	typename FilterType::Pointer Filter = FilterType::New();
	for(auto nim=vect.begin();nim != vect.end(); ++nim){
		Filter->SetInput( std::distance(vect.begin(), nim) ,*nim);
	};

	Filter->Update();


	std::cout << "images correctly paked" << std::endl;
	//obtain the mask in the image FOV
	return Filter->GetOutput();


}





template <class Pxl>
vnl_vector<Pxl> stdVectorToVnl(std::vector<Pxl> & X){

	vnl_vector<Pxl> x(X.size());

	typename std::vector<Pxl>::iterator it;


	for (it= X.begin(); it != X.end(); it++){
		x(std::distance(X.begin(), it))=*it;
	};


	return x;

};






