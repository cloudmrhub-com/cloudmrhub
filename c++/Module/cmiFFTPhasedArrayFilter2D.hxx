#ifndef __iFFTPhasedArrayFilter2D_hxx
#define __iFFTPhasedArrayFilter2D_hxx

#include "cmiFFTPhasedArrayFilter2D.h"




#include "itkImageAlgorithm.h"

#include "itkVectorImage.h"

#include "itkVectorIndexSelectionCastImageFilter.h"

#include "itkComplexToComplexFFTImageFilter.h"
#include "itkFFTShiftImageFilter.h"

#include "itkComposeImageFilter.h"

#include "itkMultiplyImageFilter.h"

namespace cm
{

template< class TImage>
void iFFTPhasedArrayFilter2D< TImage>
::GenerateData()
 {

	/*this is the type that will be used for apply 2dfft */
		typedef itk::Image<InputImageInnerPixelType, 2 >    ScalarImageType;
		/*this is the type that will be used for apply 2dfft */
		typedef typename ScalarImageType::Pointer ScalarImagePointerType;
		ScalarImagePointerType tmpIm = ScalarImageType::New();


	//
		typename TImage::ConstPointer s = this->GetInput();
		/* de reference to a non const pointer the input*/
	//	InputImageTypePointer cs=ConstPointerToPointer<InputImageType>(s);

//		InputImageTypePointer output = this->GetOutput();
//		output->SetBufferedRegion(output->GetRequestedRegion());
//		output->Allocate();


		int nChan=s->GetNumberOfComponentsPerPixel();
		
		
		typedef itk::VectorIndexSelectionCastImageFilter< InputImageType,ScalarImageType> SelectFilterType;

		typedef itk::ComplexToComplexFFTImageFilter<ScalarImageType> invFFTFilterType;
		typename invFFTFilterType::TransformDirectionType T=invFFTFilterType::TransformDirectionType::INVERSE;

		typedef itk::FFTShiftImageFilter< ScalarImageType, ScalarImageType > FFTShiftFilterType;

		typedef itk::MultiplyImageFilter<ScalarImageType,ScalarImageType,ScalarImageType> Multtype;


	typedef itk::ComposeImageFilter<ScalarImageType> ComposeFilterType;
	typename ComposeFilterType::Pointer cFilter = ComposeFilterType::New();
	
	
	//get the number of pixel
	typename InputImageType::SizeType P=s->GetLargestPossibleRegion().GetSize();

	long int NumberOfPixels=1;
	for (auto s=0; s<P.Dimension;s++)
	{
		NumberOfPixels*=P[s];
	}
	
	
	
	
	

		/* then we build the tensor image of the*/
		std::vector<ScalarImagePointerType> tmp;
omp_set_num_threads(omp_get_num_procs());
	#pragma omp parallel for
		for (int t=0;t<nChan; t++){
			
			
			typename SelectFilterType::Pointer sfilter = SelectFilterType::New();
	sfilter->SetInput(s);
	sfilter->SetIndex(t);
	sfilter->Update();
	
	typename invFFTFilterType::Pointer ifilter = invFFTFilterType::New();
	
	ifilter->SetTransformDirection(T);
	ifilter->SetInput(sfilter->GetOutput());

	typename  FFTShiftFilterType::Pointer fftShiftFilter = FFTShiftFilterType::New();
	fftShiftFilter->SetInput( ifilter->GetOutput() );
	fftShiftFilter->Update();
	tmpIm= fftShiftFilter->GetOutput();

			switch(this->GetNormalize()){
			case NormalizeIFFTType::YES:
			
			
		typename Multtype::Pointer mfilter = Multtype::New();

	mfilter->SetInput(tmpIm);
	mfilter->SetConstant(sqrt((InputImageInnerPixelType)NumberOfPixels));
	mfilter->Update();
				tmpIm=mfilter->GetOutput();
			
			break;
			case NormalizeIFFTType::NO:
			//nothing to do:)
			std::cout<<"nothing to do"<<std::endl;
			break;
		}
			cFilter->SetInput(  t,tmpIm);


		};

		/* compose it */
		
		cFilter->Update();
		this->GetOutput()->Graft(cFilter->GetOutput());




 }



 }// end namespace


#endif
































