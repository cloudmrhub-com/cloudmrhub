#ifndef __cmiFFTPhasedArrayFilter_hxx
#define __cmiFFTPhasedArrayFilter_hxx

#include "itkVectorImage.h"
#include "cmiFFTPhasedArrayFilter.h"
#include "itkImageAlgorithm.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "cmiFFTPhasedArrayFilter2D.h"

#include "itkExtractImageFilter.h"

#include "itkTileImageFilter.h"

namespace cm
{



template< class VectorImageType>
void iFFTPhasedArrayFilter< VectorImageType>::ifft2D()
 {
/** The function 
 * sliced the 3D kspace in the third dimension*/

std::cout << "start recontruction as sum of squares 2D";

	
	typename VectorImageType::ConstPointer inputImage = this->GetInput();


//std::cout<<inputImage;
//FULL=tFilter->GetOutput();
//FULL->SetOrigin(r->GetOutput()->GetOrigin());
//FULL->SetDirection(r->GetOutput()->GetDirection());
//FULL->SetSpacing(r->GetOutput()->GetSpacing());




/** and define the *SliceImageType*/
	typedef typename itk::VectorImage<typename VectorImageType::InternalPixelType, 2 >  VectorSliceImageType;
	typedef typename itk::Image<typename VectorImageType::InternalPixelType, 2 >    ScalarSliceImageType;

	
	typename VectorImageType::IndexType test3;
	test3[0]=10;
	test3[1]=10;
	test3[2]=0;
	
	
	typename VectorImageType::RegionType inputRegion = inputImage->GetLargestPossibleRegion();
	typename VectorImageType::SizeType size = inputRegion.GetSize();
	
	/** then for each slice */
	int sl= size[2];
	size[2] = 0; // we extract along z direction
	typename VectorImageType::IndexType start = inputRegion.GetIndex();



	typename VectorSliceImageType::SizeType sizeslice;
		sizeslice[0]=size[0];
		sizeslice[1]=size[1];



	typename VectorSliceImageType::IndexType startslice;
	startslice[0]=start[0];
	startslice[1]=start[1];


	using ExtractFilterType = itk::ExtractImageFilter< VectorImageType, VectorSliceImageType >;


	typename VectorImageType::Pointer output = this->GetOutput();
	output->SetBufferedRegion(output->GetRequestedRegion());
	output->Allocate();

//	omp_set_num_threads(80);
#pragma omp parallel for
	for (auto t=0; t<sl;t++)
	{

		
		typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
		extractFilter->SetDirectionCollapseToGuess();

		const unsigned int sliceNumber =t;
		start[2] = sliceNumber;
		typename VectorImageType::RegionType desiredRegion;
		desiredRegion.SetSize(  size  );
		desiredRegion.SetIndex( start );


		typename VectorSliceImageType::RegionType sliceddesiredRegion;
				sliceddesiredRegion.SetSize(  sizeslice  );
				sliceddesiredRegion.SetIndex( startslice );


		extractFilter->SetExtractionRegion( desiredRegion );
		extractFilter->InPlaceOn();
		extractFilter->SetInput( inputImage );
		extractFilter->Update();

	

/** we have just extracted the 2D KSpace*/


/** we make the ifft 2D*/
		typename VectorSliceImageType::Pointer CIFFT;

		using ifftType= cm::iFFTPhasedArrayFilter2D<VectorSliceImageType>;

		typename ifftType::Pointer ifft=ifftType::New();
		ifft->SetInput(extractFilter->GetOutput());
		ifft->Update();
		CIFFT=ifft->GetOutput();


	
				/* instantiate the runners!!*/
		itk::ImageRegionConstIterator<VectorSliceImageType> it(CIFFT,sliceddesiredRegion);
		itk::ImageRegionIterator<VectorImageType> ot(output,desiredRegion);

		/* prepare */
		it.GoToBegin();
		ot.GoToBegin();


		/* the coil array values of the pixel*/
		VectorImagePixelType v;


		/* the output pixel type*/
		VectorImagePixelType o;



		/* BANG!!*/
		while( !it.IsAtEnd() )
		{
			
			ot.Set(it.Get());
			++it;
			++ot;
		}


	}


//std::cout<<output;




 }




template< class VectorImageType>
void iFFTPhasedArrayFilter< VectorImageType>
::ifft3D()
 {

	typedef typename VectorImageType::PixelType    ScalarImagePixelType;

	typename VectorImageType::ConstPointer s = this->GetInput();

	typename VectorImageType::Pointer output = this->GetOutput();
	output->SetBufferedRegion(output->GetRequestedRegion());
	output->Allocate();


	typename VectorImageType::Pointer CIFFT;

	using ifftType= cm::iFFTPhasedArrayFilter2D<VectorImageType>;

	typename ifftType::Pointer ifft=ifftType::New();
	ifft->SetInput(s);
	CIFFT=ifft->GetOutput();
	CIFFT->Update();



	/* instantiate the runners!!*/
	itk::ImageRegionConstIterator<VectorImageType> it(CIFFT,CIFFT->GetLargestPossibleRegion());
	itk::ImageRegionIterator<VectorImageType> ot(output,output->GetLargestPossibleRegion());

	/* prepare */
	it.GoToBegin();
	ot.GoToBegin();


	/* the coil array values of the pixel*/
	VectorImagePixelType v;


	/* the output pixel type*/
	ScalarImagePixelType o;



	/* BANG!!*/
	/* BANG!!*/
		while( !it.IsAtEnd() )
		{
			
			ot.Set(it.Get());
			++it;
			++ot;
		}


 }




template< class VectorImageType>
void iFFTPhasedArrayFilter< VectorImageType>
::GenerateData()
 {
	switch(this->GetKSpaceDimension()){
	case cm::KSpaceAcquisitionDimension::TRIDIMENSIONAL:	this->ifft3D();
	break;
	case cm::KSpaceAcquisitionDimension::BIDIMENSIONAL:
	this->ifft2D();
	break;
	}
 }









}// end namespace


#endif
