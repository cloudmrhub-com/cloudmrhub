#include "itkImage.h"
#include "itkComposeImageFilter.h"
#include "itkVectorImage.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkVectorImage.h"
#include <stdlib.h>     /* srand, rand */


namespace
{
using PT = float;
//	using PixelType = std::complex< PT >;
using PixelType = std::complex<PT>;
using VectorImageType = itk::VectorImage<PixelType, 3>;
using ScalarImageType = itk::Image<PixelType, 3>;
typedef itk::ImageFileWriter<VectorImageType> WFilter;
using ConstIteratorType = itk::ImageRegionConstIterator< ScalarImageType >;
using VConstIteratorType = itk::ImageRegionConstIterator< VectorImageType >;
using VectorIteratorType = itk::ImageRegionIterator< VectorImageType >;
using ScalarIteratorType = itk::ImageRegionIterator< ScalarImageType >;

} // namespace

static void
CreateImage(ScalarImageType::Pointer image,ScalarImageType::InternalPixelType type);

int main(int argc, char *argv[])
{
  ScalarImageType::Pointer image0 = ScalarImageType::New();
  CreateImage(image0,(ScalarImageType::InternalPixelType) (5,5));

  ScalarImageType::Pointer image1 = ScalarImageType::New();
  CreateImage(image1,(ScalarImageType::InternalPixelType) (2,2));

//  ScalarImageType::Pointer image2 = ScalarImageType::New();
//  CreateImage(image2,(ScalarImageType::InternalPixelType) (3,2));

  using ImageToVectorImageFilterType = itk::ComposeImageFilter<ScalarImageType>;
  ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();
  imageToVectorImageFilter->SetInput(0, image0);
  imageToVectorImageFilter->SetInput(1, image1);
//  imageToVectorImageFilter->SetInput(2, image2);
  imageToVectorImageFilter->Update();

  VectorImageType::Pointer vectorImage = imageToVectorImageFilter->GetOutput();


  WFilter::Pointer filter = WFilter::New();
  	filter->SetInput(vectorImage);
  	filter->SetFileName(argv[1]);
  	filter->Update();


//  	VectorImageType::SizeType size =vectorImage->GetLargestPossibleRegion().GetSize();
//  	VectorImageType::IndexType PT;
//  	  for (long long y=0; y<size[1];y++){ //for each phase encoding
//  	  		for (long long x=0; x<size[0];x++){
//
//  	  			PT[0]=x;
//  	  			PT[1]=y;
//  	  			PT[2]=0;
//  	  			std::cout<<vectorImage->GetPixel(PT)<<", ";
//
//  	  		}
//  	  	std::cout<<std::endl;
//  	  	};
//  	std::cout<<"\n\n"<<std::endl;
//
//
//
//  	 VConstIteratorType inputIt(vectorImage, vectorImage->GetLargestPossibleRegion());
//
//  	  		      inputIt.GoToBegin();
//
//
//  	  		      while( !inputIt.IsAtEnd() )
//  	  		        {
//  	  		   std::cout<<inputIt.Get()  <<std::endl;
//  	  		        ++inputIt;
//  	  		        }




  return EXIT_SUCCESS;
}

void
CreateImage(ScalarImageType::Pointer image,ScalarImageType::InternalPixelType t)
{
  ScalarImageType::IndexType start;
  start.Fill(0);

  ScalarImageType::SizeType size;
  size[0]=10;
  size[1]=5;
  size[2]=1;

  ScalarImageType::RegionType region(start, size);

  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer((ScalarImageType::PixelType)t);
  
//  ScalarImageType::IndexType PT;
//  for (long long y=0; y<size[1];y++){ //for each phase encoding
//  		for (long long x=0; x<size[0];x++){
//
//  			PT[0]=x;
//  			PT[1]=y;
//  			PT[2]=0;
//  			std::cout<<image->GetPixel(PT)<<", ";
//
//  		}
//  		std::cout<<std::endl;
//  		  	  	};
//  		  	std::cout<<"\n\n"<<std::endl;



  ScalarIteratorType inputIt(image, region);

//
//  ScalarImageType::PixelType P;
//  ScalarImageType::PixelType R;
//  		      inputIt.GoToBegin();
//
//  		    srand (time(NULL));
//  		      while( !inputIt.IsAtEnd() )
//  		        {
//  		    	  P=inputIt.Get();
//  		   std::cout<<P<<"-> ";
//
//
//  		   R=((float) rand(),(float) rand());
//  		   R= P + R;
//
//  		 inputIt.Set(R);
//  		 std::cout<<inputIt.Get()  <<std::endl;
//  		        ++inputIt;
//  		        }




  }
