#ifndef __cmRootSumOfSquareImageFilter_hxx
#define __cmRootSumOfSquareImageFilter_hxx

#include "itkVectorImage.h"
#include "cmRootSumOfSquareImageFilter.h"
#include "itkImageAlgorithm.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "cmiFFTPhasedArrayFilter.h"

#include "itkExtractImageFilter.h"

#include "itkTileImageFilter.h"

namespace cm
{




template< class TImage,class TOImage>
void RootSumOfSquareImageFilter< TImage,TOImage>
::GenerateData()
 {
			using ifftType= cm::iFFTPhasedArrayFilter<TImage>;
			typename ifftType::Pointer ifft=ifftType::New();
			ifft->SetInput(this->GetInput());
			ifft->SetKSpaceDimension(this->GetKSpaceDimension());
			ifft->Update();

			typename TOImage::Pointer output = this->GetOutput();
			output->SetBufferedRegion(output->GetRequestedRegion());
			output->Allocate();

			//this->GetOutput()->Gtaft(ifft->GetOutput());

					itk::ImageRegionConstIterator<TImage> it(ifft->GetOutput(),ifft->GetOutput()->GetLargestPossibleRegion());
					itk::ImageRegionIterator<TOImage> ot(output,output->GetLargestPossibleRegion());

					/* prepare */   
							it.GoToBegin();
							ot.GoToBegin();


							/* the coil array values of the pixel*/
							VectorImagePixelType v;


							/* the output pixel type*/
							ScalarImagePixelType o;



			/* BANG!!*/
					while( !it.IsAtEnd() )
					{
						v=it.Get();
						/* reset counter and calculate the sos*/
						o=(ScalarImagePixelType) (0.0,0.0);

						/*sum of squares*/
						for (auto t=0;t<v.GetSize();t++)
						{
							o+=std::pow(std::abs(v.GetElement(t)),2);
						}

						/* root*/
						o=std::sqrt(o);
						ot.Set(o);
						++it;
						++ot;
					}





 }









}// end namespace


#endif
