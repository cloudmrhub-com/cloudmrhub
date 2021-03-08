#ifndef __cmReconstructorRootSumOfSquares_hxx
#define __cmReconstructorRootSumOfSquares_hxx


#include "cmReconstructorRootSumOfSquares.h"
#include "itkImageAlgorithm.h"
#include "itkImageRegionIterator.h"
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/algo/vnl_matrix_inverse.h"

#include "cmRootSumOfSquareImageFilter.h"
#include "itkPasteImageFilter.h"


namespace cm
{

template< class TImage,class TOImage>
void ReconstructorRootSumOfSquares< TImage,TOImage>
::ThreadedGenerateData(const VectorImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
 {


	VectorImageTypePointer KSPACEDATAIFFT=this->GetInputIFFT();

	const int Kdimension =KSPACEDATAIFFT->GetNumberOfComponentsPerPixel();

	ScalarImageTypePointer OUT=this->GetOutput();



				itk::ImageRegionIterator<VectorImageType> it(KSPACEDATAIFFT,outputRegionForThread);
				itk::ImageRegionIterator<ScalarImageType> io(OUT,outputRegionForThread);

				it.GoToBegin();
				io.GoToBegin();





				 VectorImagePixelType I;
				 ChannelArrayType IV(Kdimension,1);

				 ChannelArrayType OV(1,1);


		/* BANG!!*/
				while( !it.IsAtEnd() )
				{
					I=it.Get();
					/* reset counter and calculate the sos*/
					OV=(ScalarImagePixelType) (0.0,0.0);

					/*sum of squares*/
					for (auto t=0;t<I.GetSize();t++)
					{
						OV(0,0)+=std::pow(std::abs(I.GetElement(t)),2);
					}

					/* root*/
					OV=std::sqrt(OV(0,0));
					io.Set(OV(0,0));
					++it;
					++io;
				}




 }




}// end namespace


#endif





