#ifndef __cmSNRUnitsReconstructorRootSumOfSquares_hxx
#define __cmSNRUnitsReconstructorRootSumOfSquares_hxx


#include "cmSNRUnitsReconstructorRootSumOfSquares.h"
#include "itkImageAlgorithm.h"
#include "itkImageRegionIterator.h"
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/algo/vnl_matrix_inverse.h"

#include "cmRootSumOfSquareImageFilter.h"



namespace cm
{

template< class TImage,class TOImage>
void SNRUnitsReconstructorRootSumOfSquares< TImage,TOImage>
::ThreadedGenerateData(const VectorImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
 {



	VectorImageTypePointer KSPACEDATAIFFT=this->GetInputIFFT();
	ScalarImageTypePointer OUT=this->GetOutput();
	ChannelArrayType iCM=this->GetInverseNoiseCovariance();
		const int Kdimension =KSPACEDATAIFFT->GetNumberOfComponentsPerPixel();


				itk::ImageRegionIterator<VectorImageType> iif(KSPACEDATAIFFT,outputRegionForThread);
				itk::ImageRegionIterator<ScalarImageType> io(OUT,outputRegionForThread);

//				/* prepare */
				iif.GoToBegin();
				io.GoToBegin();





				 VectorImagePixelType I;
				 ChannelArrayType IV(Kdimension,1);

				 ChannelArrayType OV(1,1);
				 ChannelArrayType REG(1,1);
				 REG(0,0)= 2.0;


				/* BANG!!*/
								while( !iif.IsAtEnd() )
								{
									I=iif.Get();

									/* reset counter and calculate the sos*/

									for(auto y=0;y<Kdimension;y++)
									{
										IV(y,0)=I.GetElement(y);

									}

									OV=IV.conjugate_transpose()* IV;
									io.Set((ScalarImagePixelType)std::abs(std::sqrt( OV(0,0)* REG(0,0))));
									++iif;
									++io;
								}



 }




}// end namespace


#endif







