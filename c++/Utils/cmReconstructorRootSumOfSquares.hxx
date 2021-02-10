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



namespace cm
{

template< class TImage,class TOImage>
void ReconstructorRootSumOfSquares< TImage,TOImage>
::GenerateData()
 {
	std::cout<<"RSS recon"<<std::endl;

	typedef typename cm::RootSumOfSquareImageFilter<TImage,TOImage> P2;
	typename P2::Pointer pp2= P2::New();


	pp2->SetInput(this->GetInput());
	pp2->SetKSpaceDimension(this->GetKSpaceDimension());
	pp2->Update();


	this->GetOutput()->Graft(pp2->GetOutput());



 }



}// end namespace


#endif







