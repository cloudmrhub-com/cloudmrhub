#ifndef __cmVectorImageRootSumOfSquaresReconstructor_hxx
#define __cmVectorImageRootSumOfSquaresReconstructor_hxx


#include "cmVectorImageRootSumOfSquaresReconstructor.h"
#include "itkImageAlgorithm.h"
#include "itkImageRegionIterator.h"
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/algo/vnl_matrix_inverse.h"

#include "cmSumOfSquareImageFilter.h"



namespace cm
{

template< class TImage,class TOImage>
void VectorImageRootSumOfSquaresReconstructor< TImage,TOImage>
::GenerateData()
 {
	std::cout<<"RSS recon"<<std::endl;

	typedef typename cm::SumOfSquareImageFilter<TImage,TOImage> P2;
	typename P2::Pointer pp2= P2::New();


	pp2->SetInput(this->GetInput());
	pp2->SetKSpaceDimension(this->GetKSpaceDimension());
	pp2->Update();


	this->GetOutput()->Graft(pp2->GetOutput());



 }



}// end namespace


#endif







