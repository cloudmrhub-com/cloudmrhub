#ifndef __cmPseudoMultipleReplicasFilter_hxx
#define __cmPseudoMultipleReplicasFilter_hxx


//#include "itkChildTreeIterator.h"
#include "cmPseudoMultipleReplicasFilter.h"
#include "itkProgressReporter.h"
#include "itkImageAlgorithm.h"

//#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkProgressReporter.h"
#include "vnl/vnl_sample.h"


#include "itkImageDuplicator.h"

#include "itkImageRegionIterator.h"
namespace cm
{

template< class ScalarImageType>
typename cm::PseudoMultipleReplicasFilter<ScalarImageType>::VectorImageTypePointer PseudoMultipleReplicasFilter< ScalarImageType>::getPseudoReplicaKspace(VectorImageTypePointer S)
{

	//duplicate the kspace data;

	using DuplicatorType = itk::ImageDuplicator< VectorImageType>;
	typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(S);
	duplicator->Update();



	int NC=S->GetNumberOfComponentsPerPixel();

	//std::cout<<NC;
	VectorImageTypePointer SpN=duplicator->GetOutput();




	typename itk::VariableLengthVector<ScalarImagePixelType>V(NC);
	ScalarImagePixelType n;
	vnl_vector<ScalarImagePixelType>D(NC);


	itk::ImageRegionIterator< VectorImageType > it(SpN, SpN->GetLargestPossibleRegion());
	for (it.Begin(); !it.IsAtEnd(); ++it)
	{

		V=it.Get();

		for (auto v=0;v<V.GetSize();v++)
		{

			D(v)=ScalarImagePixelType(std::sqrt(0.5)*vnl_sample_normal(0.0,1.0),vnl_sample_normal(0.0,1.0));


		}



		//vcl_cout<<"\n"<<vcl_endl;
		D=this->GetCorrelationNoiseFactor()*D;


		for (auto v=0;v<V.GetSize();v++)
				{
			V.SetElement(v,V.GetElement(v)+D(v));
				}
		it.Set(V);

	}






	return SpN;


}

template< class ScalarImageType>
void PseudoMultipleReplicasFilter< ScalarImageType>::GetN()

//void PseudoMultipleReplicasFilter< ScalarImageType>::GenerateData()
{

	std::cout<< vnl_sample_normal(0.0,1.0)<<std::endl;


}



}// end namespace


#endif
