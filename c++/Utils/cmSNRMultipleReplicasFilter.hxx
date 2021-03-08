#ifndef __cmSNRMultipleReplicasFilter_hxx
#define __cmSNRMultipleReplicasFilter_hxx


//#include "itkChildTreeIterator.h"
#include "cmSNRMultipleReplicasFilter.h"
#include "itkProgressReporter.h"
#include "itkImageAlgorithm.h"

//#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "mymath.h"
#include "itkProgressReporter.h"

//#include "vnl_sample.h"
namespace cm
{

template< class ScalarImageType>
void SNRMultipleReplicasFilter< ScalarImageType>::ThreadedGenerateData(const ScalarImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)

//void SNRMultipleReplicasFilter< ScalarImageType>::GenerateData()
{

	//  // Support progress methods/callbacks
	//  itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels() );
	//
	//  // Iterate through the thread region

	ScalarImageIteratorType outputIt(this->GetOutput(), outputRegionForThread);

	std::cout<<this->GetOutput()<<std::endl;


	// Grab input image

	// Support progress methods/callbacks
	itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels() );


	for ( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt)
	{

		//std::cout<<this->getReplicasonPoint(outputIt.GetIndex())[0]<<std::endl;

		//std::vector<float>X=this->getReplicasonPointReal(outputIt.GetIndex());

		std::vector<ScalarImageType::InternalPixelType> X=this->getReplicasonPoint(outputIt.GetIndex());



//		std::cout<<CALC::mean<ScalarImagePixelType>(this->getReplicasonPoint(outputIt.GetIndex()))<<std::endl;

		outputIt.Set(CALC::mean(X)/CALC::stdvN(X));
		progress.CompletedPixel();
	} // Next pixel
}





template< class ScalarImageType>
//void SNRMultipleReplicasFilter< ScalarImageType>::ThreadedGenerateData(const ScalarImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
std::vector<float> SNRMultipleReplicasFilter< ScalarImageType>::getReplicasonPointReal(ScalarImageIndexType p)
{
	std::vector<float> V(m_stack.size());

	for(auto a=0;a<this->m_stack.size();a++)
	{
		V.at(a)=(float)std::real(m_stack.at(a)->GetPixel(p));

	}

	return V;
}


template< class ScalarImageType>
//void SNRMultipleReplicasFilter< ScalarImageType>::ThreadedGenerateData(const ScalarImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
std::vector<ScalarImageType::InternalPixelType> SNRMultipleReplicasFilter< ScalarImageType>::getReplicasonPoint(ScalarImageIndexType p)
{
	VectorofScalarInternalImagePixelType V(m_stack.size());

	for(auto a=0;a<this->m_stack.size();a++)
	{
		V.at(a)=m_stack.at(a)->GetPixel(p);

	}

	return V;
}

template< class ScalarImageType>
void SNRMultipleReplicasFilter< ScalarImageType>::pushReconstructedImage(ScalarImageTypePointer im)
{

	this->m_stack.push_back(im);

	std::cout<<"stacksize is"<<this->m_stack.size()<<std::endl;
	//std::cout<<"PROVA:"<<vnl_sample_norm(0,20)<<std::endl;


}




}// end namespace


#endif
