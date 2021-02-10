#ifndef __cmReconstructormSENSE_hxx
#define __cmReconstructormSENSE_hxx


#include "cmReconstructormSENSE.h"
#include "itkImageAlgorithm.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"

#include "vnl/algo/vnl_matrix_inverse.h"

#include "ImageUtils.h"

#include "itkShrinkImageFilter.h"
namespace cm
{

template< class VectorImageType,class ScalarImageType>
void ReconstructormSENSE< VectorImageType,ScalarImageType>
::GenerateData()
 {
switch(this->GetKSpaceDimension()){
case cm::KSpaceAcquisitionDimension::BIDIMENSIONAL:
	senseRecon2D();
	break;
}
 }



template< class VectorImageType,class ScalarImageType>
void ReconstructormSENSE<VectorImageType,ScalarImageType>::senseRecon2D()
 {


	VectorImageTypePointer SENSITIVITYMAP=this->GetSensitivityMap();

		typename VectorImageType::Pointer IF=this->GetInputIFFT();//TODO rewrtite it

		typename VectorImageType::IndexType acceleration=this->GetAcceleration();
		const int Kdimension =this->GetInput()->GetNumberOfComponentsPerPixel();
		cm::VectorImagePixelElelmentType c(Kdimension);
		cm::PixelType RTOT=(1/((float)acceleration[0]*acceleration[1]),0.0);
		c.Fill(RTOT);


//		typedef itk::MultiplyImageFilter<VectorImageType,VectorImageType,VectorImageType> Ftype;
//		typename Ftype::Pointer mu = Ftype::New();
//
//			mu->SetInput(IF);
//			mu->SetConstant(c);
//			mu->Update();
//
//			return mu->GetOutput();
//
//
//		this->SetInputIFFT(mu->GetOutput());
//
//		IF=this->GetInputIFFT();

		//FOLDE IT
    	  using ShrinkImageFilterType = itk::ShrinkImageFilter<VectorImageType ,VectorImageType>;
    	  typename ShrinkImageFilterType::Pointer shrinkFilter = ShrinkImageFilterType::New();
    	    shrinkFilter->SetInput(IF);
    	    shrinkFilter->SetShrinkFactor(0, acceleration[0]);
    	    shrinkFilter->SetShrinkFactor(1, acceleration[1]);
    	    shrinkFilter->Update();
		this->SetInputIFFT(shrinkFilter->GetOutput());
		this->CalculateInverseCovariance();



		typename VectorImageType::Pointer ALIASEDKSPACEDATAIFFT=this->GetInputIFFT();
		VectorImageRegionType aliasedRegion=ALIASEDKSPACEDATAIFFT->GetLargestPossibleRegion();

		itk::ImageRegionIteratorWithIndex<VectorImageType> aliasedRunner(ALIASEDKSPACEDATAIFFT,aliasedRegion);



		ScalarImageTypePointer output=this->GetOutput();
		output->SetBufferedRegion(output->GetRequestedRegion());
		output->Allocate();
		aliasedRunner.GoToBegin();


		typedef typename ScalarImageType::IndexType OIndex;
			OIndex	Oi;
			OIndex	POi;
			typedef typename VectorImageType::IndexType Iindex;
			Iindex Ii;

			//
			std::vector<OIndex>VOi;
			//	//* Pixel that  needs the computation */
			int P;

			VectorImagePixelType v;
			VectorImagePixelType v2;

			int AF=acceleration[0];
			int AP=acceleration[1];
			ChannelArrayType sens(Kdimension,AF*AP);
			ChannelArrayType sPinv(AF*AP,Kdimension);
			ChannelArrayType imFold(Kdimension,1);
			ChannelArrayType pp(AF*AP,1);

			std::cout<<"AF"<<AF<<"AP"<<AP<<std::endl;
			while( !aliasedRunner.IsAtEnd() )
				{


				try{
							Ii= aliasedRunner.GetIndex();
						}catch (const std::exception& e)
						{
							std::cerr << e.what();
						}






					/* Gather the aliased pixels into the sensitivity matrix */
					P=0;

					VOi.clear();
					for (auto tf=0; tf<AF;tf++)
					{
						for (auto tp=0;tp<AP;tp++)
						{

							//ouput point ans sens position
							Oi={Ii[0]+(tf*aliasedRegion.GetSize()[0]),Ii[1]+(tp*aliasedRegion.GetSize()[1]),Ii[2]};
							//array of the index point
							VOi.push_back(Oi);
							//
							/* sensitivity map */
							v2=SENSITIVITYMAP->GetPixel(Oi);



							for (int t=0;t<v2.GetSize();t++)
							{

								sens(t,P)=v2.GetElement(t);

							}
							P++;
							//vcl_cout<<sens;

						}
					}



					/* accelerated data*/
					v=aliasedRunner.Get();

					//std::cout<<"v\t"<<v<<"\n\n";

					for (int t=0;t<v.GetSize();t++)
					{

						imFold(t,0)=v.GetElement(t)/RTOT;

					}

					/* %psuedoinverse of coil sensitivity matrix */
					sPinv=vnl_svd_inverse(sens);
					//vcl_cout<<"sPinv\t"<<sPinv;
					/*do the reconstruction of this set of pixels */
					pp=sPinv*imFold;

				//	std::cout<< Ii<<"\t";
					for (auto t=0;t<AF*AP; t++){
				//		std::cout<<VOi.at(t);
						output->SetPixel(VOi.at(t),pp(t,0));
					};

				//	std::cout<<std::endl;
					++aliasedRunner;
				};





 }


//template< class VectorImageType,class ScalarImageType>
//void ReconstructormSENSE< VectorImageType,ScalarImageType>::ThreadedGenerateData(const VectorImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
// {
//
//
//
//
//
//
//
//
// }









//template< class VectorImageType,class ScalarImageType>
//void ReconstructormSENSE< VectorImageType,ScalarImageType>
//::senseRecon2D(const VectorImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
// {
//
//
//
//
//
//
//
//
//	VectorImageTypePointer SENSITIVITYMAP=this->GetSensitivityMap();
//	VectorImageTypePointer KSPACEDATAIFFTALIASED=this->GetInputIFFT();
//	const int Kdimension =KSPACEDATAIFFTALIASED->GetNumberOfComponentsPerPixel();
//	ScalarImageTypePointer OUT=this->GetOutput();
//    ChannelArrayType iCM=this->GetInverseNoiseCovariance();
//	const int Sdimension =SENSITIVITYMAP->GetNumberOfComponentsPerPixel();
//
//
//	VectorImageRegionType aliasedRegion=ALIASEDKSPACEDATAIFFT->GetLargestPossibleRegion();



	//cm::VectorImageTypePointer SENSITIVITYMAP=this->GetSensitivityMap();
//
//std::cout<<"got the sensitivity map for sense reconstruction\\n\n";
//cm::VectorImageTypePointer ALIASEDKSPACEDATAIFFT=this->GetInputIFFT();
//std::cout<<"got aliased ifft for sense reconstruction\\n\n";
//
//cm::ChannelArrayType iCM=this->GetInverseNoiseCovariance();
//
//std::cout<<"got inverse covariance for sense reconstruction\n\n";
//
//ScalarImageTypePointer output=this->GetOutput();
//output->SetBufferedRegion(output->GetRequestedRegion());
//output->Allocate();
//
//std::cout<<"sensitivitymao\n\n";
//
//const int Kdimension =ALIASEDKSPACEDATAIFFT->GetNumberOfComponentsPerPixel();
//const int Sdimension =SENSITIVITYMAP->GetNumberOfComponentsPerPixel();
//
//
//
//VectorImageRegionType aliasedRegion=ALIASEDKSPACEDATAIFFT->GetLargestPossibleRegion();
//
//
////aliased data
//itk::ImageRegionIteratorWithIndex<VectorImageType> aliasedRunner(ALIASEDKSPACEDATAIFFT,aliasedRegion);
//
//
//
//int AF=this->GetFrequencyAcceleration();
//int AP=this->GetPhaseAcceleration0();
//
//
//	aliasedRunner.GoToBegin();
//	//	snr.GoToBegin();
//
//	typedef typename ScalarImageType::IndexType OIndex;
//	OIndex	Oi;
//	OIndex	POi;
//	typedef typename VectorImageType::IndexType Iindex;
//	Iindex Ii;
//
//	//
//	std::vector<OIndex>VOi;
//	//	//* Pixel that  needs the computation */
//	int P;
//
//	VectorImagePixelType v;
//	VectorImagePixelType v2;
//
//	ChannelArrayType sens(Kdimension,AF*AP);
//	ChannelArrayType sPinv(AF*AP,Kdimension);
//	ChannelArrayType imFold(Kdimension,1);
//	ChannelArrayType pp(AF*AP,1);
//
//	std::cout<<"AF"<<AF<<"AP"<<AP<<std::endl;
//	while( !aliasedRunner.IsAtEnd() )
//		{
//
//
//		try{
//					Ii= aliasedRunner.GetIndex();
//				}catch (const std::exception& e)
//				{
//					std::cerr << e.what();
//				}
//
//
//
//
//
//
//			/* Gather the aliased pixels into the sensitivity matrix */
//			P=0;
//
//			VOi.clear();
//			for (auto tf=0; tf<AF;tf++)
//			{
//				for (auto tp=0;tp<AP;tp++)
//				{
//
//					//ouput point ans sens position
//					Oi={Ii[0]+(tf*aliasedRegion.GetSize()[0]),Ii[1]+(tp*aliasedRegion.GetSize()[1]),Ii[2]};
//					//array of the index point
//					VOi.push_back(Oi);
//					//
//					/* sensitivity map */
//					v2=SENSITIVITYMAP->GetPixel(Oi);
//
//
//
//					for (int t=0;t<v2.GetSize();t++)
//					{
//
//						sens(t,P)=v2.GetElement(t);
//
//					}
//					P++;
//					//vcl_cout<<sens;
//
//				}
//			}
//
//
//
//			/* accelerated data*/
//			v=aliasedRunner.Get();
//
//			//std::cout<<"v\t"<<v<<"\n\n";
//
//			for (int t=0;t<v.GetSize();t++)
//			{
//
//				imFold(t,0)=v.GetElement(t);
//
//			}
//
//			/* %psuedoinverse of coil sensitivity matrix */
//			sPinv=vnl_svd_inverse(sens);
//			//vcl_cout<<"sPinv\t"<<sPinv;
//			/*do the reconstruction of this set of pixels */
//			pp=sPinv*imFold;
//
//		//	std::cout<< Ii<<"\t";
//			for (auto t=0;t<AF*AP; t++){
//		//		std::cout<<VOi.at(t);
//				output->SetPixel(VOi.at(t),pp(t,0));
//			};
//
//		//	std::cout<<std::endl;
//			++aliasedRunner;
//		};





// }

//template< class VectorImageType,class ScalarImageType>
//void ReconstructormSENSE< VectorImageType,ScalarImageType>
//::GenerateOutputInformation()
// {
	//we want the output image to be of the same size of the sensibility map

	// Call the superclass' implementation of this method
//	Superclass::GenerateOutputInformation();

//
//
//	// Get pointers to the input and output
//	VectorImageTypePointer sensPtr =   this->GetSensitivityMap();
//	ScalarImageTypePointer outputPtr = this->GetOutput();
//
//
//
//	// We need to compute the output spacing, the output image size, and the
//	// output image start index
//	const typename VectorImageType::SpacingType &
//	inputSpacing = sensPtr->GetSpacing();
//	const typename VectorImageType::SizeType &   inputSize =
//			sensPtr->GetLargestPossibleRegion().GetSize();
//	const typename VectorImageType::IndexType &  inputStartIndex =
//			sensPtr->GetLargestPossibleRegion().GetIndex();
//	const typename VectorImageType::PointType &
//	inputOrigin = sensPtr->GetOrigin();
//
//	typename ScalarImageType::SpacingType outputSpacing;
//	typename ScalarImageType::SizeType outputSize;
//	typename ScalarImageType::IndexType outputStartIndex;
//	typename ScalarImageType::PointType outputOrigin;
//	//typename ScalarImageType::DirectionType outputDirection;
//
//	typename VectorImageType::SpacingType inputOriginShift;
//
//
//	outputSpacing= inputSpacing;
//	outputSize = inputSize;
//	outputStartIndex = inputStartIndex;
//
//
//
//
//
//	outputPtr->SetSpacing(outputSpacing);
//	outputPtr->SetOrigin(inputOrigin);
//	outputPtr->SetDirection(sensPtr->GetDirection());
//
//	typename ScalarImageType::RegionType outputLargestPossibleRegion;
//	outputLargestPossibleRegion.SetSize(outputSize);
//	outputLargestPossibleRegion.SetIndex(outputStartIndex);
//
//	//std::cout<<outputLargestPossibleRegion;
//
//	outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
//	//std::cout<<"qui tutto apposto"<<std::endl;
// }









/*
  GenerateInputRequesteRegion
 */
//template< class VectorImageType,class ScalarImageType>
//void ReconstructormSENSE< VectorImageType,ScalarImageType>
//::GenerateInputRequestedRegion()
// {
	// Call the superclass' implementation of this method
//	Superclass::GenerateInputRequestedRegion();

//	// Get pointers to the input and output
//	VectorImageTypePointer inputPtr =
//			const_cast< VectorImageType * >( this->GetInput() );
//	ScalarImageTypePointer outputPtr = this->GetOutput();
//
//
//
//	// We need to compute the input requested region (size and start index)
//	unsigned int i;
//	const typename ScalarImageType::SizeType & outputRequestedRegionSize =
//			outputPtr->GetRequestedRegion().GetSize();
//	const typename ScalarImageType::IndexType & outputRequestedRegionStartIndex =
//			outputPtr->GetRequestedRegion().GetIndex();
//
//	typename VectorImageType::SizeType inputRequestedRegionSize;
//	typename VectorImageType::IndexType inputRequestedRegionStartIndex;
//	typedef double SizeValueType;
//
//	/*
//	  inputRequestedSize = (outputRequestedSize / ExpandFactor) + 1)
//	  The extra 1 above is to take care of edge effects when streaming.
//	 */
//
//	inputRequestedRegionSize[0] =
//			(SizeValueType)std::ceil( (double)outputRequestedRegionSize[0]
//																		/ (double)m_FrequencyAcceleration ) + 1;
//
//	inputRequestedRegionSize[1] =
//			(SizeValueType)std::ceil( (double)outputRequestedRegionSize[1]
//																		/ (double)m_PhaseAcceleration0 ) + 1;
//
//	inputRequestedRegionSize[2] =
//			(SizeValueType)std::ceil( (double)outputRequestedRegionSize[2]
//																		/ (double)m_PhaseAcceleration1 ) + 1;
//
//
//	inputRequestedRegionStartIndex[0] =
//			(SizeValueType)std::floor( (double)outputRequestedRegionStartIndex[0]
//																			   / (double)m_FrequencyAcceleration );
//
//
//	inputRequestedRegionStartIndex[1] =
//			(SizeValueType)std::floor( (double)outputRequestedRegionStartIndex[1]
//																			   / (double)m_PhaseAcceleration0 );
//
//	inputRequestedRegionStartIndex[2] =
//			(SizeValueType)std::floor( (double)outputRequestedRegionStartIndex[2]
//																			   / (double)m_PhaseAcceleration1 );
//
//
//
//	typename VectorImageType::RegionType inputRequestedRegion;
//	inputRequestedRegion.SetSize(inputRequestedRegionSize);
//	inputRequestedRegion.SetIndex(inputRequestedRegionStartIndex);
//
//	// Make sure the requested region is within largest possible.
//	inputRequestedRegion.Crop( inputPtr->GetLargestPossibleRegion() );
//
//	// Set the input requested region.
//	inputPtr->SetRequestedRegion(inputRequestedRegion);
// }








}// end namespace


#endif
