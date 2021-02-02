/*
 * ISMRMRDToITKFilter.h
 *
 *  Created on: Jan 31, 2019
 *      Author: montie01
 */
//do not parelelize the reader!!
#ifndef __cmISMRMRDToITKImageFilter_txx
#define __cmISMRMRDToITKImageFilter_txx

#include "cmISMRMRDToITKImageFilter.h"

#include "cmSumOfSquareImageFilter.h"
#pragma once
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/dataset.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/version.h"

#include "itkPoint.h"
#include "itkPointSet.h"
#include "itkBoundingBox.h"
#include "utils.h"

#include <math.h>
#include <iostream>
namespace cm
{


template< class TImage>
std::vector<TImage::Pointer> ISMRMRDToITKImageFilter< TImage>
::getReplicas(std::string Type){

	std::cout<<"Repetition"<<std::endl;

	std::vector<TImage::Pointer> OUT(this->getNumberOfRepetition());


	//the acqisition

	//must be instantiated with a filename
	ISMRMRD::Dataset ismrmdataset(this->GetIsmrmrdFileName().c_str(),"dataset", false);

	//create the image vector
	typename VectorImageType::SizeType SIZE;
	VectorImagePointerType	IM;


	for (auto t=0; t<this->getNumberOfRepetition(); t++)
	{
	if(this->GetremoveFrequencyOversampling()){
		IM=this->getNewReconVectorImage();
		SIZE= this->getReconSize();
	}else{
		IM=this->getNewEncodedVectorImage();
		SIZE= this->getEncodedSize();
	}
OUT.at(t)=IM;
	}




	//ORIGIN is a point cloud i've to find the minimum later
	using CoordType = float;
	const unsigned int Dimension = VectorImageType::ImageDimension;


	using PointSetType = itk::PointSet< CoordType, Dimension >;

	using PointIdentifier = typename PointSetType::PointIdentifier;
	using PointType = typename PointSetType::PointType;
	using PointsContainerPointer = typename PointSetType::PointsContainerPointer;
	typename PointSetType::Pointer   pointSet = PointSetType::New();
	PointsContainerPointer  points = pointSet->GetPoints();

	unsigned long int pointId = 0;


	typename VectorImageType::DirectionType DIRECTIONS;

	unsigned short int numberofcoil=this->getNumberOfCoils();

	unsigned short int coilnumber;


	long int distance;






	//just ask once:)
	int osfactor=1;
	if(this->GetremoveFrequencyOversampling()){
		osfactor=2;
	};



	VectorImageInnerPixelType NR=1;
	if(this->GetaverageImageThroughAverages()){
		NR=(VectorImageInnerPixelType)((float)this->getNumberOfAverage());
	};


	//omp_set_num_threads(omp_get_num_procs());
	//
	// #pragma omp parallel for private(vectorValue)
	ISMRMRD::ISMRMRD_EncodingCounters counter;
	for (auto p=0; p<ismrmdataset.getNumberOfAcquisitions();p++)
	{
		typename VectorImageType::PixelType   vectorValue(numberofcoil);
		ISMRMRD::Acquisition acq;
		ismrmdataset.readAcquisition(p, acq);

		if (this->isTheAcquisitionRequestedImageReplicas(acq) && !(this->isTheAcquisitionNoise(acq)))
		{
//			this->AcquisitionInfo(acq);
			ISMRMRD::ISMRMRD_EncodingCounters counter;
			counter=acq.idx();



			IM=OUT.at(counter.repetition);


			//since it is bidimensional acquisition the 1nd index represent the frequency encoding and will change in the for loop
			typename VectorImageType::IndexType index;

			index={0,counter.kspace_encode_step_1,counter.slice};
			//for all the frequency encoded

			for(auto idx=0;idx<SIZE[0]*osfactor;idx=idx+osfactor){
				index[0]=idx/osfactor;
				vectorValue=IM->GetPixel(   index  );

				for(auto ch=0;ch<numberofcoil;ch++){
					vectorValue[ch] +=  (acq.data(idx,ch)/NR);
				}
				IM->SetPixel(   index,   vectorValue  );

			}

			DIRECTIONS=this->getAcquisitionDirections(acq);

			points->InsertElement( pointId, this->getAcquisitionOrigin(acq) );
			pointId++;

			OUT.at(counter.repetition)=IM;
		}
	};


	IM=OUT.at(0);
	//set the directions with the last acq
	IM->SetDirection(DIRECTIONS);

	using BoundingBoxType = itk::BoundingBox< PointIdentifier, Dimension, CoordType >;
	typename  BoundingBoxType::Pointer boundingBox = BoundingBoxType::New();
	boundingBox->SetPoints(points);
	boundingBox->ComputeBoundingBox();


	typename VectorImageType::PointType min=boundingBox->GetMinimum();
	IM->SetOrigin(min);

	if(this->GetremoveFrequencyOversampling()){
		typename VectorImageType::PointType oVSP;

		std::cout<<"remove OS origin was"<<min<<std::endl;


		typename VectorImageType::IndexType ind={floor(SIZE[0]/2),0,0};
		IM->TransformIndexToPhysicalPoint(ind,oVSP);


		for (auto n=0;n<3;n++){
			min[n]=min[n]+(std::sqrt(std::pow((oVSP[n]-min[n]),2)));

		}
	}



	for(int c=0;c<OUT.size();c++){
		OUT.at(c)->SetDirection(DIRECTIONS);
		OUT.at(c)->SetOrigin(min);
	}

	return OUT;
}





template< class TImage>
void ISMRMRDToITKImageFilter< TImage>
::AcquisitionInfo(ISMRMRD::Acquisition acq)
 {
	//std::cout<<"N data element:"<<acq.getNumberOfDataElements()<<std::endl;
	//std::cout<<"available chan:" <<acq.available_channels()<<std::endl;




	typename VectorImageType::PointType position,readir,phadir,sldir,bed;
	//		position=*acq.position();
	//		std::cout<<"position:" <<position<<std::endl;
	//
	//		readir=*acq.read_dir();
	//		std::cout<<"read_dir:" <<readir<<std::endl;
	//
	//		phadir=*acq.phase_dir();
	//		std::cout<<"phase_dir:" <<phadir<<std::endl;
	//
	//		sldir=*acq.slice_dir();
	//		std::cout<<"sl_dir:" <<sldir<<std::endl;
	//
	//		bed=*acq.patient_table_position();
	//		std::cout<<"bed:" <<bed<<std::endl;
	//
	//
	//		std::cout<<" Data size: "<<acq.getDataSize()<<std::endl;
	//
	//		std::cout<<" Data size: "<<acq.getNumberOfDataElements()<<std::endl;

	ISMRMRD::ISMRMRD_EncodingCounters counter;
	counter=acq.idx();


	ISMRMRD::AcquisitionHeader H;

	H=acq.getHead();

	std::cout<<"step1"<<counter.kspace_encode_step_1<<","
			<<"step2"<<counter.kspace_encode_step_2<<","
			<<"slice"<<counter.slice<<","
			<<"cnt"<<counter.contrast<<","
			<<"ph"<<counter.phase<<","
			<<"rep"<<counter.repetition<<","
			<<"set"<<counter.set<<","
			<<"segm"<<counter.segment<<"\t"
			<<"flag"<<H.flags<<std::endl;




	std::cout<<H.isFlagSet(ISMRMRD::ISMRMRD_AcquisitionFlags::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT)<<std::endl;




	int EE=0;
	for(auto a=acq.data_begin();a<acq.data_end();a++){
		EE++;
		//std::cout<<*a<<std::endl;
	}
	std::cout<<EE<<std::endl;


 }
template< class TImage>
std::vector<float>  ISMRMRDToITKImageFilter< TImage>
::getSliceOrder()
 {

	std::cout<<"sorting slices"<<std::endl;
	//must be instantiated with a filename
	 ISMRMRD::Dataset ismrmdataset(this->GetIsmrmrdFileName().c_str(),"dataset", false);

	 std::vector<float> slices;

	 typename VectorImageType::PointType position;

	 //cicle on each acquisition and get vector z
	 for (auto p=0; p<ismrmdataset.getNumberOfAcquisitions();p++)
	  {
		 std::cout<<"iteration: "<<p<<" over "<<ismrmdataset.getNumberOfAcquisitions();
	 	 ISMRMRD::Acquisition acq;
	 	 ISMRMRD::AcquisitionHeader O;
	 	 //read the acquisition
	 	 ismrmdataset.readAcquisition(p, acq);
	 	 this->AcquisitionInfo(acq);
	 	 if (this->isTheAcquisitionRequestedImage(acq) && !(this->isTheAcquisitionNoise(acq)))
	 	 {

	 		 O=acq.getHead();

	 		 position=O.position;
	 		 std::cout<<"---position: "<<position[2]<<std::endl;
	 		 slices.push_back(position[2]);
	 	 }
	  }


	 std::sort( slices.begin(), slices.end() );

	 slices.erase( std::unique( slices.begin(), slices.end() ), slices.end() );
	 std::cout<<"slices ordered"<<std::endl;

	 return slices;

 };

template< class TImage>
typename TImage::Pointer  ISMRMRDToITKImageFilter< TImage>
::bidimensionalRead()
 {	std::cout<<"bidimensional"<<std::endl;
 //the acqisition

 //must be instantiated with a filename
 ISMRMRD::Dataset ismrmdataset(this->GetIsmrmrdFileName().c_str(),"dataset", false);

 //create the image vector
 typename VectorImageType::SizeType SIZE;
 VectorImagePointerType	IM;

 //just ask once:)
 int osfactor=1;


 if(this->GetremoveFrequencyOversampling() && this->canRemoveFrequencyOS()){
	 IM=this->getNewReconVectorImage();
	 SIZE= this->getReconSize();
	 osfactor=2;
 }else{
	 IM=this->getNewEncodedVectorImage();
	 SIZE= this->getEncodedSize();
	 osfactor=1;
 }




 //ORIGIN is a point cloud i've to find the minimum later
 using CoordType = float;
 const unsigned int Dimension = VectorImageType::ImageDimension;


 using PointSetType = itk::PointSet< CoordType, Dimension >;

 using PointIdentifier = typename PointSetType::PointIdentifier;
 using PointType = typename PointSetType::PointType;
 using PointsContainerPointer = typename PointSetType::PointsContainerPointer;
 typename PointSetType::Pointer   pointSet = PointSetType::New();
 PointsContainerPointer  points = pointSet->GetPoints();

 unsigned long int pointId = 0;


 typename VectorImageType::DirectionType DIRECTIONS;

 unsigned short int numberofcoil=this->getNumberOfCoils();

 unsigned short int coilnumber;


 long int distance;







std::vector<float> SLICEPOSITION=this->getSliceOrder();





for( unsigned i = 0; i < SLICEPOSITION.size(); ++i ) {
std::cout<<"Position: "<<i<<" that's the position"<<SLICEPOSITION[i]<<std::endl;
};

 VectorImageInnerPixelType NR=1;
 if(this->GetaverageImageThroughAverages()){
	 NR=(VectorImageInnerPixelType)((float)this->getNumberOfAverage());
 };

 std::cout<<"number of averages used:"<<	NR <<"over:"<<this->getNumberOfAverage()<<std::endl;

 std::cout<<"N of acq:"<<ismrmdataset.getNumberOfAcquisitions()<<std::endl;
 typename VectorImageType::PointType position;



 //omp_set_num_threads(omp_get_num_procs());
 //
 // #pragma omp parallel for private(vectorValue)
 for (auto p=0; p<ismrmdataset.getNumberOfAcquisitions();p++)
 {


	 typename VectorImageType::PixelType   vectorValue(numberofcoil);
	 ISMRMRD::Acquisition acq;

	 ismrmdataset.readAcquisition(p, acq);


	 this->AcquisitionInfo(acq);
	 if (this->isTheAcquisitionRequestedImage(acq) && !(this->isTheAcquisitionNoise(acq)))
	 {
		 ISMRMRD::ISMRMRD_EncodingCounters counter;
		 counter=acq.idx();
		 //since it is bidimensional acquisition the 1nd index represent the frequency encoding and will change in the for loop
		 typename VectorImageType::IndexType index;
		 //to fix the interlieved acquisition

		 ISMRMRD::AcquisitionHeader O;
		 O=acq.getHead();

		position=O.position;

		 int slicepos = std::distance(SLICEPOSITION.begin(), std::find(SLICEPOSITION.begin(), SLICEPOSITION.end(),position[2] ));
		 //counter.slice
		 index={0,counter.kspace_encode_step_1,slicepos};

		 //for all the frequency encoded

		 for(auto idx=0;idx<SIZE[0]*osfactor;idx=idx+osfactor){
			 index[0]=idx/osfactor;
			 vectorValue=IM->GetPixel(   index  );

			 for(auto ch=0;ch<numberofcoil;ch++){
				 vectorValue[ch] +=  (acq.data(idx,ch)/NR);
			 }
			 IM->SetPixel(   index,   vectorValue  );

		 }

		 DIRECTIONS=this->getAcquisitionDirections(acq);

		 points->InsertElement( pointId, this->getAcquisitionOrigin(acq) );
		 pointId++;
	 }
 };


 //set the directions with the last acq
 IM->SetDirection(DIRECTIONS);

 using BoundingBoxType = itk::BoundingBox< PointIdentifier, Dimension, CoordType >;
 typename  BoundingBoxType::Pointer boundingBox = BoundingBoxType::New();
 boundingBox->SetPoints(points);
 boundingBox->ComputeBoundingBox();


 typename VectorImageType::PointType min=boundingBox->GetMinimum();
 IM->SetOrigin(min);

 if(this->GetremoveFrequencyOversampling() && this->canRemoveFrequencyOS()){
	 typename VectorImageType::PointType oVSP;

	 std::cout<<"remove OS origin was"<<min<<std::endl;


	 typename VectorImageType::IndexType ind={floor(SIZE[0]/2),0,0};
	 IM->TransformIndexToPhysicalPoint(ind,oVSP);


	 for (auto n=0;n<3;n++){
		 min[n]=min[n]+(std::sqrt(std::pow((oVSP[n]-min[n]),2)));

	 }



	 IM->SetOrigin(min);
 };


 return IM;
 };

template< class TImage>
typename TImage::Pointer  ISMRMRDToITKImageFilter< TImage>
::tridimensionalRead()
 {	std::cout<<"tridimensinoal"<<std::endl;
 //the acqisition
 ISMRMRD::Acquisition acq;
 //must be instantiated with a filename
 ISMRMRD::Dataset ismrmdataset(this->GetIsmrmrdFileName().c_str(),"dataset", false);

 //create the image vector
 typename VectorImageType::SizeType SIZE;
 VectorImagePointerType	IM;

 if(this->GetremoveFrequencyOversampling()){
	 IM=this->getNewReconVectorImage();
	 SIZE= this->getReconSize();
 }else{
	 IM=this->getNewEncodedVectorImage();
	 SIZE= this->getEncodedSize();
 }


IM->Allocate();

 //ORIGIN is a point cloud i've to find the minimum later
 using CoordType = float;
 const unsigned int Dimension = VectorImageType::ImageDimension;


 using PointSetType = itk::PointSet< CoordType, Dimension >;

 using PointIdentifier = typename PointSetType::PointIdentifier;
 using PointType = typename PointSetType::PointType;
 using PointsContainerPointer = typename PointSetType::PointsContainerPointer;
 typename PointSetType::Pointer   pointSet = PointSetType::New();
 PointsContainerPointer  points = pointSet->GetPoints();

 unsigned long int pointId = 0;


 typename VectorImageType::DirectionType DIRECTIONS;

 unsigned short int numberofcoil=this->getNumberOfCoils();
 typename VectorImageType::IndexType index;
 unsigned short int coilnumber;
 typename VectorImageType::PixelType   vectorValue(numberofcoil);

 long int distance;

 ISMRMRD::ISMRMRD_EncodingCounters counter;




 //just ask once:)
 int osfactor=1;
 if(this->GetremoveFrequencyOversampling()){
	 osfactor=2;
 };



 VectorImageInnerPixelType NR=1;
 if(this->GetaverageImageThroughAverages()){
	 NR=(VectorImageInnerPixelType)((float)this->getNumberOfAverage());
 };

 std::cout<<"number of averages used:"<<NR<<std::endl;
 std::cout<<"number of averages in the image:"<<this->getNumberOfAverage()<<std::endl;

 for (auto p=0; p<ismrmdataset.getNumberOfAcquisitions();p++)
 {
	 ismrmdataset.readAcquisition(p, acq);
	 this->AcquisitionInfo(acq);
	 if (this->isTheAcquisitionRequestedImage(acq) && !(this->isTheAcquisitionNoise(acq)))
	 {
		 counter=acq.idx();
		 //since it is bidimensional acquisition the 1nd index represent the frequency encoding and will change in the for loop
		 index={0,counter.kspace_encode_step_1,counter.kspace_encode_step_2};
		 //for all the frequency encoded

		 // std::cout<<index<<std::endl;

		 for(auto idx=0;idx<SIZE[0]*osfactor;idx=idx+osfactor){
			 index[0]=idx/osfactor;
			 vectorValue=IM->GetPixel(   index  );
			 //				 if (index[0]==0 && index[1]==0 && index[2]==0 ){
			 //				 std::cout<<p<<","<<IM->GetPixel(   index  )<<std::endl;
			 //				 std::cout<<p<<","<< NR<<"\n"<<std::endl;
			 //				 std::cout<<p<<","<< acq.data(idx,1)<<"\n"<<std::endl;
			 //
			 //
			 //				 }
			 for(auto ch=0;ch<numberofcoil;ch++){
				 vectorValue[ch] +=  (acq.data(idx,ch)/NR);
			 }
			 IM->SetPixel(   index,   vectorValue  );

		 }

		 DIRECTIONS=this->getAcquisitionDirections(acq);

		 points->InsertElement( pointId, this->getAcquisitionOrigin(acq) );
		 pointId++;
	 }
 };


 //set the directions with the last acq
 IM->SetDirection(DIRECTIONS);

 using BoundingBoxType = itk::BoundingBox< PointIdentifier, Dimension, CoordType >;
 typename  BoundingBoxType::Pointer boundingBox = BoundingBoxType::New();
 boundingBox->SetPoints(points);
 boundingBox->ComputeBoundingBox();


 typename VectorImageType::PointType min=boundingBox->GetMinimum();
 IM->SetOrigin(min);

 if(this->GetremoveFrequencyOversampling()){

	 std::cout<<min<<std::endl;

	 typename VectorImageType::PointType oVSP;

	 typename VectorImageType::IndexType ind={floor(SIZE[0]/2),0,0};
	 IM->TransformIndexToPhysicalPoint(ind,oVSP);
	 min[0]=min[0]+std::sqrt(std::pow((oVSP[0]-min[0]),2));
	 std::cout<<min<<std::endl;

	 IM->SetOrigin(min);
 };


 return IM;
 };

template< class TImage>
void ISMRMRDToITKImageFilter< TImage>
::GenerateData()
 {



	switch(m_KSpaceDimension){
	case cm::KSpaceAcquisitionDimension::BIDIMENSIONAL:
		//		if (this->GetaverageImageThroughAverages()){
		//	//		this->GetOutput()->Graft(this->bidimensionalReadWithAverages());
		//		}else{
		//			this->GetOutput()->Graft(this->bidimensionalRead());
		//		};

		this->GetOutput()->Graft(this->bidimensionalRead());

		break;
	case cm::KSpaceAcquisitionDimension::TRIDIMENSIONAL:

		typename TImage::Pointer IM=this->tridimensionalRead();
//		std::cout<<IM;
//		IM->GetPixel({{20,20,20}});
//
//		typedef typename cm::SumOfSquareImageFilter<TImage,itk::Image<vcl_complex<float>,3>> SOS;
//		typename SOS::Pointer sos= SOS::New();
//		sos->SetInput(IM);
//		sos->Update();
//		sos->SetKSpaceDimension(cm::KSpaceAcquisitionDimension::TRIDIMENSIONAL);



		this->GetOutput()->Graft(IM);
		break;
	};


 }

}// end namespace


#endif



//acq.
//
//	ismrmdataset.readAcquisition(0, acq);
//
//	this->AcquisitionInfo(acq);
//
//
//
//	uint16_t nCoils = acq.active_channels();
//	std::cout<<"N data element:"<<acq.getNumberOfDataElements()<<std::endl;
//	std::cout<<"available chan:" <<acq.available_channels()<<std::endl;
//
//	typename VectorImageType::PointType position;
//	position=*acq.position();
//	std::cout<<"position:" <<position<<std::endl;
//
//	std::cout<<"read_dir:" <<*acq.read_dir()<<std::endl;
//	//		std::cout<<"N" <<acq.getNumberOfDataElements()/acq.available_channels()<<"must be equal to "<<Fr<<std::endl;
//
//
//
//
//	//that's the frequency vector of std::complex
//	std::vector<VectorImageInnerPixelType>KS(Fr);
//	//And its iterator
//	typename std::vector<VectorImageInnerPixelType>::iterator L;
//
//
//	// Pack it's an array of 2d images of std complex (one per coil)
//	std::vector<ScalarImagePointerType>Pack;
//	//this is asingle coil image
//	typename ScalarImageType::IndexType P;
////
//
//	std::cout<<(unsigned int)l.slice->maximum+1 << "\n\n";
//
//
//	//** that's the encoded Kspace*/
//	typename ScalarImageType::SizeType size={e_space.matrixSize.x,e_space.matrixSize.y,(unsigned int)l.slice->maximum+1};
//	std::cout<<size.GetElement(2);
//	typename ScalarImageType::IndexType start={0,0,0};
//
//	typename ScalarImageType::RegionType region;
//	region.SetSize( size );
//	region.SetIndex( start );
//	float s0=e_space.fieldOfView_mm.x/(float)e_space.matrixSize.x;
//	float s1=e_space.fieldOfView_mm.y/(float)e_space.matrixSize.y;
//	float s2=e_space.fieldOfView_mm.z/(float)e_space.matrixSize.z;
//
//	typename ScalarImageType::SpacingType spacing;
//	spacing[0]=s0;
//	spacing[1]=s1;
//	spacing[2]=s2;
//
//
//	/* channel*/
//	for (uint16_t c=0; c<nCoils; c++)
//	{
//		//std::cout<<"coil"<<c<<std::endl;
//		ScalarImagePointerType coil = ScalarImageType::New();
//
//		coil->SetRegions( region );
//		coil->Allocate();
//		coil->SetSpacing( spacing );
//
//		Pack.push_back(coil);
//
//		for (auto p=0; p<ismrmdataset.getNumberOfAcquisitions();p++)
//		{
//
//
//			ismrmdataset.readAcquisition(p, acq);
//
//
//
//
//			memcpy(&KS[0],&acq.data(0, c), Fr*sizeof(ScalarImageInnerPixelType));
//
//			P[1]=acq.idx().kspace_encode_step_1;
//			P[2]=acq.idx().slice;
//			for(L=KS.begin();L!=KS.end();L++)
//			{
//
//				P[0]=std::distance(KS.begin(), L);
//				Pack[c]->SetPixel(P,*L);
//				//std::cout<<P<<" value "<< *L<<std::endl;
//			}
//
//
//
//		}
//
//
//	}
//
//
//
//
//
//VectorImagePointerType IM=VectorImageType::New();
//
//IM=composeVectorImage<ScalarImageType,VectorImageType>(Pack);
//
//
//this->GetOutput()->Graft(IM);


