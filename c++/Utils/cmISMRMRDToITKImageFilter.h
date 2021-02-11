/*
 * ISMRMRDToITKFilter.h
 *
 *  Created on: Jan 31, 2019
 *      Author: montie01
 */



#ifndef __cmISMRMRDToITKImageFilter_h
#define __cmISMRMRDToITKImageFilter_h

#pragma once
#include "cm.h"
#include "itkImageSource.h"

#include <set>
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/dataset.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/version.h"


namespace cm
{

/** \class ISMRMRDToITKImageFilter
 * \brief This class read an <a href="http://ismrmrd.github.io/#c-example-applications">ISRMRMRD v1</a> file and create a 3D kspace image.
 *
 * \author Eng. Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 */


template< class TImage>
class ISMRMRDToITKImageFilter:public itk::ImageSource< TImage >
{
public:

	/* Standard class typedefs. */
	typedef ISMRMRDToITKImageFilter         Self;
	typedef itk::ImageSource< TImage > Superclass;
	typedef itk::SmartPointer< Self >  Pointer;

	/* Method for creation through the object factory. */
	itkNewMacro(Self);

	/* Run-time type information (and related methods). */
	itkTypeMacro(ISMRMRDToITKImageFilter, ImageSource);



	/*  the main image type*/
	typedef TImage                   VectorImageType;
	typedef typename TImage::Pointer VectorImagePointerType;
	typedef typename TImage::PixelType VectorImagePixelType;
	typedef typename TImage::InternalPixelType VectorImageInnerPixelType;
	typedef typename VectorImageType::InternalPixelType::value_type VectorImageComplexInternalPixelType;// (float or double)


	/* the main image type*/
	typedef typename itk::Image<VectorImageInnerPixelType, VectorImageType::ImageDimension >                   ScalarImageType;
	typedef typename ScalarImageType::Pointer ScalarImagePointerType;
	typedef typename ScalarImageType::PixelType ScalarImagePixelType;
	typedef typename ScalarImageType::InternalPixelType ScalarImageInnerPixelType;





	itkGetMacro(IsmrmrdFileName,std::string);
	//this function is created directly from the class
	//itkSetMacro(IsmrmrdFileName,std::string);

	itkGetMacro(removeFrequencyOversampling,bool);
	itkSetMacro(removeFrequencyOversampling,bool);

	itkGetMacro(averageImageThroughAverages,bool);
	itkSetMacro(averageImageThroughAverages,bool);


	itkGetMacro(requestedEncoding,uint16_t);
	itkSetMacro(requestedEncoding,uint16_t);

	itkGetMacro(requestedAverage,uint16_t);
	itkSetMacro(requestedAverage,uint16_t);
	
	




	itkGetMacro(requestedContrast,uint16_t);
	itkSetMacro(requestedContrast,uint16_t);

	itkGetMacro(requestedRepetition,uint16_t);
	itkSetMacro(requestedRepetition,uint16_t);

	itkGetMacro(requestedSet,uint16_t);
	itkSetMacro(requestedSet,uint16_t);

	itkGetMacro(requestedSegment,uint16_t);
	itkSetMacro(requestedSegment,uint16_t);


   	   //Kspacedimension
    	itkGetMacro(KSpaceDimension,cm::KSpaceAcquisitionDimension);
    	itkSetMacro(KSpaceDimension,cm::KSpaceAcquisitionDimension);

	
    std::vector<float> getSliceOrder();

	VectorImagePointerType  bidimensionalRead();

	VectorImagePointerType  tridimensionalRead();

	std::vector<VectorImagePointerType>  getReplicas(std::string W);

	void AcquisitionInfo(ISMRMRD::Acquisition acq);


	long int getResonanceFrequency(){
		return m_hdr.experimentalConditions.H1resonanceFrequency_Hz;
	};


	uint16_t getNumberOfRepetition(){
		return m_encodingLimits.repetition->maximum+1;
	};

	uint16_t getNumberOfAverage(){
		return m_encodingLimits.average->maximum+1;
	};

	uint16_t getNumberOfContrast(){
		return m_encodingLimits.contrast->maximum+1;
	};

	uint16_t getNumberOfSegment(){
		return m_encodingLimits.segment->maximum+1;
	};

	uint16_t getNumberOfSet(){
		return m_encodingLimits.set->maximum+1;
	};

	uint16_t getNumberOfSlice(){
		return m_encodingLimits.slice->maximum+1;
	};

	uint16_t getNumberOfPhase(){

		return m_encodingLimits.phase->maximum+1;
	};

	uint16_t getNumberOfencodingstep0(){
		return m_encodingLimits.kspace_encoding_step_0->maximum+1;
	};

	uint16_t getNumberOfencodingstep1(){
		return m_encodingLimits.kspace_encoding_step_1->maximum+1;
	};
	uint16_t getNumberOfencodingstep2(){
		return m_encodingLimits.kspace_encoding_step_2->maximum+1;
	};

	uint16_t getNumberOfCoils(){
		return m_hdr.acquisitionSystemInformation.get().receiverChannels.get();
	};

	float getFieldStrength(){

		return m_hdr.acquisitionSystemInformation.get().systemFieldStrength_T.get();
	};


	typename VectorImageType::DirectionType getAcquisitionDirections(ISMRMRD::Acquisition acq){
		typename VectorImageType::DirectionType direction;
		direction.SetIdentity();

		//typename VectorImageType::PointType position,readdir,phasedir,slicedir,bed;

		float *readdir;
		readdir=acq.read_dir();
		for(auto a=0;a<3;a++){
			direction[0][a]=readdir[a];
		}

		float *phasedir;
		phasedir=acq.phase_dir();
		for(auto a=0;a<3;a++){
			direction[1][a]=phasedir[a];
		}

		float *slicedir;
		slicedir=acq.slice_dir();
		for(auto a=0;a<3;a++){
			direction[2][a]=slicedir[a];
		}

		return direction;





	};

	bool canRemoveFrequencyOS(){
		bool O=false;

		typename VectorImageType::SizeType R= this->getReconSize();
		typename VectorImageType::SizeType E= this->getEncodedSize();

		//if the x component of reconsize*2 is equal to encodesize it means i can remove the os
		if (std::floor(E[0]/2)==R[0])
		{
			O=true;
		}else{
			O=false;
		}
return O;

	}

	typename VectorImageType::PointType getAcquisitionOrigin(ISMRMRD::Acquisition acq){
		typename VectorImageType::PointType point;


		//typename VectorImageType::PointType position,readdir,phasedir,slicedir,bed;

		float *position;
		position=acq.position();
		for(auto a=0;a<3;a++){
			point[a]=position[a];
		}


		return point;
	};

	typename VectorImageType::SizeType getEncodedSize(){

		switch(m_KSpaceDimension){
		case cm::KSpaceAcquisitionDimension::BIDIMENSIONAL:
			return {m_encodedSpace.matrixSize.x,m_encodedSpace.matrixSize.y,(unsigned int)this->getNumberOfSlice()};
		case cm::KSpaceAcquisitionDimension::TRIDIMENSIONAL:
			return {m_encodedSpace.matrixSize.x,m_encodedSpace.matrixSize.y,m_encodedSpace.matrixSize.z};
		default: return {0,0,0};
		}
	};

	typename VectorImageType::SizeType getReconSize(){

		switch(m_KSpaceDimension){
		case cm::KSpaceAcquisitionDimension::BIDIMENSIONAL:
			return {m_reconSpace.matrixSize.x,m_reconSpace.matrixSize.y,(unsigned int)this->getNumberOfSlice()};
		case cm::KSpaceAcquisitionDimension::TRIDIMENSIONAL:
			return {m_reconSpace.matrixSize.x,m_reconSpace.matrixSize.y,m_reconSpace.matrixSize.z};
		default: return {0,0,0};
		}
	};


	std::vector<float> getEncodedFOV(){

		switch(m_KSpaceDimension){
		case cm::KSpaceAcquisitionDimension::BIDIMENSIONAL:
			return {m_encodedSpace.fieldOfView_mm.x,m_encodedSpace.fieldOfView_mm.y,(unsigned int)this->getNumberOfSlice()*m_encodedSpace.fieldOfView_mm.z};
		case cm::KSpaceAcquisitionDimension::TRIDIMENSIONAL:
			return {m_encodedSpace.fieldOfView_mm.x,m_encodedSpace.fieldOfView_mm.y,m_encodedSpace.fieldOfView_mm.z};
		default: return {0,0,0};
		}
	};

	std::vector<float> getReconFOV(){

		switch(m_KSpaceDimension){
		case cm::KSpaceAcquisitionDimension::BIDIMENSIONAL:
			return {m_reconSpace.fieldOfView_mm.x,m_reconSpace.fieldOfView_mm.y,(unsigned int)this->getNumberOfSlice()*m_reconSpace.fieldOfView_mm.z};
		case cm::KSpaceAcquisitionDimension::TRIDIMENSIONAL:
			return {m_reconSpace.fieldOfView_mm.x,m_reconSpace.fieldOfView_mm.y,m_reconSpace.fieldOfView_mm.z};
		default: return {0,0,0};
		}
	};



	typename VectorImageType::SpacingType getEncodedResolution(){
		std::vector<float> FOV=this->getEncodedFOV();
		typename VectorImageType::SizeType MATRIX=this->getEncodedSize();

		typename VectorImageType::SpacingType SP;
		SP.Fill(0);

		for (auto o=0;o<FOV.size();o++){
			//			SP[o]=(VectorImageType::SpacingType::ValueType)FOV[o]/(VectorImageType::SpacingType::ValueType)MATRIX[o];
			SP[o]=(double)FOV[o]/(double)MATRIX[o];

		}
		return SP;
	};

	typename VectorImageType::SpacingType getReconResolution(){

		std::vector<float> FOV=this->getReconFOV();
		typename VectorImageType::SizeType MATRIX=this->getReconSize();

		typename VectorImageType::SpacingType SP;
		SP.Fill(0);

		for (auto o=0;o<FOV.size();o++){

			SP[o]=(double)FOV[o]/(double)MATRIX[o];
		}
		return SP;

	};


	VectorImagePointerType getNewEncodedVectorImage(){

		VectorImagePointerType IM=VectorImageType::New();

		typename VectorImageType::IndexType origin={0,0,0};

		typename VectorImageType::RegionType region;
		region.SetSize( this->getEncodedSize() );
		region.SetIndex( origin );
		IM->SetRegions( region );
		IM->SetSpacing( this->getEncodedResolution() );
		IM->SetNumberOfComponentsPerPixel(this->getNumberOfCoils());
		IM->Allocate();

		return IM;


	};


	VectorImagePointerType getNewReconVectorImage(){

		VectorImagePointerType IM=VectorImageType::New();

		typename VectorImageType::IndexType origin={0,0,0};

		typename VectorImageType::RegionType region;
		region.SetSize( this->getReconSize() );
		region.SetIndex( origin );
		IM->SetRegions( region );
		IM->SetSpacing( this->getReconResolution() );
		IM->SetNumberOfComponentsPerPixel(this->getNumberOfCoils());
		IM->Allocate();

		return IM;


	};

	bool isTheAcquisitionRequestedImage(ISMRMRD::Acquisition acq){

		ISMRMRD::ISMRMRD_EncodingCounters counter;
		counter=acq.idx();

		bool ans=false;

		if(this->GetaverageImageThroughAverages()){

			if(
					//(counter.average==this->GetrequestedAverage()) &&
					(counter.repetition==this->GetrequestedRepetition()) &&
					(counter.set==this->GetrequestedSet()) &&
					(counter.segment==this->GetrequestedSegment())
					){
				ans=true;
					}


		}else{

			if(
					(counter.contrast==this->GetrequestedContrast()) &&
					(counter.average==this->GetrequestedAverage()) &&
					(counter.repetition==this->GetrequestedRepetition()) &&
					(counter.set==this->GetrequestedSet()) &&
					(counter.segment==this->GetrequestedSegment())
			){
				ans=true;

			}
		}

		return ans;
	};




	bool isTheAcquisitionRequestedImageReplicas(ISMRMRD::Acquisition acq){

		ISMRMRD::ISMRMRD_EncodingCounters counter;
		counter=acq.idx();

		bool ans=false;

		if(this->GetaverageImageThroughAverages()){

			if(
					//(counter.average==this->GetrequestedAverage()) &&
					//(counter.repetition==this->GetrequestedRepetition()) &&
					(counter.set==this->GetrequestedSet()) &&
					(counter.segment==this->GetrequestedSegment())
					){
				ans=true;
					}


		}else{

			if(
					(counter.contrast==this->GetrequestedContrast()) &&
					(counter.average==this->GetrequestedAverage()) &&
					//(counter.repetition==this->GetrequestedRepetition()) &&
					(counter.set==this->GetrequestedSet()) &&
					(counter.segment==this->GetrequestedSegment())
			){
				ans=true;

			}
		}

		return ans;
	};




	bool isTheAcquisitionNoise(ISMRMRD::Acquisition acq){

		ISMRMRD::AcquisitionHeader H;

		H=acq.getHead();
		return H.isFlagSet(ISMRMRD::ISMRMRD_AcquisitionFlags::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);

	}






	bool isAccelerated(){

		if (this->getAcceleration1()>1 || this->getAcceleration2()>1 ){
			return true;
		}else{
			return false;
		}
	};

	unsigned short int getAcceleration1(){

		return m_hdr.encoding[m_requestedEncoding].parallelImaging.get().accelerationFactor.kspace_encoding_step_1;
	};

	unsigned short int getAcceleration2(){

		return m_hdr.encoding[m_requestedEncoding].parallelImaging.get().accelerationFactor.kspace_encoding_step_2;
	};

	unsigned short int getCalibrationMode(){

		return m_hdr.encoding[m_requestedEncoding].parallelImaging.get().calibrationMode.get();
	};


	std::string getTrajectory(){


		switch(m_hdr.encoding[m_requestedEncoding].trajectory){
		case ISMRMRD::TrajectoryType::CARTESIAN: return "Cartesian";
		case ISMRMRD::TrajectoryType::EPI: return "Epi";
		case ISMRMRD::TrajectoryType::GOLDENANGLE: return "Goldenangle";
		case ISMRMRD::TrajectoryType::RADIAL: return "Radial";
		case ISMRMRD::TrajectoryType::SPIRAL: return "Spiral";
		default: return "NO";
		};



	};


	std::string getWaveForm(){


		if(m_hdr.waveformInformation.size()>0){

			switch(m_hdr.waveformInformation[m_requestedEncoding].waveformType){
			case ISMRMRD::WaveformType::ECG: return "Ecg";
			case ISMRMRD::WaveformType::GRADIENTWAVEFORM: return "GradientWaveForm";
			case ISMRMRD::WaveformType::PULSE: return "Pulse";
			case ISMRMRD::WaveformType::RESPIRATORY: return "Respiratory";
			case ISMRMRD::WaveformType::TRIGGER: return "Trigger";
			default: return "NO";
			};
		}else{
			return "no wave info";
		}


	};




	std::string getKSpaceDimensionString(){
		switch(m_KSpaceDimension){
		case cm::KSpaceAcquisitionDimension::BIDIMENSIONAL: return "Bidimensional";
		case cm::KSpaceAcquisitionDimension::TRIDIMENSIONAL: return "Tridimensional";
		default: return "NO";
		};


	};





	float getRelativeReceiverBandwidth(){

		return m_hdr.acquisitionSystemInformation.get().relativeReceiverNoiseBandwidth.get();
	};


	cm::KSpaceAcquisitionDimension getKspaceType(){
		if(this->getNumberOfencodingstep2()>1 && this->getNumberOfSlice()==1){
			return cm::KSpaceAcquisitionDimension::TRIDIMENSIONAL;
		}else{
			return cm::KSpaceAcquisitionDimension::BIDIMENSIONAL;};
	};


	void SetIsmrmrdFileName(const std::string& ismrmrdFileName)
	{
		//set the filename
		m_IsmrmrdFileName = ismrmrdFileName;
		//instantiate only with a filename
		ISMRMRD::Dataset m_ismrmrdDataset(this->GetIsmrmrdFileName().c_str(),"dataset", false);

		m_ismrmrdDataset.readHeader(m_xml);
		ISMRMRD::deserialize(m_xml.c_str(),m_hdr);

		std::cout<<"\nnumber of encoding: "<<m_hdr.encoding.size()<<std::endl;

		m_encodedSpace = m_hdr.encoding[m_requestedEncoding].encodedSpace;
		m_reconSpace=m_hdr.encoding[m_requestedEncoding].reconSpace;
		m_encodingLimits = m_hdr.encoding[m_requestedEncoding].encodingLimits;



		m_KSpaceDimension=this->getKspaceType();



	}

	const ISMRMRD::IsmrmrdHeader& getHdr() const
	{
		return m_hdr;
	}

	const std::string& getXml() const
	{
		return m_xml;
	}


	//getEncodedMatrixSlice();

protected:
	//constructor
	ISMRMRDToITKImageFilter(){
		//for multiraid
		m_requestedEncoding=0;
		
		m_requestedAverage=0;

		m_requestedContrast=0;

		m_requestedRepetition=0;

		m_requestedSet=0;

		m_requestedSegment=0;

		m_removeFrequencyOversampling=false;

		m_averageImageThroughAverages=false;

	}


	~ISMRMRDToITKImageFilter(){}

	/* Does the real work. */
	virtual void GenerateData();

private:
	ISMRMRDToITKImageFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented




	// this is the xml file that u can see in hdfview
	std::string m_xml;

	// header in the file
	ISMRMRD::IsmrmrdHeader m_hdr;

	//Encodingspace
	ISMRMRD::EncodingSpace m_encodedSpace;
	//Reconstruction Space
	ISMRMRD::EncodingSpace m_reconSpace;


	/* EncodingLimits
	Optional< Limit > 	kspace_encoding_step_0
	Optional< Limit > 	kspace_encoding_step_1
	Optional< Limit > 	kspace_encoding_step_2
	Optional< Limit > 	average
	Optional< Limit > 	slice
	Optional< Limit > 	contrast
	Optional< Limit > 	phase
	Optional< Limit > 	repetition
	Optional< Limit > 	set
	Optional< Limit > 	segment*/
	ISMRMRD::EncodingLimits m_encodingLimits;



	std::string m_IsmrmrdFileName;

	bool m_removeFrequencyOversampling;
	bool m_averageImageThroughAverages;
	
	//for multiraid
	uint16_t    m_requestedEncoding;
		
	uint16_t 	m_requestedAverage;

	uint16_t 	m_requestedContrast;

	uint16_t 	m_requestedRepetition;

	uint16_t 	m_requestedSet;

	uint16_t    m_requestedSegment;

	cm::KSpaceAcquisitionDimension m_KSpaceDimension;

};
}; //namespace ITK

#ifndef ITK_MANUAL_INSTANTIATION
#include "cmISMRMRDToITKImageFilter.hxx"
#endif

#endif // __cmISMRMRDToITKImageFilter_h
