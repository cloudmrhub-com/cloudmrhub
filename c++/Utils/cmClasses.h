#pragma once

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include <iostream>
#include <vector>
#include <functional>
#include <future>
#include <functional>
#include <type_traits>
#include "itkCastImageFilter.h"


//OR
#include "itkOrImageFilter.h"

//MINMAX
#include "itkMinimumMaximumImageCalculator.h"
//SUBTRACT,SNR
#include "itkSubtractImageFilter.h"

// multiply
#include "itkMultiplyImageFilter.h"

//json
#include "json.hpp"


#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToPhaseImageFilter.h"

#include "itkImage.h"

#include "itkVector.h"
#include "itkVariableLengthVector.h"

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>

#include "ImageUtils.h"

#include "cm.h"
#include "mropt.h"



template<class TImage, class TImage_ >
typename TImage::Pointer differenceImageSNR(typename TImage::Pointer r1, typename TImage_::Pointer r2) {

	//report images in the space non so se serve but:)
	typename TImage::Pointer im2on1 = TImage::New();
	im2on1 = roiOnImage<TImage, TImage_>(r1, r2);

	im2on1=subtractImage<TImage>(r1,im2on1);
	im2on1=multiplyImageTimesScalar<TImage>(im2on1,1.0/sqrt(2));

	std::cout << " SNR--" << std::endl;

	return im2on1;
}




template <class TI, class IT>
nlohmann::json sliceWriterjson( typename TI::Pointer cout){
	// the outputimage will be always a float one since i subdivide in real and imaginary so i don't have problem



	nlohmann::json slice;



	slice["w"]=cout->GetLargestPossibleRegion().GetSize()[1];
	slice["h"]=cout->GetLargestPossibleRegion().GetSize()[0];

	slice["type"]= "double complex";


	typedef itk::ComplexToRealImageFilter <TI,IT>itkComplexToRealImageFilterIT;
	typename itkComplexToRealImageFilterIT::Pointer reFilter= itkComplexToRealImageFilterIT::New ();
	reFilter->SetInput(cout);
	reFilter->Update();

	slice["Vr"] = imageToArray<IT>(reFilter->GetOutput());




	typedef itk::ComplexToImaginaryImageFilter <TI,IT>CTIF;
	typename CTIF::Pointer iFilter= CTIF::New ();
	iFilter->SetInput(cout);
	iFilter->Update();

	slice["Vi"] = imageToArray<IT>(iFilter->GetOutput());





	return slice;


}



template <class TI>
nlohmann::json sliceWriterjson( vnl_matrix<std::complex<TI>> &X){
	// the outputimage will be always a float one since i subdivide in real and imaginary so i don't have problem



	nlohmann::json slice;



	slice["w"]=X.rows();
	slice["h"]=X.cols();

	slice["type"]= "double complex";




	std::vector<TI>R;//(X.rows()*X.cols());
	std::vector<TI>I;//(X.rows()*X.cols());


	for (auto y=0;y<X.cols();y++)
	{
		for (auto x=0;x<X.rows();x++)
		{

			R.push_back(std::real(X(x,y)));
			I.push_back(std::imag(X(x,y)));


		}

	}

	slice["Vr"] = R;

	slice["Vi"] = I;





	return slice;


}




template <class TI>
nlohmann::json sliceWriterjson( vnl_matrix<std::complex<TI>> *X){
	// the outputimage will be always a float one since i subdivide in real and imaginary so i don't have problem



	nlohmann::json slice;



	slice["w"]=X->rows();
	slice["h"]=X->cols();

	slice["type"]= "double complex";




	std::vector<TI>R;//(X.rows()*X.cols());
	std::vector<TI>I;//(X.rows()*X.cols());


	for (auto x=0;x<X->cols();x++)
	{
		for (auto y=0;y<X->rows();y++)
		{

			R.push_back(std::real(X(x,y)));
			I.push_back(std::imag(X(x,y)));


		}

	}

	slice["Vr"] = R;

	slice["Vi"] = I;
    




	return slice;


}




//	  noisecovbis = zeros(nchan);
//	            for iCh = 1:nchan
//	                for jCh = 1:nchan
//	                    noisecovbis(iCh,jCh)=sum(sum(noise(:,:,iCh).*conj(noise(:,:,jCh))))/(size(noise,1)*size(noise,2));
//	                end
//	            end

template <class Pxl>
Pxl covarianceRiccardo(std::vector<Pxl> & X,std::vector<Pxl> & Y){
	Pxl O=0.0;
	typename std::vector<Pxl>::iterator it;
	typename std::vector<Pxl>::iterator ot;

	ot=Y.begin();
	it= X.begin();
	Pxl N=(Pxl) X.size();

	while ( it != X.end()){
		O+=*it * std::conj( *ot);
		it++;
		ot++;
	};


	return (Pxl) O/N;

};
















#pragma once
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/dataset.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/version.h"



template <class TI,class VI>
typename VI::Pointer ISMRMRDToImage(std::string IN){

	ISMRMRD::Dataset d(IN.c_str(),"dataset", false);
	std::cout<<"N acquisition: "<<d.getNumberOfAcquisitions()<<std::endl;
	std::cout<<"N Waveform: "<<d.getNumberOfWaveforms()<<std::endl;

	std::string xml;
	d.readHeader(xml);
	ISMRMRD::IsmrmrdHeader hdr;
	ISMRMRD::deserialize(xml.c_str(),hdr);
	/** this is the xml file that i can see in hdfview*/
	//std::cout<<xml<<std::endl;
	std::cout << "XML Header version: " << hdr.version << std::endl;
	std::cout << "Number of encoding spaces: " << hdr.encoding.size() << std::endl;
	ISMRMRD::EncodingSpace e_space = hdr.encoding[0].encodedSpace;

	std::cout<<"ENCODING:"<<std::endl;
	std::cout << "Matrix size: " << e_space.matrixSize.x<<", "<< e_space.matrixSize.y<<", "<< e_space.matrixSize.z << std::endl;
	std::cout << "Fov mm: " << e_space.fieldOfView_mm.x<<", "<< e_space.fieldOfView_mm.y<<", "<< e_space.fieldOfView_mm.z << std::endl;


	ISMRMRD::EncodingLimits l = hdr.encoding[0].encodingLimits;

	ISMRMRD::Optional<ISMRMRD::Limit>QQ;

	std::cout<<"limits:"<<std::endl;
	std::cout << "Repetition: " << l.repetition<<std::endl;
	std::cout << "AVG: " << l.average<<std::endl;
	QQ=l.phase.get();
	std::cout << "Phase: " << QQ<<std::endl;


	ISMRMRD::EncodingSpace r_space=hdr.encoding[0].reconSpace;
	std::cout<<"RECON:"<<std::endl;
	std::cout << "Matrix size: " <<r_space.matrixSize.x <<", "<< r_space.matrixSize.y<<", "<< r_space.matrixSize.z << std::endl;
	std::cout << "Fov mm: " << r_space.fieldOfView_mm.x<<", "<< r_space.fieldOfView_mm.y<<", "<< r_space.fieldOfView_mm.z << std::endl;

	uint16_t Fr = e_space.matrixSize.x;
	uint16_t Ph = e_space.matrixSize.y;


	//typedef itk::Image<PixelType, 2 > ScalarImageType;

	ISMRMRD::Acquisition acq;

	d.readAcquisition(0, acq);
	uint16_t nCoils = acq.active_channels();
	std::cout<<"N data element"<<acq.getNumberOfDataElements()<<std::endl;
	std::cout<<"av chan" <<acq.available_channels()<<std::endl;
	std::cout<<"N" <<acq.getNumberOfDataElements()/acq.available_channels()<<"must be equal to "<<Fr<<std::endl;


	std::vector<std::complex<float>>KS(Fr);
	std::vector<std::complex<float>>::iterator L;

	std::vector<typename TI::Pointer>Pack;
	typename TI::IndexType P;


	typename TI::SizeType size={e_space.matrixSize.x,e_space.matrixSize.y};
		typename TI::IndexType start={0,0};

		typename TI::RegionType region;
		region.SetSize( size );
		region.SetIndex( start );



		float s0=r_space.fieldOfView_mm.x/(float)r_space.matrixSize.x;
		float s1=r_space.fieldOfView_mm.y/(float)r_space.matrixSize.y;

		typename TI::SpacingType spacing;
		spacing[0]=s0;
		spacing[1]=s1;






	/** channel*/
	for (uint16_t c=0; c<nCoils; c++)
	{
		//std::cout<<"coil"<<c<<std::endl;
		typename TI::Pointer coil = TI::New();

		coil->SetRegions( region );
		coil->Allocate();
		coil->SetSpacing( spacing );

		Pack.push_back(coil);

		for (auto p=0; p<d.getNumberOfAcquisitions();p++)
		{

			d.readAcquisition(p, acq);
			memcpy(&KS[0],&acq.data(0, c), Fr*sizeof(std::complex<float>));

			P[1]=acq.idx().kspace_encode_step_1;

			for(L=KS.begin();L!=KS.end();L++)
			{

				P[0]=std::distance(KS.begin(), L);
				Pack[c]->SetPixel(P,*L);
				//				std::cout<<P<<" value "<< *L<<std::endl;
			}


		}


	}

	typename VI::Pointer IM=VI::New();

	IM=composeVectorImage<TI,VI>(Pack);


	return IM;
}






