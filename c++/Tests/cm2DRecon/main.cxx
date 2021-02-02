#include "vnl/vnl_matrix.h"

#include <boost/algorithm/string.hpp>

#include <boost/make_shared.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/locale/encoding_utf.hpp>
using boost::locale::conv::utf_to_utf;


//Filetype
#include "itkVectorImage.h"


#include <iostream>
#include <fstream>

//IO
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

#include "../../Module/cmVectorImageReconstructor.h"

#include "itkImageRegionConstIterator.h"

int main(int argc, char *argv[]) {

	/** Signal file name. */
	std::string SFN;
	/** Noise file name. */
	std::string NFN;
	/** Noise file name. */
	std::string OFN;



	/** The options */
	po::options_description desc("The software calculates 2D SNR from two vector iamges file, noise and signal: Allowed options");
	desc.add_options()("help,h", "Produce HELP message")
															("signalfilename,s",po::value<std::string>(&SFN), " Signal file name 3D Vector image complex kSpace data")
															("noisefilename,n",po::value<std::string>(&NFN), " Noise file name 3D Vector image complex kSpace data")

															("outputfilename,o",po::value<std::string>(&OFN), " Reconstructed Image file name a 3D scalar Image")
															;


	po::options_description myerror("RSS Reconstruction through RSS\n Eros Montin, PhD");

	po::variables_map vm;

	if (argc < 1) {
		std::cout << myerror << "\n";
		return 1;
	}
	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			std::cout << desc << "\n";
			return 1;
		}

	}

	catch (po::error& e) {
		std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
		std::cerr << desc << std::endl;
		return -1;

	}




	using PT = float;
	using PixelType = std::complex< PT >;
//	using PixelType = vcl_complex< PT >;
	const unsigned int      	Dimension = 3;

	typedef itk::VectorImage<PixelType, Dimension >  VectorImageType;
	typedef itk::Image<PixelType, Dimension >    ScalarImageType;
	typedef itk::Image<PT, Dimension >    ResultsImageType;

	typedef itk::ImageFileReader<VectorImageType> ReaderType;

	/**and then read the noise */



	// let's read the image
	ReaderType::Pointer rs=ReaderType::New();
	rs->SetFileName(SFN);
	rs->Update();

	// let's read the image
	ReaderType::Pointer rn=ReaderType::New();
	rn->SetFileName(NFN);
	rn->Update();

	std::cout << typeid(rn->GetOutput()).name() << '\n';

	using VConstIteratorType = itk::ImageRegionConstIterator< VectorImageType >;

	 VConstIteratorType inputIt(rn->GetOutput(), rn->GetOutput()->GetLargestPossibleRegion());

	  	  		      inputIt.GoToBegin();


	  	  		      while( !inputIt.IsAtEnd() )
	  	  		        {
	  	  		   std::cout<<inputIt.Get()  <<std::endl;
	  	  		        ++inputIt;
	  	  		        }



	// 16x3 matrix, elements not initialize


	typedef cm::VectorImageReconstructor<VectorImageType,ScalarImageType> ReconstructorFilter;
	ReconstructorFilter::Pointer filter = ReconstructorFilter::New();
	filter->SetInput(rs->GetOutput());
	filter->SetNoiseKSpace(rn->GetOutput());
	filter->SetKSpaceDimension(cm::BIDIMENSIONAL);


//	std::cout<<filter->CalculateNoiseBW()<<std::endl;

//	std::cout<<filter->GetNoiseCoefficientMatrix()<<std::endl;
//
//	std::cout<<filter->GetNoiseCovarianceMatrix()<<std::endl;

//	if(vm.count("outputimage"))
//	{
//
//
//
//		typedef itk::ImageFileWriter<ScalarImageType >  W3Type;
//		typename W3Type::Pointer w3 = W3Type::New();
//		w3->SetInput(filter->GetOutput());
//		w3->SetFileName(OFN);
//		w3->Update();
//
//
//	}


	return EXIT_SUCCESS;
}
