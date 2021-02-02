#include "../../utils/NYi.h"

#include <boost/algorithm/string.hpp>

#include <boost/make_shared.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/locale/encoding_utf.hpp>
using boost::locale::conv::utf_to_utf;


#include "itkVectorImage.h"

#include <iostream>
#include <fstream>

#include "../../Module/cmISMRMRDToITKImageFilter.h"
int main(int argc, char *argv[]) {

	/** Signal file name. */
	std::string SFN;
	/** Output Image file name. */
	std::string ON;
	/** Remove Oversampling True by default. */
	int OV=1;
	/** Average kSpace through averages. */
	int A=1;



	/** The options */
	po::options_description desc(" Allowed options");
	desc.add_options()("help,h", "Produce HELP message")
															("input,i",po::value<std::string>(&SFN), "Input K-Space raw data file (ISMRMRD v1)")
															("output,o", po::value<std::string>(&ON), " Output file name")
															("oversampling,v", po::value<int>(&OV), " Remove Oversampling (1 (dflt) for true, 0 for no) ")
															("averages,a", po::value<int>(&A), " Average K Space though Averages (1 (dflt) for true, 0 for no) ")
															;


	po::options_description myerror("Transforms ISMRMRD v1 file to image:\n\n\n\n Eros Montin, PhD");

	po::variables_map vm;

	if (argc < 2) {
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
	const unsigned int      	Dimension = 3;


	//					using PixelType = float;

	typedef itk::VectorImage<PixelType, Dimension >  VectorImageType;
	typedef itk::Image<PixelType, Dimension >    ScalarImageType;
	typedef itk::Image<PT, Dimension >    ResultsImageType;


	typedef cm::ISMRMRDToITKImageFilter<VectorImageType> RF;
	RF::Pointer rs=RF::New();

	/**read the image */
	rs->SetIsmrmrdFileName(SFN);
	rs->SetremoveFrequencyOversampling(true);
	rs->SetaverageImageThroughAverages(true);


	if(OV==0) //remove Oversampling
	{
		rs->SetremoveFrequencyOversampling(false);
	}

	if(A==0) // average over averages
		{
		rs->SetaverageImageThroughAverages(false);
		}

	rs->Update();
		VectorImageType::Pointer s=rs->GetOutput();


		if(vm.count("output"))
		{

			writeImage<VectorImageType>(s,ON);

		}




		return EXIT_SUCCESS;
}
