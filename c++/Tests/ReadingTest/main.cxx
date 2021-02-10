#include "../../Utils/cmRootSumOfSquareImageFilter.h"
#include <omp.h>

#include <boost/make_shared.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/locale/encoding_utf.hpp>
using boost::locale::conv::utf_to_utf;


#include "itkVectorImage.h"

#include "itkComplexToModulusImageFilter.h"

#include "../../Utils/cmISMRMRDToITKImageFilter.h"

#include "../../Utils/cmiFFTPhasedArrayFilter.h"



#include "itkExtractImageFilter.h"

#include "itkTileImageFilter.h"

#include "itkImageFileWriter.h"

#include "../../Utils/cm.h"

#include "../../Utils/ImageUtils.h"

int main(int argc, char *argv[]) {

	/** Signal file name. */

	std::string IN;
	/** Output file name. */
	std::string ON;
	/** Output file name. */
	std::string FON;


	/** */
	bool AVG =false;
	bool ROS = false;

	/** The options */
	po::options_description desc("Allowed options");
	desc.add_options()("help,h", "Produce HELP message")
																									("IN,s",po::value<std::string>(&IN), "Signal File Name H5")
																									("FON,f", po::value<std::string>(&FON), " Output File Name not sliced")
																									("AVG,a", po::value<bool>(&AVG), " Averaged Kspace")
																									("ROS,r", po::value<bool>(&ROS), " Remove OS")
																									;



	po::options_description myerror("Read an ISMRMRD file and calculate the Root Sum of Squares");

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
	using PixelType = vcl_complex< PT >;
	const unsigned int      	Dimension = 3;


	//					using PixelType = float;

	typedef itk::VectorImage<PixelType, Dimension >  VectorImageType;
	typedef itk::Image<PixelType, Dimension >    ScalarImageType;
	typedef itk::Image<PT, Dimension >    ResultsImageType;


	typedef itk::VectorImage<PixelType, Dimension-1 >  VectorSliceImageType;
	typedef itk::Image<PixelType, Dimension-1 >    ScalarSliceImageType;


	typedef cm::ISMRMRDToITKImageFilter<VectorImageType> RF;
	RF::Pointer r=RF::New();
	r->SetIsmrmrdFileName(IN);
	r->SetremoveFrequencyOversampling(ROS);
	r->SetaverageImageThroughAverages(AVG);

	r->Update();








//			typedef cm::iFFTPhasedArrayFilter<VectorImageType> TEST;
//			TEST::Pointer coil= TEST::New();

//			coil->SetInput(r->GetOutput());
//			coil->SetKSpaceDimension(r->GetKSpaceDimension());
//			coil->Update();

//			writeImage<VectorImageType>(coil->GetOutput(),"coil.mha");




	if(vm.count("FON"))
	{











		typedef cm::RootSumOfSquareImageFilter<VectorImageType,ScalarImageType> P2;
		P2::Pointer pp2= P2::New();

		pp2->SetInput(r->GetOutput());
		pp2->SetKSpaceDimension(r->GetKSpaceDimension());
		pp2->Update();



		typedef itk::ComplexToModulusImageFilter<ScalarImageType,ResultsImageType> MFilterType;
		MFilterType::Pointer mfilter = MFilterType::New();
		mfilter->SetInput(pp2->GetOutput());
		mfilter->Update();


		std::cout<<r->getNumberOfRepetition()<<std::endl;
			std::cout<<r->getNumberOfCoils()<<std::endl;
			std::cout<<r->getRelativeReceiverBandwidth()<<std::endl;
			std::cout<<r->getFieldStrength()<<std::endl;
			std::cout<<r->getNumberOfSlice()<<std::endl;
			std::cout<<r->getTrajectory()<<std::endl;
			std::cout<<r->isAccelerated()<<std::endl;

			std::cout<<r->getAcceleration1()<<std::endl;
			std::cout<<r->getAcceleration2()<<std::endl;
			std::cout<<r->getWaveForm()<<std::endl;
			std::cout<<r->getKSpaceDimensionString()<<std::endl;
			std::cout<<r->getEncodedSize()<<std::endl;

			std::cout<<r->getEncodedResolution()<<std::endl;
			std::cout<<r->getReconResolution()<<std::endl;

		typedef itk::ImageFileWriter<ResultsImageType>  W3Type;
		typename W3Type::Pointer w3 = W3Type::New();
		w3->SetInput(mfilter->GetOutput());
		w3->SetFileName(FON);
		w3->Update();


	}






	return EXIT_SUCCESS;
}
