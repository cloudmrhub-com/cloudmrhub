#include "../../lNYi/cm2DSNRArrayCombiningMethodRss.h"

#include "../../lNYi/NYi.h"

#include <boost/algorithm/string.hpp>

#include <boost/make_shared.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/locale/encoding_utf.hpp>
using boost::locale::conv::utf_to_utf;


#include "itkVectorImage.h"
#include "lNYi/json.hpp"

#include <iostream>
#include <fstream>

#include "../../lNYi/cm3DISMRMRDToITKImageFilter.h"
#include "vnl/vnl_matrix.h"
#include "../../lNYi/cmSumOfSquareImageFilter.h"
#include "itkExtractImageFilter.h"
int main(int argc, char *argv[]) {

	/** Signal file name. */
	std::string SFN;
	/** Noise file name. */
	std::string NFN;

	/** FA map file name. */
	std::string FAN;
	/** Output Image file name. */
	std::string ON;
	/**Json Output file name. */
	std::string JO;
	/**noise typee. */
	std::string NType="noiseFile";


	/** */
	/** Use covariance matrix set to false. */
	int UC=0;
	bool UCM=false;
	float NB=0;




	/** The options */
	po::options_description desc("The software calculates 2D SNR from two ismrmrd file, noes and signal: Allowed options");
	desc.add_options()("help,h", "Produce HELP message")
													("signalfilename,s",po::value<std::string>(&SFN), "Signal Filename")
													("noiselfilename,n",po::value<std::string>(&NFN), " Noise Filename")
													("flipanglefilename,f",po::value<std::string>(&FAN), " Flip Angle Map Filename")
													("noisetype,t",po::value<std::string>(&NType), " Noise Type (noiseFile,selfMulti,selfSingle)")
													("outputmap,o", po::value<std::string>(&ON), " Output Name")
													("usecovariancematrix,c", po::value<int>(&UC), " Use covariance matrix default false")
													("nbwidth,b", po::value<float>(&NB), " Noise BW")
													("JO,j", po::value<std::string>(&JO), " Json Out")
													;


	po::options_description myerror("compute the SNR as rss (-u 0) or rss psi (-u 1) psi\n\n\n\n Eros Montin, PhD");

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



	std::cout<< "UC"<<UC<<"is UC equal to one:"<< UC+1 <<" and UCM"<<UCM<<std::endl;
	if (UC==1){
		UCM=true;
	}

	std::cout<< "UC"<<UC<<" and UCM"<<UCM<<std::endl;



	using PT = float;
	using PixelType = std::complex< PT >;
	const unsigned int      	Dimension = 3;


	//					using PixelType = float;

	typedef itk::VectorImage<PixelType, Dimension >  VectorImageType;
	typedef itk::Image<PixelType, Dimension >    ScalarImageType;
	typedef itk::Image<PT, Dimension >    ResultsImageType;

	typedef itk::Image<PT, Dimension-1 >    ResultsSliceImageType;
	typedef itk::VectorImage<PixelType, Dimension-1 >  VectorSliceImageType;
	typedef itk::Image<PixelType, Dimension-1 >    ScalarSliceImageType;



	typedef itk::cm3DISMRMRDToITKImageFilter<VectorImageType> RF;
	RF::Pointer rs=RF::New();

	/**read the image */
	rs->SetIsmrmrdFileName(SFN);
	rs->Update();
	VectorImageType::Pointer s=rs->GetOutput();


	/**and then read the noise */


	if (boost::iequals(NType,"noiseFile"))
	{
	    // Strings are identical
	}
	RF::Pointer rn=RF::New();
	rn->SetIsmrmrdFileName(NFN);
	rn->Update();
	VectorImageType::Pointer n=rn->GetOutput();





	//		const VectorImageType * inputImage = s->GetOutput();



	VectorImageType::RegionType inputRegion = s->GetLargestPossibleRegion();
	VectorImageType::SizeType size = inputRegion.GetSize();

	VectorImageType::RegionType nRegion = n->GetLargestPossibleRegion();
	VectorImageType::SizeType nsize = nRegion.GetSize();


	int sl= size[2];
	int nsl= nsize[2];
	size[2] = 0; // we extract along z direction
	nsize[2] = 0; // we extract along z direction
	VectorImageType::IndexType start = inputRegion.GetIndex();
	VectorImageType::IndexType nstart = nRegion.GetIndex();


	nlohmann::json jo;

	jo["version"]="20180831";
	jo["author"]="eros.montin@gmail.com";
	jo["type"]="ACM";

	nlohmann::json images=nlohmann::json::array();
	nlohmann::json SNRslice=nlohmann::json::array();
	nlohmann::json COVslice=nlohmann::json::array();

	for (auto t=0; t<sl;t++)
	{
		using ExtractFilterType = itk::ExtractImageFilter< VectorImageType, VectorSliceImageType >;
		ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
		extractFilter->SetDirectionCollapseToSubmatrix();



		// r->GetOutput()->Print(std::cout,0);

		// set up the extraction region [one slice]

		const unsigned int sliceNumber =t;
		start[2] = sliceNumber;
		VectorImageType::RegionType desiredRegion;
		desiredRegion.SetSize(  size  );
		desiredRegion.SetIndex( start );


		extractFilter->SetExtractionRegion( desiredRegion );
		extractFilter->SetDirectionCollapseToSubmatrix();
		extractFilter->InPlaceOn();
		extractFilter->SetInput( s );
		extractFilter->Update();




		ExtractFilterType::Pointer extractFiltern = ExtractFilterType::New();
		extractFiltern->SetDirectionCollapseToSubmatrix();





		if (t<nsl)
		{
			nstart[2] = t;
		}else{
			nstart[2] = 0;
		}

		VectorImageType::RegionType ndesiredRegion;
		ndesiredRegion.SetSize(  nsize  );
		ndesiredRegion.SetIndex( nstart );


		extractFiltern->SetExtractionRegion( ndesiredRegion );
		extractFiltern->SetDirectionCollapseToSubmatrix();
		extractFiltern->InPlaceOn();
		extractFiltern->SetInput( n );
		extractFiltern->Update();






		typedef itk::cm2DSNRArrayCombiningMethodRss<VectorSliceImageType,ScalarSliceImageType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(extractFilter->GetOutput());
		filter->SetNoise(extractFiltern->GetOutput());
		filter->SetUseCovarianceMatrix(UCM);
		filter->Update();




		if(vm.count("outputmap"))
		{


			typedef itk::Image<PT, Dimension-1 > MAGImageType;
			typedef itk::ComplexToModulusImageFilter<ScalarSliceImageType,MAGImageType> MFilterType;
			MFilterType::Pointer mfilter = MFilterType::New();
			mfilter->SetInput(filter->GetOutput());
			mfilter->Update();


			writeImage<MAGImageType>(mfilter->GetOutput(),std::to_string(t) + ON);

		}


		if(vm.count("JO"))
		{

			if (t==0)
			{



				vnl_matrix<PixelType>L;
				L=filter->GetCovarianceMatrix();
				//L.print(o);
				COVslice.push_back(sliceWriterjson<PT>(L));
			}



			/** this class can be thought as an array*/
			SNRslice.push_back(sliceWriterjson<ScalarSliceImageType,ResultsSliceImageType>(filter->GetOutput()));
		}





	}


	if(vm.count("JO"))
	{



		images.push_back({{"slice",SNRslice},{"imageName","SNR"}});
		images.push_back({{"slice",COVslice},{"imageName","Noise Covariance"}});
		//images.push_back({{"slice",SOSslice},{"imageName","SOS Image"}});

		jo["images"]=images;
		std::ofstream o(JO);
		o << jo << std::endl;
	}








	return EXIT_SUCCESS;
}
