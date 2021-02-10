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

#include "../../Utils/cmReconstructor.h"
#include "../../Utils/cmReconstructorRootSumOfSquares.h"

#include "itkImageRegionConstIterator.h"

#include "../../Utils/cmISMRMRDToITKImageFilter.h"


#include "../../Utils/cmClasses.h"

#include <iostream>   // std::cout
#include <string>     // std::string, std::stoi

#include "../../Utils/cm.h"
#include "../../Utils/mropt.h"

int main(int argc, char *argv[]) {

	/** Signal file name. */
	std::string SFN;
	/** Noise file name. */
	std::string NFN;

	/** Noise file name. */
	std::string OFN;

	/**Noise typee. */
	std::string NType="noiseFile";

	/**Reconstruciton type. */
	std::string RECON="RSS";

	/**Noise Bandwidth evaluation method. */
	std::string NBWType="readfromnoise";
	bool NBW=true;
	float NBWv= std::nanf("1");

	bool UCM=false;


	/** Body coil. */
	std::string BC="no";


	/**noise typee. */
	std::string SENS="inner";


	// acceleration
	int frequencyAcceleration=1;
	int phase0Acceleration=1;
	int phase1Acceleration=std::nanf("1");//this part will be usefull for the 3D acceleration

	//autocalibrations
	int frequencyAutocalibrationLines=1;
	int phase0AutocalibrationLines=1;
	int phase1AutocalibrationLines=std::nanf("1"); //this part will be usefull for the 3D acceleration




	/** The options */
	po::options_description desc("The software calculates 2D SNR from two ismrmrd file, noise and signal: Allowed options");
	desc.add_options()("help,h", "Produce HELP message")
																	("signalfilename,s",po::value<std::string>(&SFN), " Signal Filename (H5)")
																	("noiselfilename,n",po::value<std::string>(&NFN), " Noise Filename (H5)")
																	("noisetype,t",po::value<std::string>(&NType), " Noise Type (noiseFile,selfMulti,selfSingle)")
																	("outputimage,o",po::value<std::string>(&OFN), " OUTPUTIMAGE file")
																	("noisebandwidthvalue,b",po::value<float>(&NBWv), " noise BW value")
																	("noisebandwidthmethod,t",po::value<std::string>(&NBWType), " How do you want your bandwidth to be evaluated? (readfromnoise,calculatefromnoise,setvalue,no)")
//																	("usecovariancematrix,U",po::value<bool>(&UCM), " use noisecovariance matrix (false)")
																	("reconstructiontype,r",po::value<std::string>(&RECON), " Reconstruction type (RSS)")
																	("coilsensitivitymethod,S", po::value<std::string>(&SENS), "coilsensitivitymethod (INNER)")
																	//  														    ("bc,B", po::value<std::string>(&BC), " the BC id s set the sensemethod is set to BC")
																	("frequencyacceleration,A",po::value<int>(&frequencyAcceleration), " frequency acceleration (1)")
																	("phaseacceleration,B",po::value<int>(&phase0Acceleration), " phase acceleration (1) ");
																	("frequencyautocalibrationlines,L",po::value<int>(&frequencyAutocalibrationLines), " frequency acceleration (1)");
	   																("phaseautocalibration,M",po::value<int>(&phase0AutocalibrationLines), " phase acceleration (0) ");

	;


	po::options_description myerror("test main class SNR array combining\n Eros Montin, PhD");

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






	// let's read the image
	mroptimum::H5ReaderType ::Pointer rs=mroptimum::H5ReaderType::New();
	rs->SetIsmrmrdFileName(SFN);
	rs->SetaverageImageThroughAverages(true);
	rs->SetremoveFrequencyOversampling(true);
	rs->Update();


	mroptimum::H5ReaderType::Pointer rn=mroptimum::H5ReaderType::New();

	if (boost::iequals(NType,"noiseFile"))
	{
		rn->SetIsmrmrdFileName(NFN);
		rn->SetremoveFrequencyOversampling(false);
		rn->SetaverageImageThroughAverages(false);
		rn->Update();
	}

	//Filter Area
	mroptimum::ReconstructorType::Pointer filter=mroptimum::selectReconstructor(RECON);


	filter->SetInput(rs->GetOutput());
	filter->SetKSpaceDimension(rs->GetKSpaceDimension());
	filter->SetNoiseKSpace(rn->GetOutput());
	filter->SetNoiseBandWidthCorrection(true);


	//	Noise Bandwidth area
	if (boost::iequals(NBWType,"readfromnoise")){
		filter->SetNoiseBandWidth(rn->getRelativeReceiverBandwidth());
	}

	if (boost::iequals(NBWType,"calculatefromnoise")){
		//nothing to do already done in the class
	}
	if (boost::iequals(NBWType,"setvalue"))
	{filter->SetNoiseBandWidth(NBWv);}

	if (boost::iequals(NBWType,"no")){
		filter->SetNoiseBandWidthCorrection(false);
		filter->SetNoiseBandWidth(1.00);
	}


	//Sensitivity area

	if(filter->GetHasSensitivity())
	{
		if (boost::iequals(SENS,"INNER"))
		{

			filter->SetSensitivityMapCalculationMode(cm::SensitivityMapCalculation::INNER);


		}
		//		if(boost::iequals(SENS,"BC") || vm.count("bc"))
		//		{
		//
		//
		//			filter->SetSensitivityMapCalculationMode(cm::SensitivityMapCalculation::BODYCOIL);
		//
		//
		//			if (vm.count("bc")) {
		//
		//
		//				filter->SetBodyCoilSource(rb->GetOutput());
		//
		//			}
		//
		//		}
	}


	//Aceleration
	if(filter->GetHasAcceleration())
	{

		cm::VectorImageType::IndexType ACL={frequencyAutocalibrationLines, phase0AutocalibrationLines, phase1AutocalibrationLines};
		cm::VectorImageType::IndexType ACC={frequencyAutocalibrationLines, phase0AutocalibrationLines, phase1AutocalibrationLines};

		filter->SetAcceleration(ACC);
		filter->SetAutocalibrationLines(ACL);
	}


	filter->Update();

	if(vm.count("outputimage"))
	{



		typedef itk::ImageFileWriter<cm::ScalarImageType >  W3Type;
		typename W3Type::Pointer w3 = W3Type::New();
		w3->SetInput(filter->GetOutput());
		w3->SetFileName(OFN);
		w3->Update();


	}


	return EXIT_SUCCESS;
}
