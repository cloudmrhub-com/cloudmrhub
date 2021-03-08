
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



#include "vnl/vnl_matrix.h"
#include "../../Utils/cmSNRMultipleReplicasFilter.h"



//#include "cm.h";

int main(int argc, char *argv[]) {

	/** the Kspace data. */
	std::string IM;
	/** Second Image. */
	std::string OFN;
	/** Second Image. */
	std::string ITERATEON="repetitions";

	/** Second Image. */
	std::string RECON="RSS";


	/** Noise file name. */
	std::string NFN;

	/**noise typee. */
	std::string NType="noiseFile";

	bool NBW=false;
	float NBWv= std::nanf("1");

	bool UCM=false;

	long int numberofreplicas=2;

	std::string SENS="inner";

	std::string BC="no";

	/** The options */
	po::options_description desc("The software calculates 3D SNR from two images (reconstructed!!) ");
	desc.add_options()("help,h", "Produce HELP message")
								("KSPACE replicas IMAGE,s",po::value<std::string>(&IM), "Kspace Replicas Filename")
								("noiselfilename,n",po::value<std::string>(&NFN), " Noise Filename (H5)")
								("ITERATEON,i", po::value<std::string>(&ITERATEON), " Iterate on (repetitions,averages)")
								("RECON,r", po::value<std::string>(&RECON), " Iterate on (RSS,B1)")
								("outputimage,o",po::value<std::string>(&OFN), " OUTPUTIMAGE file")
								("noisetype,t",po::value<std::string>(&NType), " Noise Type (noiseFile,selfMulti,selfSingle)")
								("numberofreplica,R", po::value<long int>(&numberofreplicas), " number of replicas to use")

								;



	po::options_description myerror("compute the SNR of two reconstructed images on a roi");

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


	typedef cm::ISMRMRDToITKImageFilter<VectorImageType> RF;

	vnl_matrix<VectorImageType::InternalPixelType> A;

	if (boost::iequals(RECON,"RSS"))
	{

	}else{

		if (boost::iequals(NType,"noiseFile"))
		{
			RF::Pointer rn=RF::New();
			rn->SetIsmrmrdFileName(NFN);
			rn->Update();
			A=vectorImageToVNLMatrix<VectorImageType>(rn->GetOutput());


		}
	}
	/**read the replicas */

	RF::Pointer rs=RF::New();
	rs->SetIsmrmrdFileName(IM);
	rs->SetaverageImageThroughAverages(true);
	rs->SetremoveFrequencyOversampling(true);
	//	rs->Update();

	std::vector<VectorImageType::Pointer> K=rs->getReplicas(ITERATEON);

	std::cout<<ITERATEON<<std::endl;



	if( vm.count("numberofreplica"))
		{

		}else{
				numberofreplicas=K.size();
		}
	std::cout<< "number of replicas: "<<numberofreplicas<<std::endl;



	//instantiate the MR filter
	typedef cm::MultipleReplicasFilter<ScalarImageType> SNRFilterType;
	SNRFilterType::Pointer MR=SNRFilterType::New();

	typedef cm::ArrayCombiningMethodRSS<VectorImageType,ScalarImageType> ReconstructionFilterTypeRSS;
	typedef cm::ArrayCombiningMethodB1<VectorImageType,ScalarImageType> ReconstructionFilterTypeB1;


	RF::Pointer rb=RF::New();
	if( vm.count("bc"))
	{
		rb->SetIsmrmrdFileName(BC);
		rb->SetaverageImageThroughAverages(true);
		rb->SetremoveFrequencyOversampling(true);
		rb->Update();
	}

	omp_set_num_threads(omp_get_num_procs()*2);
#pragma omp parallel for shared(MR)
	for (int t=0; t<numberofreplicas;t++){


		ScalarImageType::Pointer IM =ScalarImageType::New();

		if (boost::iequals(RECON,"RSS"))
		{

			ReconstructionFilterTypeRSS::Pointer reconfilter = ReconstructionFilterTypeRSS::New();
			reconfilter->SetInput(K.at(t));
			reconfilter->SetKSpaceDimension(rs->GetKSpaceDimension());
			reconfilter->SetKSpaceDimension(rs->GetKSpaceDimension());
			reconfilter->Update();
			//images.push_back(reconfilter->GetOutput());


			//	writeImage<ScalarImageType>(reconfilter->GetOutput(),"it_"+ std::to_string(t)+".mha" );


			IM=reconfilter->GetOutput();

		}else if (boost::iequals(RECON,"B1")){

			ReconstructionFilterTypeB1::Pointer filter = ReconstructionFilterTypeB1::New();
			filter->SetInput(K.at(t));
			filter->SetKSpaceDimension(rs->GetKSpaceDimension());
			filter->SetNoise(A);
			if (std::isnan(NBWv)){
				filter->SetNoiseBandWidthCorrection(NBW);
			}else{
				filter->SetNoiseBandWidthCorrection(true);
				filter->SetNoiseBandWidth(NBWv);
			}
			filter->SetSensitivityMapSource(K.at(t)); //inner



			if (boost::iequals(SENS,"INNER"))
			{

				filter->SetSensitivityMapCalculationMode(cm::SensitivityMapCalculation::INNER);


			}
			if(boost::iequals(SENS,"BC") || vm.count("bc"))
			{


				filter->SetSensitivityMapCalculationMode(cm::SensitivityMapCalculation::BODYCOIL);


				if (vm.count("bc")) {



					filter->SetBodyCoilSource(rb->GetOutput());

				}

			}




			filter->Update();
			IM=filter->GetOutput();
		}


		//PUSH the image on the stack
		if(t==0){
			MR->SetInput(IM);
		}

		MR->pushReconstructedImage(IM);
	}




	writeImage<ScalarImageType>(MR->GetOutput(),OFN);





	return EXIT_SUCCESS;
}
