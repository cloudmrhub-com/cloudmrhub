#pragma once
#include "cmISMRMRDToITKImageFilter.h"
#include "cmReconstructor.h"
#include "cmReconstructorRootSumOfSquares.h"
#include "cmReconstructorB1Weighted.h"
#include "cmReconstructormSENSE.h"

#ifndef CM_H
#define CM_H
#include "cm.h"
#endif 

#include <boost/algorithm/string.hpp>

namespace mroptimum
{
	using H5ReaderType = cm::ISMRMRDToITKImageFilter<cm::VectorImageType>;
	using ReconstructorRSSType =cm::ReconstructorRootSumOfSquares<cm::VectorImageType,cm::ScalarImageType>;
	using ReconstructorType =cm::Reconstructor<cm::VectorImageType,cm::ScalarImageType>;
	using ReconstructorB1Type =cm::ReconstructorB1Weighted<cm::VectorImageType,cm::ScalarImageType>;
	using ReconstructormSENSEType =cm::ReconstructormSENSE<cm::VectorImageType,cm::ScalarImageType>;

ReconstructorType::Pointer selectReconstructor(std::string IN){

	const std::string str = boost::algorithm::to_lower_copy(IN);


	ReconstructorType::Pointer filter;

	if (boost::iequals(str,"rss")){

		ReconstructorRSSType::Pointer recon= ReconstructorRSSType::New();

		filter= recon;
		}

//	if (boost::iequals(IN,"ORIGIN")){
//
//		ReconstructorType::Pointer recon= ReconstructorType::New();
//
//		filter= recon;
//		}

	if (boost::iequals(str,"b1")){

		ReconstructorB1Type::Pointer recon= ReconstructorB1Type::New();

		filter= recon;
		}

	if (boost::iequals(str,"sense")){

			ReconstructormSENSEType::Pointer recon= ReconstructormSENSEType::New();

			filter= recon;
			}


if (boost::iequals(str,"msense")){

		ReconstructormSENSEType::Pointer recon= ReconstructormSENSEType::New();

		filter= recon;
		}

return filter;


}

}
