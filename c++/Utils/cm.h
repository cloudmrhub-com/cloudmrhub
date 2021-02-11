#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageRegionIterator.h"
//

namespace cm
{
#pragma once
	enum KSpaceAcquisitionDimension {
		BIDIMENSIONAL = 2,
		TRIDIMENSIONAL = 3
	};


#pragma once
	enum SensitivityMapCalculation {
		INNER = 1,
		BODYCOIL = 2
	};


	using ComplexType = float;

	using PixelType = std::complex<ComplexType>;

	using VectorImageType = itk::VectorImage<PixelType, 3>;
	using VectorImageInternalPixelType = VectorImageType::InternalPixelType;
	using VectorImagePixelType = VectorImageType::PixelType;
	using VectorImageTypePointer = VectorImageType::Pointer;
	using VectorImagePixelElelmentType = itk::VariableLengthVector<PixelType>;
	using VectorImageRegionIteratorType = itk::ImageRegionIterator<VectorImageType>;


	using ScalarImageType = itk::Image<PixelType, 3>;
	using ScalarImageTypePointer = ScalarImageType::Pointer;
	using ScalarImagePixelType = ScalarImageType::PixelType;
	using ScalarImageTypePointer = ScalarImageType::Pointer;
	using ScalarImageRegionIteratorType = itk::ImageRegionIterator<ScalarImageType>;

	using RealScalarImageType = itk::Image<ComplexType, 3>;

	using  ChannelArrayType=  vnl_matrix<VectorImageInternalPixelType>;




}

//https://itk.org/Doxygen/html/SphinxExamples_2src_2Core_2Common_2FilterAndParallelizeImageRegion_2Code_8cxx-example.html#_a10

