#ifndef __cmPhasedArraySensitivityMapFilter_hxx
#define __cmPhasedArraySensitivityMapFilter_hxx

#include "cmPhasedArraySensitivityMapFilter.h"
#include "itkImageAlgorithm.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "utils.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

namespace cm
{



template< class VectorImageType>
void PhasedArraySensitivityMapFilter< VectorImageType>
::rescaleOutputVectorImage(){


	VectorImageTypePointer OUT =this->GetOutput();

	vnl_matrix<typename VectorImageType::InternalPixelType>O;
	O=vectorImageToVNLMatrix<VectorImageType>(OUT);



	VectorImagePixelType MAX;
	MAX.SetSize(O.rows());
	//MAX.fill(0);
	VectorImagePixelType tmp;
	tmp.SetSize(O.rows());
	VectorImagePixelType tmp2;
	tmp2.SetSize(O.rows());



	vnl_vector<VectorImageInternalPixelType> VO(O.cols());
	//there's no way to order complex number so i'm comparing the abs part (like matlab)
	vnl_vector<float> VOR(O.cols());
	unsigned m;

	typename VectorImageType::IndexType Ii;

	for(auto t=0;t<O.rows();t++){
		VO=O.get_row(t);
		for(auto p=0;p<VO.size();p++){
				//compute the abs
			VOR(p)=std::abs(VO(p));

		}
		// find the maximum for the channel
		m=VOR.arg_max();
		//set the element maximum into the array that later we will use to normalize the coil sens
		MAX.SetElement(t,VO(m));

	}



	itk::ImageRegionIteratorWithIndex<VectorImageType> io(OUT,OUT->GetLargestPossibleRegion());

	io.GoToBegin();
	// actual calc
	while( !io.IsAtEnd() )
	{
		//get the coil array
		tmp=io.Get();

		try{
			Ii= io.GetIndex();
		}catch (const std::exception& e)
		{
			std::cerr << e.what();
		}
	//divide by the max of that coil
		for (int t=0;t<tmp.GetSize();t++)
		{
			tmp2=tmp.GetElement(t)/MAX.GetElement(t);


		}
		io.Set(tmp2);
		io++;
	}



}





}// end namespace


#endif
