#ifndef __cmPrewhiteningPhasedArray_hxx
#define __cmPrewhiteningPhasedArray_hxx

#include "cmPrewhiteningPhasedArray.h"
#include "itkImageAlgorithm.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cblas.h>

#include "vnl/vnl_matlab_write.h"

#include <typeinfo>



namespace cm
{




template< class TImage,class TOImage>
void PrewhiteningPhasedArray< TImage,TOImage>
::GenerateData()
 {

ChannelArrayType psi=this->GetNoiseCovarianceMatrix();


VectorImageTypePointer KSpace =this->GetInput();

typename TImage::SizeType size=KSpace->GetLargestPossibleRegion().GetSize();




	unsigned long long NP=size[0]*size[1];
	short int NC=KSpace->GetNumberOfComponentsPerPixel();



	//to get the chol lower we have to pass from VNL to gsl:)
		gsl_matrix_complex *A;
		A=gsl_matrix_complex_alloc(NC, NC);

		VectorImageInternalPixelType G;

		gsl_complex TEMP;

		for (auto x=0;x<NC;x++)
		{
			for (auto y=0;y<NC;y++)
			{
				G=psi(x,y);
				GSL_SET_COMPLEX(&TEMP, G.real(),G.imag());
				gsl_matrix_complex_set(A,x,y,TEMP);
			}
		}


		//	//https://www.gnu.org/software/gsl/manual/html_node/Cholesky-Decomposition.html
		//	//factorization is made on the matrix
			gsl_linalg_complex_cholesky_decomp(A);


			// my LU (L) implementation:) that get back the data to VNL
			gsl_complex z;
			//back to vnl
			for (auto x=0;x<NC;x++)
			{
				for (auto y=0;y<NC;y++)
				{
					if (y<=x)
					{
						z=gsl_matrix_complex_get(A,x,y);
						VectorImageInternalPixelType mycomplex(GSL_REAL(z),GSL_IMAG(z));
						psi(x,y)=mycomplex;
					}else{
						VectorImageInternalPixelType mycomplex(0,0);
						psi(x,y)=mycomplex;
					}
				}
			}



			//invert the correlation
			psi=vnl_matrix_inverse<VectorImageInternalPixelType>(psi);


// reshape te input as NCxnpixels

	ChannelArrayType Signal(NC,NP * size[2],0.0);

	ChannelArrayType sliceSignal;

	using SliceType = itk::VectorImage<PixelType, 2>;
	typename SliceType::Pointer axial;

#pragma omp parallel for shared(Signal)
	for (int slice=0; slice<size[2];slice++){
		axial=axiallySliceThisImage<TImage,SliceType>(KSpace,slice);
		sliceSignal=vectorImageToVNLMatrix<SliceType>(axial);
		Signal.update(sliceSignal,0,NP*slice);
	}




//	 finally!! prewhite the input
	Signal=psi*Signal;



	VectorImageType::IndexType Oi;

	VectorImageTypePointer Out=this->GetOutput();

	Out->SetBufferedRegion(Out->GetRequestedRegion());
	Out->Allocate();

	VectorImagePixelElelmentType F(NC);


	//reshapeit
	long int Arrindex=0;
#pragma parallel for shared(Out)
	for (auto sl=0;sl<size[2];sl++)
	{
		for (auto x=0;x<size[0]; x++)
		{

			for (auto y=0;y<size[1]; y++)
			{
				for (auto c=0;c<NC;c++)
				{

					F.SetElement(c, Signal(c,sl*NP + y*size[0]+x));

				}
				Oi={x,y,sl};
				Out->SetPixel(Oi,F);

			}
		}

 }
//	this->GetOutput()->Graft(Out);

 }





//
//
//
//
//	typename TImage::ConstPointer s = this->GetInput();
//	/** de reference to a non const pointer the input*/
//	InputImageTypePointer cs=ConstPointerToPointer<InputImageType>(s);
//
//
//
//
//	OutputImageTypePointer output = this->GetOutput();
//	output->SetBufferedRegion(output->GetRequestedRegion());
//	output->Allocate();
//
//
//	typename TImage::Pointer noise = this->GetNoise();
//
//
//
//	/* the number of channel */
//	int nChan=s->GetNumberOfComponentsPerPixel();
//
//
//
//
//
//	/* calculate noise correlation matrix */
//	std::vector<InputImageInnerPixelType> o;
//
//
//	typename InputImageType::SizeType sz;
//	sz=noise->GetLargestPossibleRegion().GetSize();
//
//	int NP= sz[0]*sz[1];
//
//	vnl_matrix<InputImageInnerPixelType>C(nChan,NP);
//	vnl_matrix<InputImageInnerPixelType>CT(NP,nChan);
//	vnl_matrix<InputImageInnerPixelType>S(nChan,NP);
//
//
//
//  	/* calculate signal correlation matrix */
//  		std::vector<InputImageInnerPixelType> oss;
//	for (int t;t<nChan; t++){
//		/* get channels */
//
//
//		o=VectorImageElementAsArray<InputImageInnerPixelType,TImage>(t,noise);
//
//		oss=VectorImageElementAsArray<InputImageInnerPixelType,TImage>(t,cs);
//
//		for(auto it = o.begin(); it != o.end(); ++it) {
//
//			int p=std::distance(o.begin(), it);
//
//			C(t,p)=*it;
//
//		}
//
//
//		for(auto it = oss.begin(); it != oss.end(); ++it) {
//
//			int p=std::distance(oss.begin(), it);
//
//			S(t,p)=*it;
//
//		}
//
//
//	};
//
//
//
//
//
//	/*	Noise Correlation Matrix Definition*/
//	vnl_matrix<InputImageInnerPixelType> RN(nChan,nChan);
//
//	CT=C.conjugate_transpose();
//	RN=C*CT;
//	RN=RN/(InputImageInnerPixelType)NP;
//
//
//
//	//to cet the chol lower we have to pass to gsl:)
//	gsl_matrix_complex *A;
//	A=gsl_matrix_complex_alloc(nChan, nChan);
//
//	InputImageInnerPixelType G;
//
//	gsl_complex TEMP;
//
//	for (auto x=0;x<nChan;x++)
//	{
//		for (auto y=0;y<nChan;y++)
//		{
//			G=RN(x,y);
//			GSL_SET_COMPLEX(&TEMP, G.real(),G.imag());
//
//			gsl_matrix_complex_set(A,x,y,TEMP);
//			std::cout<<abs(G)<<",";
//
//		}
//		std::cout<<"\n";
//	}
//
//
//
//
//	//https://www.gnu.org/software/gsl/manual/html_node/Cholesky-Decomposition.html
//	//factorization is made on the matrix
//	gsl_linalg_complex_cholesky_decomp(A);
//
//
//	// my LU (L) implementation:)
//	gsl_complex z;
//	//back to vnl
//	for (auto x=0;x<nChan;x++)
//	{
//		for (auto y=0;y<nChan;y++)
//		{
//
//
//			if (y<=x)
//			{
//				z=gsl_matrix_complex_get(A,x,y);
//				InputImageInnerPixelType mycomplex(GSL_REAL(z),GSL_IMAG(z));
//				RN(x,y)=mycomplex;
//			}else{
//				InputImageInnerPixelType mycomplex(0,0);
//				RN(x,y)=mycomplex;
//
//			}
//
//
//
//
//
//		}
//	}
//
//
//
//	//invert the correlation right before
//	RN=vnl_matrix_inverse<InputImageInnerPixelType>(RN);
//
//
//	// finally!! prewhit the input
//	S=RN*S;
//
//
//
//	typename InputImageType::SizeType szs;
//	szs=s->GetLargestPossibleRegion().GetSize();
//
//	typename OutputImageType::IndexType Oi;
//
//
//	InputImagePixelType F(nChan);
//
//	//reshapeit
//	long int Arrindex=0;
//
//		for (auto x=0;x<szs[0]; x++)
//		{
//
//			for (auto y=0;y<szs[1]; y++)
//			{
//
//
//				for (auto c=0;c<nChan;c++)
//				{
//
//					F.SetElement(c, S(c,y*szs[0]+x));
//
//				}
//
//				Oi={x,y};
//
//
//				output->SetPixel(Oi,F);
//
//
//			}
//		}
//
//
//
//
//
//
//
//
// }






}// end namespace


#endif
