#ifndef __cmPhasedArraySensitivityMapFilterBodyCoil_hxx
#define __cmPhasedArraySensitivityMapFilterBodyCoil_hxx

#include "cmPhasedArraySensitivityMapFilterBodyCoil.h"
#include "itkImageAlgorithm.h"
//#include "NYi.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "cmiFFTPhasedArrayFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "vnl/vnl_vector.h"

namespace cm
{

template< class VectorImageType>
void PhasedArraySensitivityMapFilterBodyCoil<VectorImageType>::ThreadedGenerateData(const VectorImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
{



	VectorImageTypePointer COILIFFT=this->GetIfft();
	VectorImageTypePointer OUT=this->GetOutput();
	VectorImageTypePointer BC=this->GetBodyCoil();



	//ifftBC
	VectorImageTypePointer BCIFFT;
	BCIFFT=GetBCIfft();


	itk::ImageRegionIterator<VectorImageType> icoil(COILIFFT,outputRegionForThread);

	itk::ImageRegionIterator<VectorImageType> ibc(BCIFFT,outputRegionForThread);
	itk::ImageRegionIterator<VectorImageType> io(OUT,outputRegionForThread);

	//				/* prepare */
	icoil.GoToBegin();

	ibc.GoToBegin();
	io.GoToBegin();
	//
	//
	int nChan=COILIFFT->GetNumberOfComponentsPerPixel();
	/** the number of channel */
	int nChanBC=BCIFFT->GetNumberOfComponentsPerPixel();

	/* the coil array values of the pixel*/
	VectorImagePixelType coil;
	coil.SetSize(nChan);
	coil.Fill(0);


	/* the output pixel type*/
	VectorImagePixelType bc;
	bc.SetSize(nChanBC);
	bc.Fill(0);



	/* the output pixel type*/
	VectorImagePixelType divisor;
	divisor.SetSize(nChan);
	divisor.Fill(0);

	VectorImageInternalPixelType repeatingdivisor;
	VectorImageInternalPixelType ji =VectorImageInternalPixelType(1,1);
	VectorImageInternalPixelType N=std::sqrt(VectorImageInternalPixelType(2.0,0));

	/* the output*/
	VectorImagePixelType out;
	out.SetSize(COILIFFT->GetNumberOfComponentsPerPixel() );







	///* actual calc*/
	while( !icoil.IsAtEnd() )
	{
		//get the coil
		coil=icoil.Get();
		//get the BC
		bc=ibc.Get();

		switch(nChanBC){
		break;
		case 2:
			repeatingdivisor=(bc.GetElement(1)+ ( ji*bc.GetElement(2)))/N;
			break;
		default:
			repeatingdivisor=bc.GetElement(1);

		}





		for (int t=0;t<coil.GetSize();t++)
		{

			out.SetElement(t,(VectorImageInternalPixelType)(coil.GetElement(t)/repeatingdivisor));

		};







		io.Set(out);
		++io;
		++icoil;
		++ibc;
	}













}
template< class VectorImageType>
void PhasedArraySensitivityMapFilterBodyCoil< VectorImageType>
::AfterThreadedGenerateData(){


	//	VectorImageTypePointer OUT =this->GetOutput();
	//
	//	vnl_matrix<typename VectorImageType::InternalPixelType>O;
	//	O=vectorImageToVNLMatrix<VectorImageType>(OUT);
	//
	//
	//
	//	VectorImagePixelType MAX;
	//	MAX.SetSize(O.rows());
	//	//MAX.fill(0);
	//	VectorImagePixelType tmp;
	//	tmp.SetSize(O.rows());
	//	VectorImagePixelType tmp2;
	//	tmp2.SetSize(O.rows());
	//
	//
	//
	//	vnl_vector<VectorImageInternalPixelType> VO(O.cols());
	//	//there's no way to order complex number so i'm comparing the real part (like matlab)
	//	vnl_vector<float> VOR(O.cols());
	//	unsigned m;
	//
	//	typename VectorImageType::IndexType Ii;
	//
	//	for(auto t=0;t<O.rows();t++){
	//		VO=O.get_row(t);
	//		for(auto p=0;p<VO.size();p++){
	//
	//			VOR(p)=std::abs(VO(p));
	//
	//		}
	//
	//		m=VOR.arg_max();
	//
	//		MAX.SetElement(t,VO(m));
	//
	//}
	//
	//
	//
	//itk::ImageRegionIteratorWithIndex<VectorImageType> io(OUT,OUT->GetLargestPossibleRegion());
	//
	//io.GoToBegin();
	/////* actual calc*/
	//while( !io.IsAtEnd() )
	//{
	//	//get the coil
	//	tmp=io.Get();
	//
	//	try{
	//		Ii= io.GetIndex();
	//	}catch (const std::exception& e)
	//	{
	//		std::cerr << e.what();
	//	}
	//
	//	for (int t=0;t<tmp.GetSize();t++)
	//	{
	//		tmp2=tmp.GetElement(t)/MAX.GetElement(t);
	//
	//
	//	}
	//	io.Set(tmp2);
	//	io++;
	//}

	this->rescaleOutputVectorImage();

}

template< class VectorImageType>
void PhasedArraySensitivityMapFilterBodyCoil< VectorImageType>
::BeforeThreadedGenerateData(){
	//we calculate the sos and the ifft:)

	VectorImageTypePointer BCC=GetBCIfft();
	std::cout<<"start 2D ifft BC\n";
	typedef cm::iFFTPhasedArrayFilter<VectorImageType> IFFTFilter0;
	typename IFFTFilter0::Pointer ifFilter0= IFFTFilter0::New();//=this->GetIfft();
	ifFilter0->SetKSpaceDimension(this->GetKSpaceDimension());
	ifFilter0->SetInput(this->GetBodyCoil());
	ifFilter0->Update();
	BCC=ifFilter0->GetOutput();
	this->SetBCIfft(BCC);



	VectorImageTypePointer IF=this->GetIfft();

	typedef cm::iFFTPhasedArrayFilter<VectorImageType> IFFTFilter;
	typename IFFTFilter::Pointer ifFilter= IFFTFilter::New();//=this->GetIfft();

	ifFilter->SetKSpaceDimension(this->GetKSpaceDimension());





	ifFilter->SetInput(this->GetInput());
	ifFilter->Update();

	IF=ifFilter->GetOutput();
	this->SetIfft(IF);




}




}// end namespace


#endif
