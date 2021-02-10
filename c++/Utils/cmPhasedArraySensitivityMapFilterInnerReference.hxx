#ifndef __cmPhasedArraySensitivityMapFilterInnerReference_hxx
#define __cmPhasedArraySensitivityMapFilterInnerReference_hxx

#include "cmPhasedArraySensitivityMapFilterInnerReference.h"
#include "itkImageAlgorithm.h"


#include "cmRootSumOfSquareImageFilter.h"

#include "cmiFFTPhasedArrayFilter.h"


#include "itkVectorImage.h"
namespace cm
{

template< class VectorImageType>
void PhasedArraySensitivityMapFilterInnerReference< VectorImageType>
::ThreadedGenerateData(const VectorImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
{
//before the threaded we calculated the SOS and the coil sensitivity now we are building the image as COIL/ rep(SOS)

	VectorImageTypePointer COILIFFT=this->GetIfft();
	ScalarImageTypePointer SOS=this->GetSos();
	VectorImageTypePointer OUT=this->GetOutput();

				itk::ImageRegionIterator<VectorImageType> it0(COILIFFT,outputRegionForThread);
				itk::ImageRegionIterator<ScalarImageType> it1(SOS,outputRegionForThread);
				itk::ImageRegionIterator<VectorImageType> ot(OUT,outputRegionForThread);

//				/* prepare */
				it0.GoToBegin();
				it1.GoToBegin();
				ot.GoToBegin();
//
//
				/* the coil array values of the pixel*/
				VectorImagePixelType coil;


				/* the output pixel type*/
				ScalarImagePixelType sos;

				/* the output*/
				VectorImagePixelType out;

				out.SetSize(COILIFFT->GetNumberOfComponentsPerPixel() );


				/* BANG!!*/
				while( !it0.IsAtEnd() )
				{
					coil=it0.Get();
					sos=it1.Get();

					/* reset counter and calculate the sos*/


					/*sum of squares*/
					for (auto t=0;t<coil.GetSize();t++)
					{

					out[t]=coil.GetElement(t)/sos;

					}


					/* root*/
					ot.Set(out);
					++it0;
					++it1;
					++ot;
				}







}

template< class VectorImageType>
void PhasedArraySensitivityMapFilterInnerReference< VectorImageType>
::BeforeThreadedGenerateData(){
//we calculate the sos and the ifft:)

	std::cout<<"start reconstruction SOS\n";
	typedef cm::RootSumOfSquareImageFilter<VectorImageType,ScalarImageType> P2;
			typename P2::Pointer pp2= P2::New();



			ScalarImageTypePointer SOS=this->GetSos();

			VectorImageTypePointer KK=this->GetInput();

			pp2->SetInput(KK);
			pp2->SetKSpaceDimension(this->GetKSpaceDimension());
			pp2->Update();
			SOS=pp2->GetOutput();



			this->SetSos(SOS);



			VectorImageTypePointer IF=this->GetIfft();

				typedef cm::iFFTPhasedArrayFilter<VectorImageType> IFFTFilter;
				typename IFFTFilter::Pointer ifFilter= IFFTFilter::New();//=this->GetIfft();
				
				ifFilter->SetKSpaceDimension(this->GetKSpaceDimension());





				ifFilter->SetInput(KK);
				ifFilter->Update();
					
					IF=ifFilter->GetOutput();
					this->SetIfft(IF);




}




}// end namespace


#endif
