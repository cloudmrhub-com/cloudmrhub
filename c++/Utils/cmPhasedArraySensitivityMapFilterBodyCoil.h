//https://itk.org/Wiki/ITK/Examples/Developer/ImageFilter



#ifndef __cmPhasedArraySensitivityMapFilterBodyCoil_h
#define __cmPhasedArraySensitivityMapFilterBodyCoil_h

#include "itkImageToImageFilter.h"
#include <typeinfo>
#include "itkImageRegionIterator.h"
#include "cm.h"
#include "cmPhasedArraySensitivityMapFilter.h"
namespace cm
{

/** \class PhasedArraySensitivityMapFilterBodyCoil
 * \brief This class calculates the sensitivity maps as described by Walsh et al. BodyCoil Reconstruction of Phased Array MR imagery \cite BodyCoilPhasedArrayReconstruction
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 *
 *
 */


template< class VectorImageType>
class PhasedArraySensitivityMapFilterBodyCoil:public PhasedArraySensitivityMapFilter<VectorImageType >
{
public:
  /* Standard class typedefs. (similar to using..)*/
  typedef PhasedArraySensitivityMapFilterBodyCoil                                                 Self;
  typedef PhasedArraySensitivityMapFilter< VectorImageType > Superclass;
  typedef itk::SmartPointer< Self >                                Pointer;
  typedef itk::SmartPointer< const Self >                                ConstPointer;

  /* Method for creation through the object factory. */
  itkNewMacro(Self);

  /* Run-time type information (and related methods). */
  itkTypeMacro(PhasedArraySensitivityMapFilterBodyCoil, PhasedArraySensitivityMapFilter);

/*  kljk the main image type*/
  
   typedef typename VectorImageType::Pointer VectorImageTypePointer;
   typedef typename VectorImageType::PixelType VectorImagePixelType;
   typedef typename VectorImageType::InternalPixelType VectorImageInternalPixelType;
   typedef typename VectorImageType::RegionType          VectorImageRegionType;
    typedef typename itk::ImageRegionIterator<VectorImageType>        VectorImageIteratorType;

	
	itkGetMacro(BodyCoil,VectorImageTypePointer);
	itkSetMacro(BodyCoil,VectorImageTypePointer);
	
		itkGetMacro(BCIfft,VectorImageTypePointer);
		itkSetMacro(BCIfft,VectorImageTypePointer);




protected:

	PhasedArraySensitivityMapFilterBodyCoil(){


	};
  ~PhasedArraySensitivityMapFilterBodyCoil(){};

  /* Does the real work. */
  //virtual void GenerateData();
  virtual void ThreadedGenerateData(const VectorImageRegionType& outputRegionForThread,
      itk::ThreadIdType threadId);
      
  //calculate the ifftt and the sos
  	void BeforeThreadedGenerateData();
    //calculate the ifftt and the sos
    	void AfterThreadedGenerateData();

private:
  PhasedArraySensitivityMapFilterBodyCoil(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

VectorImageTypePointer m_BodyCoil; 
  

	VectorImageTypePointer m_BCIfft;



};
} //namespace CM


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmPhasedArraySensitivityMapFilterBodyCoil.hxx"
#endif


#endif // __cmPhasedArraySensitivityMapFilterBodyCoil_h
