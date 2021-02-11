//https://itk.org/Wiki/ITK/Examples/Developer/ImageFilter



#ifndef __cmPrewhiteningPhasedArray_h
#define __cmPrewhiteningPhasedArray_h

#include "itkImageToImageFilter.h"
#include "cmReconstructor.h"

namespace cm
{

/** \class PrewhiteningPhasedArray
 * \brief This class Prewhite 2D k-space phased array fully sampled raw data.
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
* \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 */


template< class TImage,class TOImage>
class PrewhiteningPhasedArray:public Reconstructor< TImage, TOImage >
{
public:
  /* Standard class typedefs. (similar to using..)*/
  typedef PrewhiteningPhasedArray                                                 Self;
  typedef itk::ImageToImageFilter< TImage, TOImage > Superclass;
  typedef itk::SmartPointer< Self >                                Pointer;
  typedef itk::SmartPointer< const Self >                                ConstPointer;

  /* Method for creation through the object factory. */
  itkNewMacro(Self);

  /* Run-time type information (and related methods). */
  itkTypeMacro(PrewhiteningPhasedArray, Reconstructor);



protected:
  PrewhiteningPhasedArray(){}
  ~PrewhiteningPhasedArray(){}

  /* Does the real work. */
  virtual void GenerateData();

private:
  PrewhiteningPhasedArray(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented




};
} //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmPrewhiteningPhasedArray.hxx"
#endif


#endif // __PrewhiteningPhasedArray_h
