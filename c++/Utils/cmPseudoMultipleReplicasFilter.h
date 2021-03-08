//https://itk.org/Wiki/ITK/Examples/Developer/ImageFilter



#ifndef __cmPseudoMultipleReplicasFilter_h
#define __cmPseudoMultipleReplicasFilter_h

#include "itkImageToImageFilter.h"
#include "itkNumericTraits.h"

#include "cmMultipleReplicasFilter.h"

#include "itkImageRegionIteratorWithIndex.h"

#include "itkVectorImage.h"

#include "vnl/algo/vnl_complex_eigensystem.h"
namespace cm
{

/** \class PseudoMultipleReplicasFilter
 * \brief This class calculates SNR for Pseudo Multiple replicas.
 * replicas are stored in a std::vectorobj
 *
 * \author Eros Montin,PhD \email eros.montin@nyulangone.org
 * \author Prof. Riccardo Lattanzi,PhD \email riccardo.lattanzi@nyulangone.org
 *
 * \note This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536, and it was performed under the rubric of the Center for Advanced Imaging Innovation and Research (CAI2R, www.cai2r.net), a NIBIB Biomedical Technology Resource Center (P41 EB017183). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
 *
 */


template< class ScalarImageType>
class PseudoMultipleReplicasFilter:public MultipleReplicasFilter< ScalarImageType >
{
public:
	/* Standard class typedefs. (similar to using..)*/
	typedef PseudoMultipleReplicasFilter                                                 Self;
	typedef MultipleReplicasFilter< ScalarImageType> Superclass;
	typedef itk::SmartPointer< Self >                                Pointer;
	typedef itk::SmartPointer< const Self >                                ConstPointer;

	/* Method for creation through the object factory. */
	itkNewMacro(Self);

	/* Run-time type information (and related methods). */
	itkTypeMacro(MultipleReplicasFilter, MultipleReplicasFilter);




	typedef typename ScalarImageType::Pointer ScalarImageTypePointer;
	typedef typename ScalarImageType::PixelType ScalarImagePixelType;
	typedef typename ScalarImageType::InternalPixelType ScalarImageInternalPixelType;
	typedef typename ScalarImageType::IndexType ScalarImageIndexType;

	typedef typename std::vector<ScalarImageTypePointer> 	VectorofScalarImageType;
	typedef typename std::vector<ScalarImageInternalPixelType> 	VectorofScalarInternalImagePixelType;
   typedef typename itk::ImageRegionIteratorWithIndex<ScalarImageType>        ScalarImageIteratorType;
	typedef typename ScalarImageType::RegionType          ScalarImageRegionType;



		typedef itk::VectorImage <ScalarImagePixelType, 3 >  VectorImageType;
		typedef typename VectorImageType::Pointer VectorImageTypePointer;
		typedef typename VectorImageType::PixelType VectorImagePixelType;
		typedef typename VectorImageType::InternalPixelType VectorInternalPixelType;
		typedef typename VectorImageType::IndexType VectorImageIndexType;


		void GetN();

		itkSetMacro(CorrelationNoiseFactor,vnl_matrix<ScalarImageType::InternalPixelType>);
		itkGetMacro(CorrelationNoiseFactor,vnl_matrix<ScalarImageType::InternalPixelType>);
		//itkSetMacro(NoiseCovarianceMatrix,vnl_matrix<ScalarImageType::InternalPixelType>);

		itkGetMacro(NoiseCovarianceMatrix,vnl_matrix<ScalarImageType::InternalPixelType>);
		void SetNoiseCovarianceMatrix(vnl_matrix<ScalarImageType::InternalPixelType> NC){
			//store the value in the class
			this->m_NoiseCovarianceMatrix=NC;

			//calculate the noisecorrelationonisefactor to be used in the pseudonoisegenerator

			//eigenvalues for complex is templated on double so we need to cast out covariance matrix
			vnl_matrix<vcl_complex<double>>doubletmp(NC.rows(),NC.cols());
			doubletmp=cast2DVNLMatrixfromto<ScalarImagePixelType, vcl_complex<double> >(NC);

			//CAluclate Eiganvalues
			vnl_complex_eigensystem Eigs(doubletmp,true,false);

				vnl_matrix<vcl_complex<double>>D(NC.rows(),NC.cols());
				D.fill(0.0);
				//from the formula we collect the eigenvalues and we make the sqrt
				for (auto t=0;t<Eigs.N;t++){
								D(t,t)=std::sqrt(Eigs.W(t));
							}

//				D.apply(sqrt);

				vnl_matrix<vcl_complex<double>>iR=vnl_matrix_inverse<vcl_complex<double>>(Eigs.R.transpose());

				vnl_matrix<vcl_complex<double>>TMP=Eigs.R.transpose()*D*iR;

				vnl_matrix<ScalarImagePixelType>TMP2=cast2DVNLMatrixfromto< vcl_complex<double>,ScalarImagePixelType >(TMP);
				this->SetCorrelationNoiseFactor(TMP2);

		};


		VectorImageTypePointer getPseudoReplicaKspace(VectorImageTypePointer S);


protected:
	PseudoMultipleReplicasFilter(){}
	~PseudoMultipleReplicasFilter(){}


	  /* Does the real work. */
//	  virtual void GenerateData();



private:
	PseudoMultipleReplicasFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented

	vnl_matrix<ScalarImageType::InternalPixelType> m_NoiseCovarianceMatrix;

	//set when u set the noisecovariance matrix
	vnl_matrix<ScalarImageType::InternalPixelType> m_CorrelationNoiseFactor;








};
} //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "cmPseudoMultipleReplicasFilter.hxx"
#endif


#endif // __cmPseudoMultipleReplicasFilter_h
