classdef cm2DReconSENSE<cm2DReconWithSensitivityAutocalibrated
   
    
    
    methods
        
       function this=cm2DReconSENSE()

            this.setAutocalibrationFrequency(NaN);
        end        
        
        
            function o=getAutocalibration(this)
            o=this.getAutocalibrationPhase();
        end
        
        function setAutocalibration(this,acl)
            this.setAutocalibrationPhase(acl);
        end
        
   
        
         
        
        
        
        
        
        
       
        
        function MRimage=getOutput(this)
            
        
            %Riccardo Lattanzi riccardo.lattanzi@nyulangone.org
            %SENSE_img_recon December 2020
            
              [SS]=this.getSignalSize();
             nf=SS(1);
             np=SS(2);
             nc = this.getSignalNCoils();
             R1=this.getAccelerationFrequency();
             R2=this.getAccelerationPhase();
            ACL=this.getAutocalibration();
            pw_signalrawdata=this.getPrewhitenedSignal();
            k_temp=this.getUndersampledKSpaceZeroPadded(pw_signalrawdata,R1,R2);
            img_matrix=this.get2DKSIFFT(k_temp);
            pw_sensmap=this.getCoilSensitivityMatrix();

            invRn=this.getInverseNoiseCovariancePrewhithened();
             

            Rtot = R1*R2;
            
            
            img_matrix = img_matrix*sqrt(Rtot);
            
            
            
            imfold=   this.shrinktoaliasedmatrix_2d(img_matrix,R1,R2);
            
            
            
            MRimage = zeros(nf,np);
            
           % G=figure();
            %m=zeros(nf,np);
            for irow=1:floor(nf/R1)
                for icol=1:floor(np/R2)
                    current_R1 = length(irow:floor(nf/R1):nf);
                    current_R2 = length(icol:floor(np/R2):np);
                    current_Rtot = current_R1*current_R2;
                    s = squeeze(pw_sensmap(irow:floor(nf/R1):nf,icol:floor(np/R2):np,:));   % Gather the aliased pixels into the sensitivity matrix
                    s = reshape(s,[current_Rtot nc]);
                    s = s.';
                    u = pinv(s'*invRn*s)*(s'*invRn);
                    u(isnan(u)) = 0;
                    MRimage(irow:floor(nf/R1):nf,icol:floor(np/R2):np)=reshape(u*squeeze(imfold(irow,icol,:)),[current_R1 current_R2]); % put this recon into the right places in the image
                  %  m(irow:floor(nf/R1):nf,icol:floor(np/R2):np)=m(irow:floor(nf/R1):nf,icol:floor(np/R2):np)+1;
                   % imagesc(m);colorbar(); pause(0.1);
                end
            end
            
            MRimage=this.demimicAccelerated2DImage(MRimage);
        end
        
        

        
    end
      methods (Static)
           
            function TEST=therecon()
                %instantiate the class
                TEST=cm2DReconSENSE();
                [K,nc]=TEST.getMRoptimumTestData();
                
                R1=2;R2=4;
                A1=20;
                A2=20;
                
                KS=undersampleSense2D(K,R1,R2);
                TEST.setSignalKSpace(KS);
                TEST.setNoiseCovariance(nc);
                %please set noise before set the sens
                TEST.setCoilSensitivityMatrixSourcePrewhitened(KS)
                TEST.setAccelerationFrequency(R1)
                TEST.setAccelerationPhase(R2)
                 TEST.setAutocalibrationFrequency(A1)
                TEST.setAutocalibrationPhase(A2)
                TEST.setCoilSensitivityMatrixCalculationMethod('simplesense')
                            end
        end
end


