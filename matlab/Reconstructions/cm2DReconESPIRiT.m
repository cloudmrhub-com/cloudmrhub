classdef cm2DReconESPIRiT<cm2DReconWithSensitivityAutocalibrated
   
    
    
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
        
   
        
         
        
        
        
        
        
        
       
        
        function o=getOutput(this)
            
        
            %Riccardo Lattanzi riccardo.lattanzi@nyulangone.org
            %SENSE_img_recon December 2020
            
              [SS]=this.getSignalSize();
             nf=SS(1);
             np=SS(2);
             nc = this.getSignalNCoils();
             R1=this.getAccelerationFrequency();
             R2=this.getAccelerationPhase();
            pw_signalrawdata=this.getPrewhitenedSignal();
            k_temp=this.getUndersampledKSpaceZeroPadded(pw_signalrawdata,R1,R2);
           pw_sensmap=this.getCoilSensitivityMatrix();
            sens=this.bartMy2DKSpace(pw_sensmap);
            K=this.bartMy2DKSpace(k_temp);
            o = bart('pics -r0.', K, sens);
            
       o=this.demimicAccelerated2DImage(o);
        
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


