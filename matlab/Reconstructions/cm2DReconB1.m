classdef cm2DReconB1<cm2DReconWithSensitivity
    
    
 
    
    methods
        
        
        
     function this=cm2DReconB1()
           this.setHasAcceleration(0);
            this.setHasSensitivity(1);       
        end
        
        
        function im=getOutput(this)
            %this method works on 2d images only
            pw_signalrawdata=this.getPrewhitenedSignal();
            img_matrix=this.get2DKSIFFT(pw_signalrawdata);
            pw_sensmap=this.getCoilSensitivityMatrix();

            invRn=this.getInverseNoiseCovariancePrewhithened();
             [SS]=this.getSignalSize();
             nf=SS(1);
             np=SS(2);
             nc = this.getSignalNCoils();
             im = zeros(nf,np);
            for irow=1:nf
                for icol=1:np
                    s_matrix = squeeze(pw_sensmap(irow,icol,:));
                    im(irow,icol) = (s_matrix')*invRn*squeeze(img_matrix(irow,icol,:));
                end
            end
        end
    
    
    
       
            
   
    end    
        
 
        methods (Static)
            function TEST=therecon()
                 TEST=cm2DReconB1();
                [K,nc]=TEST.getMRoptimumTestData();
                TEST.setSignalKSpace(K);
                TEST.setNoiseCovariance(nc);
                TEST.setCoilSensitivityMatrixCalculationMethod('simplesense')
                TEST.setCoilSensitivityMatrixSourcePrewhitened(K);
            end
        end
        
        
      
    
    
    
end






