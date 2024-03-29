classdef cm2DReconAdapt<cm2DRecon
    
   
    
    methods

        function this=cm2DReconAdapt()
           this.setHasAcceleration(0);
            this.setHasSensitivity(0);       
        end
        
        function [im,sens_m]=getOutput(this)
             rawdata=this.getPrewhitenedSignal();
            [SS]=this.getSignalSize();
             nf=SS(1);
             np=SS(2);
             nc = this.getSignalNCoils();
             img_matrix=this.get2DKSIFFT(rawdata);
            Rn=eye(nc);
           [im,sens_m] = adapt_array_2d(img_matrix,Rn,1);
        end
        
        
        function sens_m=getOutputSensitivityMap(this)
            [~,sens_m]=this.getOutput(this);
        end
        
        
                   function TEST=test(this)
         TEST=this.therecon();
         TEST.plotImageAfterTest(TEST.getOutput(),'recon');
         TEST.whatHappened();
               end
            
   
    end    
        
 
        methods (Static,Access=protected)
            function TEST=therecon()
                 TEST=cm2DReconAdapt();
                [K,nc]=TEST.getMRoptimumTestData();
                TEST.setSignalKSpace(K);
                TEST.setNoiseCovariance(nc);
                               
            end
    
        end
end




