classdef cm2DReconSS<cm2DRecon
    %Root sum of squares class 
    %author:eros.montin@gmail.com
    %v012021
    %46&2 just ahead of me    
    methods
        function this=cm2DReconSS()
           this.setHasAcceleration(0);
            this.setHasSensitivity(0);       
        end
        function im=getOutput(this)
% reconstruct individual coils' images and apply FFT scale factor
% iFFT scales data by 1/N, so I need to scale back, as the noise covariance
% matrix is calculated in k-space (not Fourier transformed)
             rawdata=this.getPrewhitenedSignal();
             nc = this.getSignalNCoils();
             img_matrix=this.get2DKSIFFT(rawdata);
            im =sum(img_matrix.^2,3);
        end
    end
        methods (Static)
            function TEST=therecon() % edit this to get the test working
                 TEST=cm2DReconSS();
                [K,nc]=TEST.getMRoptimumTestData();
                TEST.setSignalKSpace(K);
                TEST.setNoiseCovariance(nc);
            end
        end
    
    
end




