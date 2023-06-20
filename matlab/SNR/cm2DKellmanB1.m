classdef cm2DKellmanB1<cm2DReconB1
    properties
    end
    methods
            function SNR=getOutput(this)
            %this method works on 2d images only

            pw_signalrawdata=this.getPrewhitenedSignal();
            img_matrix=this.get2DKSIFFT(pw_signalrawdata);
            pw_sensmap=this.getCoilSensitivityMatrix();
            
            %we are using the inverse of eye because the signal and the sens matrix are already
            %prewhitened
            invRn=this.getInverseNoiseCovariancePrewhithened();
            
             [SS]=this.getSignalSize();
             nf=SS(1);
             np=SS(2);
             nc = this.getSignalNCoils();
             
            for irow = 1:nf
                for icol = 1:np
                    s_matrix = squeeze(pw_sensmap(irow,icol,:));
                    SNR(irow,icol) = sqrt(2)*abs((s_matrix')*invRn*squeeze(img_matrix(irow,icol,:)))/...
                        sqrt((s_matrix')*invRn*s_matrix);
                end
            end
            
        end

    end
    
    
    
  
    
    
    
end






