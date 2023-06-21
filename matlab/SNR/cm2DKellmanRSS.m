classdef cm2DKellmanRSS<cm2DReconRSS
    properties
    end
    methods
               function SNR=getOutput(this)
               [SS]=this.getSignalSize();
               nf=SS(1);
               np=SS(2);
                nc = this.getSignalNCoils();
                
            pw_signalrawdata=this.getPrewhitenedSignal();
            SNR = zeros(nf,np);
                img_matrix = this.get2DKSIFFT(pw_signalrawdata); %MRifft and th normalization
                 for irow = 1:nf
                     for icol = 1:np
                         SNR(irow,icol) = abs(sqrt(2*(squeeze(img_matrix(irow,icol,:))')*squeeze(img_matrix(irow,icol,:))));
                     end
                 end
            
        end
    end
end




