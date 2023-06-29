classdef cm2DGFactorSENSE<cm2DReconSENSE
    %this is an initial class derived by the code developped at cbi.
    %sensitivity.
    %http://mriquestions.com/senseasset.html
    %last update 2020 Dec
    %update 2020 Jan
    %---------------------------------------------------------------------
    %kspace is at the fully sampled size and with zeros on the accelerated
    %positions
    
    properties
        
    end
    
    methods
        
            
        
        
    
   
        
        
        
       
        
        function G=getOutput(this)
            
        
            %Riccardo Lattanzi riccardo.lattanzi@nyulangone.org
            %SENSE_img_recon December 2020
            [SS]=this.getSignalSize();
             nf=SS(1);
             np=SS(2);
             nc = this.getSignalNCoils();
             
          pw_signalrawdata=this.getPrewhitenedSignal();
            pw_sensmap=this.getCoilSensitivityMatrix();

            invRn=this.getInverseNoiseCovariancePrewhithened();
             R1=this.getAccelerationFrequency();
             R2=this.getAccelerationPhase();
            
             Rn = eye(nc); %it's prewhitened

            Rtot = R1*R2;
            
             k_temp=this.getUndersampledKSpaceZeroPadded(pw_signalrawdata,R1,R2);
            aliased_image = this.get2DKSIFFT(k_temp)*sqrt(Rtot);
            
            
          imfold = this.shrinktoaliasedmatrix_2d(aliased_image,R1,R2);
            
            G=zeros(nf,np);
            
            
            for irow=1:floor(nf/R1)
                for icol=1:floor(np/R2)
                    current_R1 = length(irow:floor(nf/R1):nf);
                    current_R2 = length(icol:floor(np/R2):np);
                    current_Rtot = current_R1*current_R2;
                    s = squeeze(pw_sensmap(irow:floor(nf/R1):nf,icol:floor(np/R2):np,:));   % Gather the aliased pixels into the sensitivity matrix
                    s = reshape(s,[current_Rtot nc]);
                    s = s.';
%                     u = pinv(s'*invRn*s)*(s'*invRn);
                    u(isnan(u)) = 0;
            G(irow:floor(nf/R1):nf,icol:floor(np/R2):np) = reshape(sqrt(abs(diag(pinv(s'*invRn*s)).*diag(s'*invRn*s))),[current_R1 current_R2]);
                end
            end
            
                G=this.demimicAccelerated2DImage(G);

        end
        
        
    end
 
end


