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
             R1=this.getAccelerationFrequency();
             R2=this.getAccelerationPhase();
            ACL=this.getAutocalibration();
          pw_signalrawdata=this.getPrewhitenedSignal();
%         k_temp              = zeros(size(pw_signalrawdata));
%         k_temp(1:R1:nf,1:R2:np,:)  = pw_signalrawdata(1:R1:nf,1:R2:np,:);
            k_temp=undersamplemSense2D(pw_signalrawdata,R1,R2,ACL);
            img_matrix=this.get2DKSIFFT(k_temp);
            pw_sensmap=this.getCoilSensitivityMatrixPrewhitened();

            invRn=this.getInverseNoiseCovariancePrewhithened();
             

            Rtot = R1*R2;
            
            
            % Since only nf*np/Rtot samples are actually sampled, we need to scale the image by the square root of the
            % total acceleration factor to maintain the same scaling as for the noise (fully sampled)
            img_matrix = img_matrix*sqrt(Rtot);
            
            
            
            imfold=   this.shrinktoaliasedmatrix_2d(img_matrix,R1,R2);
            
            
            
            G = zeros(nf,np);
            
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
            G(irow:floor(nf/R1):nf,icol:floor(np/R2):np) = reshape(sqrt(abs(diag(pinv(s'*invRn*s)).*diag(s'*invRn*s))),[current_R1 current_R2]);
                end
            end
        end
        
    G=this.demimicAccelerated2DImage(this,G);
        
    end
 
end


