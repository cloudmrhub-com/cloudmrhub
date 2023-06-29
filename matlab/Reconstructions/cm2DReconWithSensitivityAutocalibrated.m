classdef cm2DReconWithSensitivityAutocalibrated<cm2DReconWithSensitivity
    %CoilSenseitivityMap is the actual Sensitivitymap used to calculate the SNR
    %RealCoilSensitivityMap is the MROPTDATA with the sampled sensitivity
    %method
    %SensitivityCalculationMethod, a string {'simplesense','adaptive','zerofilling','bodycoil'}
    %To get SNR we neeed to fet S (sensitivity) and Noisecovariance matrix
    %v24122020
    %the main differences with the parents class is that the kspace of the
    %Source is accelerated with some ACL
    properties(Access=protected)
        AutocalibrationF
        AutocalibrationP
        AccelerationF
        AccelerationP
        MimicAcceleratedAcquisition
        PreMimicAcceleratedSize
    end
    
    methods%(Access=protected)
        %this function should be overrided on othe
        %implementations e.g. when you are using accelerated
        %methods
        
        function this=cm2DReconWithSensitivityAutocalibrated()
            this.setHasAcceleration(1);
            this.setHasSensitivity(1);
        end
        
        function coilsens_set=calculateCoilSensitivityMatrix(this)
            
            coilsens_set=this.calculateCoilSensitivityMatrixUnderSampled();
            if isempty(coilsens_set)
                coilsens_set=this.calculateCoilSensitivityMatrixFullySampled();
            end
            
        end
        
                function coilsens_set=getCoilSensitivityMatrixForFullySampledKSpaceSimpleSense(this)
            %in this case sens_matrix is the full image
            nc=this.getSignalNCoils;
            sensmap_temp = MRifft(this.getCoilSensitivityMatrixSource(),[1,2]);
            ref_img = sqrt(sum(abs(sensmap_temp).^2,3));
            coilsens_set = sensmap_temp./repmat(ref_img,[1 1 nc]);
            if(this.getMaskCoilSensitivityMatrix())
            sensmask = ref_img>mean(ref_img(:));
            sensmaskrep = repmat(sensmask,[1 1 nc]);
            coilsens_set = coilsens_set.*sensmaskrep;
            end
            %prewhitening
            coilsens_set = this.prewhiteningSignal(coilsens_set, this.getNoiseCovariance() );
            this.logIt(['sensitivity map calculated as simplesense'],'ok');
            
                end
        
        function coilsens_set=getCoilSensitivityMatrixForUnderSampledKSpaceSimpleSense(this)
            %in this case sens_matrix is the full image
            nc=this.getSignalNCoils;
            [SS]=this.getSignalSize();
            nf=SS(1);
            np=SS(2);
            ACLp=this.getAutocalibrationPhase();
            
            ACLf=this.getAutocalibrationFrequency();
            
            TOT=ACLf*ACLp; %only theoretical!
            if isempty(ACLf)  %sense
                ACLf=NaN;
                TOT=nf*ACLp;
            else
                if(isnan(ACLf)) %sense
                    ACLf=NaN;
                    TOT=nf*ACLp;
                end
                
            end
            sens_data=this.getAutocalibrationsLinesKSpaceZeroPadded(this.getSignalKSpace(), ...
                ACLf,ACLp);
            
            sensmap_temp = MRifft(sens_data,[1,2])*sqrt(TOT);
            ref_img = sqrt(sum(abs(sensmap_temp).^2,3));
            coilsens_set = sensmap_temp./repmat(ref_img,[1 1 nc]);
            
            
                        if(this.getMaskCoilSensitivityMatrix())
                
            sensmask = ref_img>mean(ref_img(:));
            sensmaskrep = repmat(sensmask,[1 1 nc]);
            coilsens_set = coilsens_set.*sensmaskrep;
            end

            coilsens_set = this.prewhiteningSignal(coilsens_set, this.getNoiseCovariance() );
            
            
            
            this.logIt(['sensitivity map calculated as simplesense'],'ok');
            
        end
        
        
        function coilsens_set=calculateCoilSensitivityMatrixUnderSampled(this)
            coilsens_set=[];
            nchan=this.getSignalKSpace();
            switch(lower(this.CoilSensitivityMatrixCalculationMethod))
                case {'bartsense'}
                    %   https://mrirecon.github.io/bart/examples.html#2
                    SOURCE=this.bartMy2DKSpace(this.getCoilSensitivityMatrixPrewhitened());
                    try
                        c=['caldir ' num2str(this.getAutocalibrationPhase())];
                        sens = bart(c, SOURCE); %sense
                        % sens = bart('slice 4 0', calib);
                        coilsens_set = this.debartMy2DKSpace(sens);
                    catch
                        coilsens_set = espirit_sensitivitymap(this.getCoilSensitivityMatrixPrewhitened(),this.getAutocalibration());
                        
                    end
                case {'espirit','espiritaccelerated'}
                    SOURCE=this.bartMy2DKSpace(this.getCoilSensitivityMatrixPrewhitened());
                    %[calib ~] = bart(['ecalib  -r ' num2str(this.getAutocalibration()) ' -m 2'], SOURCE); %sense
                    [calib, ~] = bart(['ecalib  -r ' num2str(this.getAutocalibrationPhase()) ], SOURCE); %sense
                    sens = bart('slice 4 0', calib);
                    coilsens_set = this.debartMy2DKSpace(sens);
                    
                case {'espiritv1'}
                    %in this case sens_matrix is the full image
                    this.logIT(['sensitivity map calculated as espirit'],'ok');
                    
                    coilsens_set = espirit_sensitivitymap(this.getCoilSensitivityMatrixPrewhitened(),this.getAutocalibrationPhase());
                case {'simplesenseacl','internal referenceacl','inneracl'}
                    this.logIt(['sensitivity map calculated as simplesense'],'?');
                    coilsens_set=this.getCoilSensitivityMatrixForUnderSampledKSpaceSimpleSense();
            end
        end
    end
    
    methods
        
        
        function o=getAutocalibrationFrequency(this)
            o=this.AutocalibrationF;
        end
        
        
        function setAutocalibrationFrequency(this,acl)
            this.AutocalibrationF=acl;
        end
        
        
        function o=getAutocalibrationPhase(this)
            o=this.AutocalibrationP;
        end
        
        
        function setAutocalibrationPhase(this,acl)
            this.AutocalibrationP=acl;
        end
        
        
        
        function o=getAccelerationFrequency(this)
            o=this.AccelerationF;
        end
        
        function setAccelerationFrequency(this,R1)
            this.AccelerationF=R1;
        end
        
        
        function o=getAccelerationPhase(this)
            o=this.AccelerationP;
        end
        
        function setAccelerationPhase(this,R2)
            this.AccelerationP=R2;
        end
        
        
        
        
        function o=getMimicAcceleratedAcquisition(this)
            o=this.MimicAcceleratedAcquisition;
        end
        
        
        function setMimicAcceleratedAcquisition(this,m)
            this.MimicAcceleratedAcquisition=m;
        end
        
  

       
        function K=mimicmSenseDataFromFullysampledZeroPadded(this,K,accf,accp,ACL)
            %get the Kspace zeropadded
            %it's a static method
            this.setMimicAcceleratedAcquisition(true);
            this.PreMimicAcceleratedSize=size(K);
            this.logIt('KSpace Mimicked accelerated','warning')
            if(this.needsRegridding(K,accf,accp))            
            [K]=undersample2DDataExtreme(K,accf,accp,NaN,ACL);
    this.logIt(['The matrix size interpolated to be divisible by the acceleration factor ' num2str(accf) ' x '  num2str(accf)] ,'warning')
            else
                [K]=undersample2DData(K,accf,accp,NaN,ACL);
            end
        end                
            
        
        
        function K=mimicmGrappaDataFromFullysampledZeroPadded(this,K,accf,accp,aclf,aclp)
            %get the Kspace zeropadded
            %it's a static method
           
            this.setMimicAcceleratedAcquisition(true);
            this.PreMimicAcceleratedSize=size(K);
               this.logIt('KSpace Mimicked accelerated','warning')
             if(this.needsRegridding(K,accf,accp))            
            [K]=undersample2DDataExtreme(K,accf,accp,aclf,aclp);
            this.logIt(['The matrix size interpolated to be divisible by the acceleration factor ' num2str(accf) ' x '  num2str(accf)] ,'warning')
            else
                 [K]=undersample2DData(K,accf,accp,aclf,aclp);
             end
            
             
        end
        
        function K=demimicAccelerated2DImage(this,K)
            
              if(this.needsBackRegridding(K))
                  P=this.PreMimicAcceleratedSize;
                            K=this.resize2DMatrixbyInterpolationVoxelSpace(K,[1:P(1)],[1:P(2)]);
             end
        end
        
         function o=needsBackRegridding(this,K)
         o=false;
         S=size(K);
         %was it mimicked?
         if(this.getMimicAcceleratedAcquisition)
             %if the size now is different by the one stored in PreMimic
             if(sum(S(1:2)~=this.PreMimicAcceleratedSize(1:2))>0)
                 this.logIt('demimicked ','warning')
                 o=true;
             end
         end
         end
        
        
        
    end
    
    methods(Static)
        function o=needsRegridding(K,accf,accp)
            o=false;        
            S=size(K);
            if (sum(mod(S(1:2),[accf accp]))~=0)
                o=true;
            end
        end
            
        function ACS=getAutocalibrationsLinesKSpace(K,ACLf,ACLp)
            %ACLf=naN to get all the lines o the freq directions
            nf=size(K,1);
            np=size(K,2);
            
            if isnan(ACLf)
                ACS = K(:,floor(np/2)-floor(ACLp/2)+1:floor(np/2)+floor(ACLp/2),:);
            else
                ACS = K(floor(nf/2)-floor(ACLf/2)+1:floor(nf/2)+floor(ACLf/2),floor(np/2)-floor(ACLf/2)+1:floor(np/2)+floor(ACLf/2),:);
            end
            
        end
        
        function ACS=getAutocalibrationsLinesKSpaceZeroPadded(K,ACLf,ACLp)
            %ACLf=naN to get all the lines o the freq directions
            nf=size(K,1);
            np=size(K,2);
            
            ACS=zeros(size(K));
            
            if isempty(ACLf)  %sense
                Y=floor(np/2)-floor(ACLp/2)+1:floor(np/2)+floor(ACLp/2);
                ACS(:,Y,:) = K(:,Y,:);
            else
                
                if isnan(ACLf) %sense
                    Y=floor(np/2)-floor(ACLp/2)+1:floor(np/2)+floor(ACLp/2);
                    ACS(:,Y,:) = K(:,Y,:);
                else
                    
                    
                    Y=floor(np/2)-floor(ACLp/2)+1:floor(np/2)+floor(ACLp/2);
                    X=floor(nf/2)-floor(ACLf/2)+1:floor(nf/2)+floor(ACLf/2);
                    
                    ACS(X,Y,:) = K(X,Y,:);
                    
                    
                    
                end
                
            end
        end
        
        
        function K=getUndersampledKSpace(K,accf,accp)
            %get the shrinked kspaceKspace
            %it's a static method
            [~,K]=undersample2DData(K,accf,accp,0,0);
        end
        
        
        
        
        function K=getUndersampledKSpaceZeroPadded(K,accf,accp)
            %get the Kspace zeropadded
            %it's a static method
            [K]=undersample2DData(K,accf,accp,0,0);
        end
        
        
        
        
        
        %
        %                                     function K=getUndersampledKSpaceExtreme(K,accf,accp)
        %                 %get the shrinked kspaceKspace
        %                 %it's a static method
        %                 [~,K]=undersample2DDataExtreme(K,accf,accp,0,0);
        %             end
        %
        %
        %
        %
        %             function K=getUndersampledKSpaceZeroPaddedExtreme(K,accf,accp)
        %                 %get the Kspace zeropadded
        %                 %it's a static method
        %                 this.logIt(['The matrix size interpolated to be divisible by the acceleration factor ' num2str(accf) ' x '  num2str(accf)] ,'warning')
        %                 [K]=undersample2DDataExtreme(K,accf,accp,0,0);
        %             end
        %
        %             function K=mimicmSenseDataFromFullysampledZeroPaddedExtreme(K,accf,accp,ACL)
        %                 %get the Kspace zeropadded
        %                 %it's a static method
        %                 [K]=undersample2DDataExtreme(K,accf,accp,NaN,ACL);
        %
        %             end
        %
        %                         function K=mimicmGrappaDataFromFullysampledZeroPaddedExtreme(K,accf,accp,aclf,aclp)
        %                 %get the Kspace zeropadded
        %                 %it's a static method
        %                 [K]=undersample2DDataExtreme(K,accf,accp,aclf,aclp);
        %                         end
        %
        
        
        
        
        
    end
end

