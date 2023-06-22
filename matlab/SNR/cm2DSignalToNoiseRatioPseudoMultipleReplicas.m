classdef cm2DSignalToNoiseRatioPseudoMultipleReplicas<cm2DSignalToNoiseRatioMultipleReplicas
    %main class of array combining methods, the constructor is ovewritten by
    %the class constructor
    properties
        Reconstructor
        NumberOfPseudoReplicas
        ReferenceImage
    end
    methods
        function this = cm2DSignalToNoiseRatioPseudoMultipleReplicas(Recon,x)
            %the class expects a 3D matrix composed by a tile of 2D images
            %or nothing
            this.Type='PMR';
            if nargin>0
                this.setReconstructor(Recon)
            end
            
            if nargin>1
                this.add2DStackOfImages(x)
            end
            
        end
        function setReconstructor(this,R)
            this.Reconstructor=R;
        end
        function setNumberOfPseudoReplicas(this,f)
            
            this.NumberOfPseudoReplicas=f;
        end
        
        function o =getNumberOfPseudoReplicas(this)
            o = this.NumberOfPseudoReplicas;
        end
        function setNoiseKSpace(this,f)
            % 2DKspace
            this.Reconstructor.setNoiseKSpace(f);
        end
        function o =getNoiseKSpace(this)
            % 2DKspace
            o = this.Reconstructor.getNoiseKSpace();
        end
        function setSignalKSpace(this,f)
            %2Dkspace
            this.Reconstructor.setSignalKSpace(f);
        end
        function o=getSignalKSpace(this)
            o=this.Reconstructor.getSignalKSpace();
        end
        function o=getNoiseCovariance(this)
            o= this.Reconstructor.getNoiseCovariance();
        end
        
        
        
        function o=getReferenceImage(this)
            o=this.ReferenceImage;
        end
        
        
        function setReferenceImage(this,ref)
            this.ReferenceImage=ref;
        end
        
        function [SNR,noise]=getOutput(this)
            [SNR,noise]= this.calculate();
            
        end
        
        function [OUT]=fakesomedata(this,NR)
            %NR=this.getNumberOfPseudoReplicas();
            %this method works on 2d images only
            OUT=[];
            noisecov=this.getNoiseCovariance();
            %riccardo lattanzi
            [V,D] = eig(noisecov);
            corr_noise_factor = V*sqrt(D)*inv(V); %1 i
            
            %calculate the sensitivity maps once.
            K=this.Reconstructor.getSignalKSpace();
            
            
            for r=1:NR
                n=this.getPseudoNoise(size(K),corr_noise_factor);
                
                if(this.Reconstructor.isAccelerated)
                    n=n.*(K~=0); %also the noise is accelerated:)
                end
                
                OUT=cat(4,OUT,n+K);
            end
            
            %get back th old kspace
            this.Reconstructor.setSignalKSpace(K);
            
        end
        
        
        function [o, o2]=calculate(this)
            
            
            if(~isempty(this.getImageArray()))
                %using the imageArray
                this.logIt('SNR PMr caclulated with image Array, considering the first replica the  reference','ok');
                
                %in this case the reference image must be set otside the class or is taken
                %the mean of the stack
                
            else
                
                
                try
                    NR=this.getNumberOfPseudoReplicas();
                    %this method works on 2d images only
                    noisecov=this.getNoiseCovariance();
                    %riccardo lattanzi
                    [V,D] = eig(noisecov);
                    corr_noise_factor = V*sqrt(D)*inv(V); %1 i
                    
                    %calculate the sensitivity maps once.
                    
                    
                    
                    if ( this.Reconstructor.needsSensitivity())
                        S=this.Reconstructor.getCoilSensitivityMatrix();
                        this.Reconstructor.setCoilSensitivityMatrix(S);
                    end
                    
                    
                    %reconstruct the ref image
                    K=this.Reconstructor.getSignalKSpace();
                    
                    
                    
                    for r=1:NR
                        n=this.getPseudoNoise(size(K),corr_noise_factor);
                        
                        if(this.Reconstructor.isAccelerated)
                            n=n.*(K~=0); %also the noise is acceleraed:)
                            
                        end
                        
                        this.Reconstructor.setSignalKSpace(K+n);
                        this.add2DImage(this.Reconstructor.getOutput());
                    end
                    
                    %get back the original Kspace
                    this.Reconstructor.setSignalKSpace(K);
                    
                catch
                    this.logIt('cannot calculate Pseudo MR','ko');
                    
                end
            end
            o2=this.getPseudoMultipleReplicaNoiseImage();
            o=this. getPseudoMultipleReplicaSignalImage()./o2;
            this.SNR=o;
            this.add2DImagetoExport(o,'SNR')
            this.add2DImagetoExport(o2,'STD')
            this.logIt('PSEUDO MR calculated','ok');
        end
        % for pseudomr
        function pmr_image_noise=getPseudoMultipleReplicaNoiseImage(this)
            %std and mean in lattanzi are mean(abs()) and std(abs())
            pmr_image_stack=this.getImageArray();
            pmr_image_noise = std(abs(pmr_image_stack) + max(abs(pmr_image_stack(:))),[],3);
            pmr_image_noise(pmr_image_noise < eps) = 1;
            
        end
        
        function numerator=getPseudoMultipleReplicaSignalImage(this)
            % I can have a reference image
            % i can have a reconstructor
            %i can have a stack of images ()
            numerator=this.getReferenceImage();
            if (isempty(numerator))
                numerator=abs(this.Reconstructor.getOutput());
                this.setReferenceImage(numerator);
                if isempty(numerator)
                    numerator=this.getImageArrayMEAN();
                    this.setReferenceImage(numerator);
                    if (isempty(numerator))
                        this.logIt('couldn find a reference image','error')
                        numerator=[];
                    else
                        this.logIt(['reference image was set from the reconstructor class ' class(this.Reconstructor)],'warning')
                    end
                else
                    this.logIt(['reference image was set from the reconstructor class ' class(this.Reconstructor)],'warning')
                end
            else
                this.logIt('reference image was set outside the class','warning')
            end
        end
        %ovverride
        function exportResults(this,fn)
            O.version='CLOUDMR2DACM20190409';
            O.author='eros.montin@gmail.com';
            if isempty(this.Type)
                
                O.type='DATA';
            else
                O.type=this.Type;
            end
            
            
            
            if isempty(this.subType)
                
                O.subtype='';
            else
                O.subtype=this.subType;
            end
            
            %defined in every acm method
            b=this.Reconstructor.getParams();
            b.NR=this.NR;
            
            for t=1:size(this.Exporter,1)
                if(strcmp(this.Exporter{t,1},'image2D'))
                    im.slice=this.image2DtoJson(this.Exporter{t,3});
                    im.imageName=this.Exporter{t,2};
                    O.images(t)=im;
                    clear im;
                    
                end
                
                
            end
            
            O.settings=b;
            myjsonWrite(jsonencode(O),fn);
        end
        
        
        
        
        
        
        
        
    end
    
    
    methods (Static)
        
        function N=getFakeRawDataNoise(ksize,NC,corr_noise_factor)
            %ksize (freq,phase,coil)
            F=cm2DSignalToNoiseRatioPseudoMultipleReplicas();
            if nargin<3
                corr_noise_factor=F.getNoiseCorrelationFactor(NC);
            end
            N=F.getPseudoNoise(ksize,corr_noise_factor);
        end
        
        function corr_noise_factor=getNoiseCorrelationFactor(noisecov)
            [V,D] = eig(noisecov);
            corr_noise_factor = V*sqrt(D)*inv(V); %1 i
        end
        
        function gaussian_whitenoise=getPseudoNoise(msize,corr_noise_factor)
            %msize (freq,phase,coil)
            
            %dunno why... ask Riccardo
            N=sqrt(0.5)*(randn(msize)+1i*randn(msize));
            
            %             for n=1:msize(3)
            %                 X(n,:)=reshape(N(:,:,n),1,[]);
            %             end
            %
            %
            %             X = noise_corr_coeff * X;
            %
            %             for n=1:msize(3)
            %                 o(:,:,n)=reshape(X(n,:),msize(1),msize(2));
            %             end
            
            
            nrow=msize(1);
            ncol=msize(2);
            nchan=msize(3);
            
            gaussian_whitenoise = reshape(N,[nrow*ncol nchan]);
            gaussian_whitenoise = corr_noise_factor*(gaussian_whitenoise.');
            gaussian_whitenoise = reshape((gaussian_whitenoise.'),[nrow ncol nchan]);
            
            
            
            
            
            
            
        end
        
        
        
        
        
    end
    
end


