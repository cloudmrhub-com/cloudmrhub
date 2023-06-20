classdef cm2DReconWithSensitivity<cm2DRecon
    %CoilSenseitivityMatrix
    %
    %method
    %CoilSensitivityMatrixCalculationMethod, a string {'simplesense','adaptive','zerofilling','bodycoil'}
    %To get SNR we neeed to set S (sensitivity) and Noisecovariance matrix
    %v24122020
    properties%(Access=private)
        CoilSensitivityMatrix %Smap
        CoilSensitivityMatrixSourcePrewhitened %the matrix for the sensitivity calculation to be set prewhitened
        CoilSensitivityMatrixSource %the matrix for the sensitivity calculation to be set
        CoilSensitivityMatrixCalculationMethod %'simplesense' ,'adaptive','espirit'
        CoilSensitivityMatrixSourceNCoils %number of coils set with the source
        CoilSensitivityMatrixSourceSmooth=false
        MaskCoilSensitivityMatrix=true
        
    end
    methods%(Access=protected)
        function this=cm2DReconWithSensitivity()
            this.setHasAcceleration(0);
            this.setHasSensitivity(1);
        end
        function coilsens_set=getCoilSensitivityMatrixForFullySampledKSpaceEspirit(this)
            %in this case sens_matrix is the full image
            coilsens_set = espirit_sensitivitymap(this.getCoilSensitivityMatrixSourcePrewhitened());
            this.logIt(['sensitivity map calculated as espirit'],'ok');
            
        end
        function coilsens_set=getCoilSensitivityMatrixForFullySampledKSpaceSimpleSense(this)
            %in this case sens_matrix is the full image
            nc=this.getSignalNCoils;
            sensmap_temp = MRifft(this.getSignalKSpace(),[1,2]);
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
        function coilsens_set=getCoilSensitivityMatrixForFullySampledKSpaceAdapt2D(this)
            %check if it must be prewhitened
            recon=cm2DReconAdapt();
            recon.setPrewhitenedSignal(this.getCoilSensitivityMatrixSourcePrewhitened());
            coilsens_set  = recon.getOutputSensitivityMap();
            this.logIt(['sensitivity map calculated as adaptive'],'ok');
            %[~,coilsens_set] = adapt_array_2d(MRifft(this.getCoilSensitivityMatrixPrewhitened(),[1,2]),eye(nv));
        end
        function coilsens_set=getCoilSensitivityMatrixForFullySampledKSpaceBodyCoil(this)
            %not prewhitened because the number of coils do not match the
            %signal usually 2 coils
            sens_matrix=this.getCoilSensitivityMatrixSource();
            %in this case sens_matrix is the BC
            reference_image=MRifft(this.getPrewhitenedSignal(),[1,2]);
            switch( this.getCoilSensitivityMatrixSourceNCoils)
                case 1
                    coilsens_set = reference_image/repmat(sens_matrix,[1 1 nchan]);
                case 2
                    %old                                     BCIM=(sens_matrix(:,:,1)+1i*(sens_matrix(:,:,2)))/sqrt(2);
                    center_phase_1 = angle(sens_matrix(floor(end/2),floor(end/2),1));
                    center_phase_2 = angle(sens_matrix(floor(end/2),floor(end/2),2));
                    
                    image_1 = sens_matrix(:,:,1)*exp(-1i*center_phase_1);
                    image_2 = sens_matrix(:,:,2)*exp(-1i*center_phase_2);
                    BCIM =abs(image_1 + image_2)/sqrt(2);
                    coilsens_set = reference_image./repmat(abs(BCIM),[1 1 nchan]);
                otherwise
                    %in this case 1th kspace of the bodycoil
                    coilsens_set = reference_image./repmat(sens_matrix(:,:,1),[1 1 nchan]);
            end
            %                             coilsens_set(coilsens_set>1)=1;
            this.logIt(['sensitivity map calculated as BodyCoil'],'ok');
            
            for na=1:size(coilsens_set,1)
                for nb=1:size(coilsens_set,2)
                    V=[coilsens_set(na,nb,:)];
                    PMax=prctile(reshape(V,[],1),100);
                    PMin=prctile(reshape(V,[],1),0);
                    coilsens_set(na,nb,:)= (V-PMin)./PMax;
                end
            end
        end
        function coilsens_set=calculateCoilSensitivityMatrixFullySampled(this)
            
            %                                   if (~isempty(this.CoilSensitivityMatrixSourcePrewhitened))
            %                     this.logIt(['a Source Coil sensitivity map has been correctly set so i can calculate the senstivity map'],'ok');
            switch(lower(this.CoilSensitivityMatrixCalculationMethod))
                case {'espirit'}
                    this.logIt(['sensitivity map calculated as espirit'],'?');
                    coilsens_set=this.getCoilSensitivityMatrixForFullySampledKSpaceEspirit();
                case {'simplesense','internal reference','inner'}
                    this.logIt(['sensitivity map calculated as simplesense'],'?');
                    coilsens_set=this.getCoilSensitivityMatrixForFullySampledKSpaceSimpleSense();
                case 'adaptive'
                    this.logIt(['sensitivity map calculated as aspative'],'?');
                    coilsens_set = getCoilSensitivityMatrixForFullySampledKSpaceAdapt2D();
                case 'bodycoil'
                    this.logIt(['sensitivity map calculated as aspative'],'?');
                    coilsens_set = getCoilSensitivityMatrixForFullySampledKSpaceBodyCoil();
            end
        end
        %                   end
        
        %this function should be overrided on othe
        %implementations e.g. when you are using accelerated
        %methods
        function coilsens_set=calculateCoilSensitivityMatrix(this)
            %coil sensiticity matrix calculated by prewhitening the singal
            %and then obtined the matrix
            coilsens_set=this.calculateCoilSensitivityMatrixFullySampled();
            
        end
    end
    methods
        function setCoilSensitivityMatrix(this,S)
            this.CoilSensitivityMatrix=S;
        end
        function resetCoilSensitivityMatrix(this)
            this.setSensitivityMatrix([]);
        end
        function coilsens_set=getCoilSensitivityMatrix(this)
            %the sensitivity matrix source is fully sampled!!
            this.logIt(['sensitivity has been requested'],'ok');
            if (isempty(this.CoilSensitivityMatrix)) %we don't have a coil sensitivity matrix
                this.logIt(['Coil sensitivity map has never been calculated'],'ok');
                coilsens_set=this.calculateCoilSensitivityMatrix();
                if (this.getCoilSensitivityMatrixSourceSmooth())
                    for aa=1:size(coilsens_set,3)
                        coilsens_set(:,:,aa)=medfilt2(real(coilsens_set(:,:,aa)),[3 3],'symmetric')+1i*medfilt2(imag(coilsens_set(:,:,aa)),[3 3],'symmetric');
                    end
                    
                end
                this.CoilSensitivityMatrix=coilsens_set;
                this.logIt(['start sensitivity map export'],'ok');
            else
                coilsens_set=this.CoilSensitivityMatrix;
                
            end
        end
        function setCoilSensitivityMatrixSourcePrewhitened(this,x)
            this.CoilSensitivityMatrixSourcePrewhitened=x;
            if isempty(this.getCoilSensitivityMatrixSourceNCoils())
                this.setCoilSensitivityMatrixSourceNCoils(size(x,3));
            end
        end
        function pw_S=getCoilSensitivityMatrixSourcePrewhitened(this)
            this.logIt('is there a prewhitened Source Sensitivity Matrix?','?')
            if(isempty(this.CoilSensitivityMatrixSourcePrewhitened))
                this.logIt('there is no prewhitened Sensitivity Matrix then calculate it','go')
                %get sens
                this.logIt('get the Sensitivity Matrix','?')
                S=this.getCoilSensitivityMatrixSource();
                this.logIt('get back to the calculation of the prewhitened Sensitivity Matrix','go')
                try
                    Rn= this.getNoiseCovariance();
                catch
                    
                end
                
                if ~isempty(Rn)
                    pw_S=this.prewhiteningSignal(S,Rn);
                    this.logIt('we got a prewhitened Sensitivity Matrix','ok')
                    this.setCoilSensitivityMatrixSourcePrewhitened(pw_S);
                else
                    this.errorMessage();
                    return
                end
                
            else
                this.logIt('there is a prewhitened Sensitivity Matrix and i am serving it to you','ok')
                pw_S=this.CoilSensitivityMatrixSourcePrewhitened;
            end
            
        end
        function setCoilSensitivityMatrixSource(this,IM)
            %expect a 2D coil sens map
            this.CoilSensitivityMatrixSource=IM;
            this.setCoilSensitivityMatrixSourceNCoils(size(IM,3));
        end
        
        
        function s=getCoilSensitivityMatrixSource(this)
            s=this.CoilSensitivityMatrixSource();
        end
        
        
        
                function setMaskCoilSensitivityMatrix(this,x)
            this.MaskCoilSensitivityMatrix=x;
        end
        
        function x=getMaskCoilSensitivityMatrix(this)
            x=this.MaskCoilSensitivityMatrix;
        end

        
        
        function setCoilSensitivityMatrixCalculationMethod(this,x)
            this.CoilSensitivityMatrixCalculationMethod=x;
        end
        
        function x=getCoilSensitivityMatrixCalculationMethod(this)
            x=this.CoilSensitivityMatrixCalculationMethod;
        end
        function o=getCoilSensitivityMatrixSourceNCoils(this)
            o=this.CoilSensitivityMatrixSourceNCoils();
        end
        
        
        function setCoilSensitivityMatrixSourceNCoils(this,o)
            this.CoilSensitivityMatrixSourceNCoils=o;
        end
        
        
        function setCoilSensitivityMatrixSourceSmooth(this,S)
            this.CoilSensitivityMatrixSourceSmooth=S;
        end
        
        function o=getCoilSensitivityMatrixSourceSmooth(this)
            o=this.CoilSensitivityMatrixSourceSmooth;
        end
        function o=testSensitivityMatrixvalidity(this)
            ss= size(this.getCoilSensitivityMatrixPrewhitened());
            is= size(this.getPrewhitenedSignal());
            
            if sum(ss(1:2)== is(1:2))==2
                o=true;
                this.logIt('sensitivity is valid','ok')
            else
                o=false;
                this.logIt('sensitivity is not valid (size differencies)','ok')
                
            end
        end
    end
end

