classdef cm2DRecon<cmOutput
    %Main reconstrucion class 
    %author:eros.montin@gmail.com
    %v012021
    %46&2 just ahead of me
    
    properties(Access=protected)
        SignalKSpace %freq,phase,coils
        SignalNCoils %freq,phase,coils
        SignalSize %freq,phase
        NoiseKSpace%freq,phase,coils
        NoiseNCoils %freq,phase,coils
        NoiseSize %freq,phase
        NoiseCovariance %ncoils.ncoils
        InverseNoiseCovariance %ncoils.ncoils
        SignalPrewhitened %freq,phase,coils      
        HasSensitivity
        HasAcceleration
        NoiseBandWidth
       end
    
    
    
    methods(Access=protected)
        function setHasSensitivity(this,s)
                this.HasSensitivity=s;
        end
        
            function setHasAcceleration(this,s)
                this.HasAcceleration=s;
            end    

    end    
    methods
        function this = cm2DRecon()
            %the class expects a 3D matrix composed by a tile of 2D kspaces (fxpxncoils) of a signal and a
            %noise covariance matrix.
            %if a third element is given that's the prewhitened signal
            %kspace
%             log=[];
%             if nargin>0
%                 log=[log 'a signal matrix of size: '  num2str(size(s)) ];
%                 this.setSignalKSpace(s);
%             end
%             
%             
%             if nargin>1
%                 if(~isempty(n))
%                 this.setNoiseKSpace(n)
%                 log=[log ' a noise matrix of size: ' num2str(size(n))];
%                 end
%             end
%             
%             
%             if nargin>2
%                 if(~isempty(nc))
%                 this.setNoiseCovariance(nc)
%                 log=[log ' a covariance matrix of size: ' num2str(size(nc))];
%                 end
%             end
%             
%             if nargin>3
%                 if(~isempty(pws))
%                 this.setPrewhitenedSignal(pws)
%                 log=[log ' a covariance matrix of size: ' num2str(size(pws))];
%                 end
%             end
%             
%             if isempty(log)
%                 this.logIt(['instantiated without parameters'  log],'ok');
%             else
%                 
%             this.logIt(['instantiated with' log],'ok');
%             end
         end
        
               function TEST=test(this)
         TEST=this.therecon();
         TEST.plotImageAfterTest(TEST.getOutput(),'recon');
         TEST.whatHappened();
               end
         
        function setSignalKSpace(this,f)
            %2Dkspace
            
            this.SignalKSpace=f;
            this.logIt('set signal kspace','inputfiles')
            this.setSignalNCoils(size(f,3));
            this.logIt('set signal kspace','inputfiles')
            this.setSignalSize([size(f,1) size(f,2)]);

            %changed the signal Prewhitened is no more available
            this.setPrewhitenedSignal([]);

        end

        
        function o=getSignalKSpace(this)
            o=this.SignalKSpace;
        end
        
   function setNoiseKSpace(this,f)
            % 2DKspace
            this.NoiseKSpace=f;
            this.setNoiseNCoils(size(f,3));
        end
        
        function o =getNoiseKSpace(this)
            % 2DKspace
            o = this.NoiseKSpace;
        end
        
        
          function o=getNoiseNCoils(this)
            o=this.NoiseNCoils;
          end
        
               function o=setNoiseNCoils(this,o)
            this.NoiseNCoils=o;
               end
        
        function o=getNoiseBandWidth(this)
            if(isempty( this.NoiseBandWidth))
                noise_bandwidth = mrir_noise_bandwidth(this.getNoiseKSpace(),0);
                
                this.NoiseBandWidth=noise_bandwidth;
            end
            
            o=this.NoiseBandWidth;
            
        end
        
        
        
        function[h]=plotTheReconstruciton(this)
            A=this.getOutput();
            imshow(abs(A),[]);
        end
        
        
              function setSignalNCoils(this,n)
            this.SignalNCoils=n;
            
        end
        
        
        function o= getSignalNCoils(this)
            o=this.SignalNCoils;
        end
        
        
              function setSignalSize(this,ss)
                  %array with 2 components
            this.SignalSize=ss;
            this.logIt('signal matrix set','ok');
        end
        
        
        function o= getSignalSize(this)
            o=this.SignalSize;
        end
        
        function o=getNoiseCoefficient(this)
            o=this.calculateNoiseCoefficientsMatrix(this.getNoiseCovariance);
        end
        
        function Rn=getNoiseCovariance(this)
            this.logIt('get Covariance Matrix','?');
        Rn=this.NoiseCovariance;
        
        if (isempty(Rn))
            this.logIt('Covariance Matrix need to be calculated','?');
            n=this.getNoiseKSpace();
            if(isempty(n))
                this.logIt('need to set a Covariance Matrix or a noise KSpace datas','ok');
            else
            Rn=this.calculateCovarianceMatrix(this.getNoiseKSpace(),this.getNoiseBandWidth());
            this.setNoiseCovariance(Rn)
            this.logIt('Covariance Matrix calculated','ok');
            end
        else
            this.logIt('Covariance Matrix retrieved','ok');
        end
        
        
        end

         function setNoiseCovariance(this,nc)
            this.NoiseCovariance=nc;    
            this.logIt('noise covariance matrix set','ok');
            if(isempty(this.getInverseNoiseCovariance()))
                this.setInverseNoiseCovariance(inv(nc));
            end
         end
        
         
         
        function o=getInverseNoiseCovariance(this)
        o=this.InverseNoiseCovariance;
        end

         function setInverseNoiseCovariance(this,inc)
            this.InverseNoiseCovariance=inc;    
            this.logIt('inverse noise covariance matrix set','ok');
         end
         
         
                 function invRn=getInverseNoiseCovariancePrewhithened(this)
                Rn = this.getNoiseCovariancePrewhithened();
                invRn = inv(Rn);
         
                 end

        
       function Rn=getNoiseCovariancePrewhithened(this)
                Rn = eye(this.getSignalNCoils);
        end
         
         
         
            
            
         
        
        
         function pw=getPrewhitenedSignal(this)
            if(isempty(this.SignalPrewhitened))
             pw=this.prewhiteningSignal(this.getSignalKSpace(), this.getNoiseCovariance());
             this.setPrewhitenedSignal(pw);
            else
                pw=this.SignalPrewhitened;
            end
         end
         
           function setPrewhitenedSignal(this,pw)
            this.SignalPrewhitened=pw;    
            this.logIt('prewhitened signal matrix set','ok');
            if(isempty(this.getSignalSize))
            this.setSignalSize([size(pw,1) size(pw,2)]);
            end            
            if(isempty(this.getSignalNCoils))
            this.setSignalNCoils(size(pw,3));
            end
           end

        
        
        
         function o=getHasSensitivity(this)
                o=this.HasSensitivity();
        end
        
            function o=getHasAcceleration(this)
                o=this.HasAcceleration();
            end
            
            
        function o=needsSensitivity(this)
            o= this.getHasSensitivity();     
        end
        
        
        
                   function o=isAccelerated(this)
            o= this.getHasAcceleration();
                
            end
            
            
    end
    
    methods(Static)
           function o=get2DKSIFFT(k)
            % reconstruct individual coils' images and apply FFT scale factor
            % iFFT scales data by 1/sqrt(N), so I need to scale back, as the noise covariance
            % matrix is calculated in k-space (not Fourier transformed)
            
            o=MRifft(k,[1,2])*sqrt(size(k,1)*size(k,2));
        end
        
        
        function o=bartMy2DKSpace(K,cartesianFLAG)
            
            if nargin>1
                if(noncartesianFLAG)
                    display('not yet implemented')
                else
                    o(:,:,1,:)=K;
                end
            else
                if nargin>0
                    o(:,:,1,:)=K;
                else
                    o=[];
                    display('nothing to display');
                    
                end
            end
            
        end
        
        
        function o=debartMy2DKSpace(K,cartesianFLAG)
            o=[];
            if nargin>1
                if(noncartesianFLAG)
                    display('not yet implemented')
                else
                    o=squeeze(K);
                end
            else
                if nargin>0
                    o=squeeze(K);
                else
                    
                    display('nothing to display');
                    
                end
            end
            o=double(o);
        end
        
        function pw_signalrawdata=prewhiteningSignal(signalrawdata,psi)
            %static cm2DRecon.prewhiteningSignal(signalrawdata,psi)
            L = chol(psi,'lower');
            L_inv = inv(L);
            nc=size(signalrawdata,3);
            nf=size(signalrawdata,1);
            np=size(signalrawdata,2);
            pw_signalrawdata = permute(signalrawdata,[3,2,1]); % [ncoil nphase nfreq]
            pw_signalrawdata = L_inv*pw_signalrawdata(:,:);
            pw_signalrawdata = reshape(pw_signalrawdata,nc,np,nf);
            pw_signalrawdata = permute(pw_signalrawdata,[3 2 1]); % [nfreq nphase ncoil]
            
        end
        
%         https://github.com/SyneRBI/tools/blob/fd3404ea96945f75a999bbabd00bb9ab2ea1a91b/gen_us_data.m
        function write2DCartesianKSpaceDataIsmrmDv1Slice(K7,filename)
            % It is very slow to append one acquisition at a time, so we're going
% to append a block of acquisitions at a time.
% In this case, we'll do it one repetition at a time to show off this
% feature.  Each block has nY aquisitions

            nX=size(K7,4);
            nCoils=size(K7,7);
            
            nAvg=size(K7,1);
            nCnt=size(K7,2);
            nReps=size(K7,3);
            nPh=size(K7,5);
            nSl=size(K7,6);
            nY=nPh;
            dset = ismrmrd.Dataset(filename);

            %transform in this space
            %L=permute(K7,[ 4,5,7,3,1,2,6]);
            
            %K =L(:,:,:,:,1,1,1);
acqblock = ismrmrd.Acquisition(nY);

% Set the header elements that don't change
acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = nX;
acqblock.head.center_sample(:) = floor(nX/2);
acqblock.head.active_channels(:) = nCoils;
acqblock.head.read_dir  = repmat([1 0 0]',[1 nY]);
acqblock.head.phase_dir = repmat([0 1 0]',[1 nY]);
acqblock.head.slice_dir = repmat([0 0 1]',[1 nY]);

% Loop over the acquisitions, set the header, set the data and append
for rep = 1:nReps
    for acqno = 1:nY
        
        % Set the header elements that change from acquisition to the next
        % c-style counting
        acqblock.head.scan_counter(acqno) = (rep-1)*nY + acqno-1;
        acqblock.head.idx.kspace_encode_step_1(acqno) = acqno-1; 
        acqblock.head.idx.repetition(acqno) = rep - 1;
        
        % Set the flags
        acqblock.head.flagClearAll(acqno);
        if acqno == 1
            acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', acqno);
            acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', acqno);
            acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', acqno);
        elseif acqno==size(K,2)
            acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
            acqblock.head.flagSet('ACQ_LAST_IN_SLICE', acqno);
            acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', acqno);
        end
        
        % fill the data
        acqblock.data{acqno} = squeeze(K(:,acqno,:,rep));
    end

    % Append the acquisition block
    dset.appendAcquisition(acqblock);
        
end % rep loop


%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill the xml header %
%%%%%%%%%%%%%%%%%%%%%%%%
% We create a matlab struct and then serialize it to xml.
% Look at the xml schema to see what the field names should be

header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 128000000; % 3T

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.systemVendor = 'ISMRMRD Labs';
header.acquisitionSystemInformation.systemModel = 'Virtual Scanner';
header.acquisitionSystemInformation.receiverChannels = nCoils;

% The Encoding (Required)
header.encoding.trajectory = 'cartesian';
header.encoding.encodedSpace.fieldOfView_mm.x = 256;
header.encoding.encodedSpace.fieldOfView_mm.y = 256;
header.encoding.encodedSpace.fieldOfView_mm.z = 5;
header.encoding.encodedSpace.matrixSize.x = size(K,1);
header.encoding.encodedSpace.matrixSize.y = size(K,2);
header.encoding.encodedSpace.matrixSize.z = 1;
% Recon Space
% (in this case same as encoding space)
header.encoding.reconSpace = header.encoding.encodedSpace;
% Encoding Limits
header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_0.maximum = size(K,1)-1;
header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(size(K,1)/2);
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = size(K,2)-1;
header.encoding.encodingLimits.kspace_encoding_step_1.center = floor(size(K,2)/2);
header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = nReps-1;
header.encoding.encodingLimits.repetition.center = 0;

%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

%% Write the dataset
dset.close();
        end
        
        
        function write2DCartesianKspacedatainISMRMRDv1(K,filename)
            % It is very slow to append one acquisition at a time, so we're going
            % to append a block of acquisitions at a time.
            % In this case, we'll do it one repetition at a time to show off this
            % feature.  Each block has nYsamp aquisitions
            
            
            
            
            
            %
            % f='meas_MID00024_FID188178_Multislice.dat';
            %
            % F=CLOUDMRRD(f);
            
            
            % K=F.getNoiseKSpace();
            %
            %
            %
            % size(K)
            % 1     1     1    96    96     5    16
            %   O={'1: Average',
            %       '2: contrast',
            %       '3: repetition',
            %       '4: Frequency Encode',
            %       '5: Phase Encode',
            %       '6: Slice',
            %       '7: Coils'};
            
            
            
            
            
            % filename ='test.H5';
            dset = ismrmrd.Dataset(filename);
            
            
            nX=size(K,4);
            nCoils=size(K,7);
            
            nAvg=size(K,1);
            nCnt=size(K,2);
            nRep=size(K,3);
            nPh=size(K,5);
            nSl=size(K,6);
            
            nYsamp=nPh;% ERROR nYsamp=prod(size(K))/(nX*nCoils);
            %nYsamp is number of actually
            acqblock = ismrmrd.Acquisition(nYsamp);
            
            
            % Set the header elements that doesn't change
            acqblock.head.version(:) = 1;
            acqblock.head.number_of_samples(:) = nX;
            acqblock.head.center_sample(:) = floor(nX/2);
            acqblock.head.active_channels(:) = nCoils;
            acqblock.head.available_channels(:) =nCoils;
            acqblock.head.read_dir  = repmat([1 0 0]',[1 nYsamp]);
            acqblock.head.phase_dir = repmat([0 1 0]',[1 nYsamp]);
            acqblock.head.slice_dir = repmat([0 0 1]',[1 nYsamp]);
            
            
            
            
            
            counter=0;
            
            for avg=1:nAvg
                
                
                
                
                for cnt=1:nCnt
                    for rep = 1:nRep
                        
                        for sl = 1:nSl
                            for p = 1:nPh
                                counter=counter+1;
                                
                                if (p==1)
                                    acqblock.head.flagClearAll(counter);
                                    acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1',counter);
                                    
                                    if avg == 1
                                        acqblock.head.flagSet('ACQ_FIRST_IN_AVERAGE',counter);
                                    end
                                    if rep == 1
                                        acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION',counter);
                                    end
                                    if cnt == 1
                                        acqblock.head.flagSet('ACQ_FIRST_IN_CONTRAST',counter);
                                    end
                                    if sl == 1
                                        acqblock.head.flagSet('ACQ_FIRST_IN_SLICE',counter);
                                    end
                                    
                                    
                                    
                                end
                                
                                % Set the header elements that change from acquisition to the next
                                % c-style counting
                                acqblock.head.scan_counter(counter) = counter;
                                % Note next entry is k-space encoded line number (not acqno which
                                % is just the sequential acquisition number)
                                acqblock.head.idx.kspace_encode_step_1(counter) = p-1;
                                acqblock.head.idx.repetition(counter) = rep-1;
                                acqblock.head.idx.average(counter) = avg-1;
                                
                                acqblock.head.idx.slice(counter)=sl-1;
                                acqblock.head.idx.contrast(counter)=cnt-1;
                                acqblock.head.idx.kspace_encode_step_2(counter)=0;
                                acqblock.head.idx.segment(counter)=0;
                                acqblock.head.idx.set(counter)=0;
                                
                                % fill the data
                                acqblock.data{counter} = squeeze(K(avg,cnt,rep,:,p,sl,:));
                                % Append the acquisition block
                                
                                
                                if counter==nPh
                                    
                                    
                                    acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1',counter);
                                    
                                    if avg == nAvg
                                        acqblock.head.flagSet('ACQ_LAST_IN_AVERAGE',counter);
                                    end
                                    if rep == nRep
                                        acqblock.head.flagSet('ACQ_LAST_IN_REPETITION',counter);
                                    end
                                    if cnt == nCnt
                                        acqblock.head.flagSet('ACQ_LAST_IN_CONTRAST',counter);
                                    end
                                    
                                    if sl == nSl
                                        acqblock.head.flagSet('ACQ_LAST_IN_SLICE',counter);
                                    end
                                    
                                    
                                    
                                    
                                end
                                
                                
                                
                                
                            end
                            
                            if sl ==nSl
                                acqblock.head.flagSet('ACQ_LAST_IN_SLICE',counter);
                            end
                        end
                    end
                end
            end
            acqblock.head.flagSet('ACQ_LAST_IN_MEASUREMENT',counter);
            
            dset.appendAcquisition(acqblock);
            
            
            %get the resolution ffrom somewhere
            
            r=[1 1 1];
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %% Fill the xml header %
            %%%%%%%%%%%%%%%%%%%%%%%%
            % We create a matlab struct and then serialize it to xml.
            % Look at the xml schema to see what the field names should be
            
            header = [];
            
            % Experimental Conditions (Required)
            header.experimentalConditions.H1resonanceFrequency_Hz = 128000000; % 3T
            
            % Acquisition System Information (Optional)
            header.acquisitionSystemInformation.systemVendor = 'CLOUDMR www.cloudmrhub.com';
            header.acquisitionSystemInformation.systemModel = 'CM scanner v01';
            header.acquisitionSystemInformation.receiverChannels = nCoils;
            header.acquisitionSystemInformation.systemFieldStrength_T=2.893620;
            header.acquisitionSystemInformation.relativeReceiverNoiseBandwidth=0.793;
            
            % The Encoding (Required)
            header.encoding.trajectory = 'cartesian';
            header.encoding.encodedSpace.fieldOfView_mm.x = nX*r(1);
            header.encoding.encodedSpace.fieldOfView_mm.y = nPh*r(2);
            header.encoding.encodedSpace.fieldOfView_mm.z = nSl*r(3);
            header.encoding.encodedSpace.matrixSize.x = nX;
            header.encoding.encodedSpace.matrixSize.y = nPh;
            header.encoding.encodedSpace.matrixSize.z = 1;
            % Recon Space
            % (in this case same as encoding space)
            header.encoding.reconSpace = header.encoding.encodedSpace;
            % Encoding Limits
            header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
            header.encoding.encodingLimits.kspace_encoding_step_0.maximum = nX-1;
            header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(nX/2);
            header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
            header.encoding.encodingLimits.kspace_encoding_step_1.maximum = nPh-1;
            header.encoding.encodingLimits.kspace_encoding_step_1.center = floor(nPh/2);
            
            header.encoding.encodingLimits.kspace_encoding_step_2.minimum = 0;
            header.encoding.encodingLimits.kspace_encoding_step_2.maximum = 0;
            header.encoding.encodingLimits.kspace_encoding_step_2.center = 0;
            
            header.encoding.encodingLimits.repetition.minimum = 0;
            header.encoding.encodingLimits.repetition.maximum = nRep-1;
            header.encoding.encodingLimits.repetition.center = 0;
            
            header.encoding.encodingLimits.average.minimum = 0;
            header.encoding.encodingLimits.average.maximum = nAvg-1;
            header.encoding.encodingLimits.average.center = 0;
            
            
            header.encoding.encodingLimits.phase.minimum = 0;
            header.encoding.encodingLimits.phase.maximum = nPh-1;
            header.encoding.encodingLimits.phase.center = 0;
            
            
            header.encoding.encodingLimits.contrast.minimum = 0;
            header.encoding.encodingLimits.contrast.maximum = nCnt-1;
            header.encoding.encodingLimits.contrast.center = 0;
            
            header.encoding.encodingLimits.slice.minimum = 0;
            header.encoding.encodingLimits.slice.maximum = nSl-1;
            header.encoding.encodingLimits.slice.center = 0;
            
            header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 = 1 ;
            header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2 = 1 ;
            % header.encoding.parallelImaging.calibrationMode = 'embedded' ;
            
            % Commented code below appears not necessary - saw this parameter after converting
            % a scanner file using siemens_to_ismrmrd
            % header.userParameters.userParameterLong.name = 'EmbeddedRefLinesE1' ;
            % header.userParameters.userParameterLong.value = ACShw *2  ;
            
            %% Serialize and write to the data set
            xmlstring = ismrmrd.xml.serialize(header);
            dset.writexml(xmlstring);
            
            %% Write the dataset
            dset.close();
            
        end
        
        
        
        
        function [K,nc]=getMRoptimumTestData()
            %create a fake kspace
            k=imread('fabric.png');
            k=k(:,1:480,:);
            K=MRifft(k,[1,2]);
            %create a fake noise covariance matrix
            nc=eye(3).*rand(3)/1000000;
        end
        
        
          function [KOUT]=shrinktoaliasedmatrix_2d(K,R1,R2)
            %K is a kspace 2d data (frequency,phase,coil)
            % phaseacceleration
            % set ACS lines for mSense simulation (fully sampled central k-space
            % region)
            
            s=size(K);
            nf=s(1);
            np=s(2);
            KOUT=K(1:floor(nf/R1),1:floor(np/R2),:);
            
            
          end
        
          function Rn =calculateCovarianceMatrix(noise,bandwidth)
            %freq,phase,coils
            %update 09/24/2020 Riccardo Lattanzi
            
            nchan = size(noise,3);
            if (nargin<2)
                bandwidth=1;
            end
            
            Rn = zeros(nchan);
            noise_samples = reshape(noise,[size(noise,1)*size(noise,2) nchan]);
            noise_samples = noise_samples.';
            Rn = 1/(2*size(noise_samples,2))*(noise_samples*noise_samples');
            Rn = Rn/bandwidth;
        end
        
        
 
        
        
    end
    
end



