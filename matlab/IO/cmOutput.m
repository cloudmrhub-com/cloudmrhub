classdef cmOutput<handle
    %main class of cloudmr
    %Main reconstrucion class
    %author:eros.montin@gmail.com
    %v012021
    %46&2 just ahead of me
    
    
    properties(Access=private)
        Exporter % a cell array of name and alues to be exported
        Log
        Type
        SubType
        OUTPUTLOGFILENAME
        OUTPUTFILENAME
    end
    
    
    
    methods(Access=protected)
        
        function addToExporter(this,type,name,value)
            this.Exporter=[this.Exporter;[{type},{name},{double(value)}]];
        end
        
    end
    methods
        function this = cmOutput()
            this.logIt(['class ' class(this) ' instantiated'],'start');
            try
                TurnOffWarnings();
                
            catch
                
            end
            
            
            
            
            
        end
        
        
        function setOutputLogFileName(this,L)
            this.OUTPUTLOGFILENAME=L;
        end
        
        function o=getOutputLogFileName(this)
            o=this.OUTPUTLOGFILENAME;
        end
        
        
        function setOutputFileName(this,L)
            this.OUTPUTFILENAME=L;
        end
        
        function o=getOutputFileName(this)
            o=this.OUTPUTFILENAME;
        end
        
        
        function setLog(this,L)
            this.Log=L;
        end
        
        function o=getLog(this)
            o=this.Log;
        end
        
        
        
        
        function appendLog(this,L)
            this.setLog(cat(2,this.getLog(),L));
        end
        
        function logIt(this,W,t)
            if nargin<3
                t='null';
            end
            j.time=datestr(clock);
            j.text=W;
            j.type=t;
            this.Log=[this.Log j];
            
        end
        
        
        
        
        
        
        function O=getResults(this)
            O=this.exportResults();
        end
        
        
        function O=getResultsAsImages(this)
            O=this.getImagesFromResuls(this.getResults());
        end
        
        
        
        
        
        function O=exportResults(this,fn)
            
            O.version='20201223';
            O.author='eros.montin@gmail.com';
            
            if isempty(this.Type)
                
                O.type='DATA';
            else
                O.type=this.Type;
            end
            
            
            
            if isempty(this.SubType)
                
                O.subtype='';
            else
                O.subtype=this.SubType;
            end
            
            
            
            
            for t=1:size(this.Exporter,1)
                if(strcmp(this.Exporter{t,1},'image2D'))
                    im.slice=this.image2DToJson(this.Exporter{t,3});
                    im.imageName=this.Exporter{t,2};
                    O.images(t)=im;
                    clear im;
                    
                end
                
                
            end
            if (nargin>1)
                myjsonWrite(jsonencode(O),fn);
            else
                try
                    myjsonWrite(jsonencode(O),this.getOutputFileName());
                catch
                    display('not saved to file just got back to the console');
                end
            end
        end
        
        
        
        
        
        
        function exportLog(this,fn)
            if(nargin>1)
                myjsonWrite(jsonencode(this.Log),fn);
            else
                try
                    myjsonWrite(jsonencode(this.Log),this.getOutputLogFileName());
                catch
                    display('error');
                end
                
            end
        end
        
        function whatHappened(this)
            %where's bugo?
            for t=1:numel(this.Log)
                display(this.Log(t));
            end
        end
        
        
        
        function errorMessage(this)
            this.logIt('ERROR','error');
        end
        
        function outputError(this,fn)
            this.errorMessage();
            this.exportLog(fn);
        end
        
        
        
        
        
        
        function add2DImagetoExport(this,im,name)
            %adds to results the image and its name
            % a.add2DImagetoExport(randn(20),'imname')
            if (isnumeric(im))
                this.addToExporter('image2D',name,im);
            else
                display('nothing added')
            end
        end
    end
    
    methods (Static)
        
        function Oout=image2DToJson(im)
            for sl=1:size(im,3)
                o=im(:,:,sl).';
                O.w=size(o,2);
                O.h=size(o,1);
                O.type= 'double complex';
                O.Vr=real(reshape(o,1,[]));
                O.Vi=imag(reshape(o,1,[]));
                Oout(sl)=O;
            end
        end
        
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
        
        
        
        function o=rescale01(bla)
            MI=min(bla(:));
            MA=max(bla(:));
            o= (bla(:) - MI)./ ( MA - MI );
        end
        
        
        
        
        function [O]=getJsonResultFromJsonFile(file)
            %from jsonfilename you get the struct of data
            O=readanddecodejson(file);
        end
        
        
        
        function [O]=getImagsesFromJsonResultFile(file)
            %from the json results file get the entire image
            %set in a struct variable
            data=CLOUDMROutput.getJsonResultFromJsonFile(file);
            O=CLOUDMROutput.getImagesFromResuls(data);
        end
        
        function [O] =getImagesFromResuls(data)
            %from the struct derived by the json results file got the entire image
            %set
            
            for imnumber=1:numel(data.images)
                h=data.images(imnumber).slice.h;
                w= data.images(imnumber).slice.w;
                NSL=numel(data.images(imnumber).slice);
                clear image_;
                %we switch i know!:)
                image_=NaN(w,h,NSL);
                O(imnumber).ImageName=data.images(imnumber).imageName;
                for slnumber=1:NSL
                    a=data.images(imnumber).slice(slnumber).Vr(:)+data.images(imnumber).slice(slnumber).Vi(:)*1i;
                    
                    image_(:,:,slnumber)=reshape(a,h,w).';
                end
                O(imnumber).image=image_;
                
                
            end
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
            nYsamp=prod(size(K))/(nX*nCoils);
            nAvg=size(K,1);
            nCnt=size(K,2);
            nRep=size(K,3);
            nPh=size(K,5);
            nSl=size(K,6);
            
            %nYsamp is number of actually
            acqblock = ismrmrd.Acquisition(nYsamp);
            
            
            % Set the header elements that don't change
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
            header.acquisitionSystemInformation.relativeReceiverNoiseBandwidth=0793;
            
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
        
        
        
        function TEST= test()
            %instantiate the class
            TEST=cmOutput();
            TEST.logIt('the test start now','start')
            T0=randn(5)+i*randn(5);
            T1=randn(5);
            TEST.logIt('created 2 results','ok')
            %add two results
            TEST.logIt('add first result','ok')
            try
                TEST.add2DImagetoExport(T0,'test0');
            catch
                TEST.errorMessage()
            end
            
            TEST.logIt('add second result','ok')
            try
                TEST.add2DImagetoExport(T1,'test0');
            catch
                TEST.errorMessage()
            end
            %get the results
            TEST.logIt('retrieve the results json','ok')
            try
                R= TEST.getResults();
            catch
                TEST.errorMessage();
            end
            TEST.logIt('results retrieved','ok')
            
            %get the results
            TEST.logIt('retrieve the results json','ok')
            
            TEST.logIt('retrieve the results as images','ok')
            try
                O=TEST.getResultsAsImages();
                TEST.logIt('retrieve the results json','ok')
            catch
                TEST.errorMessage();
            end
            %check the results
            TEST.logIt('tested the results ','?')
            try
                f=sum(abs(O(1).image(:)-T0(:)))+sum(abs(O(2).image(:)-T1(:)));
                TEST.logIt(['tested the results  and result error is :' num2str(f) ],'ok')
            catch
                TEST.errorMessage();
            end
            
            
            if(f==0)
                TEST.logIt(['tested the results  and result error is :' num2str(f) ],'ok')
            else
                TEST.errorMessage();
            end
            
            
            
            TEST.write2DCartesianKspacedatainISMRMRDv1(O(1).image,['~/test_mroptimum' num2str(uint16(100000*randn(1))) '.H5'])
            
            
            
            
            
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
        
        
        function [KOUT]=plotImageAfterTest(IM,tit)
            
            figure()
            subplot(121)
            imshow(IM,[]);
            title(tit);
            colorbar();
            subplot(122)
            imshow(abs(IM),[]);
            title([tit  ' abs']);
            colorbar();
            
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 0.98,[num2str(size(IM,1)) 'x' num2str(size(IM,2))])
            
            
        end
        
        
        
        function [KOUT]=plotTwoImagesAfterTest(IM1,IM2,tit1,tit2)
            
            figure()
            subplot(121)
            imshow(IM1,[]);
            title(tit1);
            colorbar();
            subplot(122)
            imshow(IM2,[]);
            title([tit2]);
            colorbar();
            
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 0.98,['IM1: ' num2str(size(IM1,1)) 'x' num2str(size(IM1,2))])
            text(0.5, 0.88,['IM2: ' num2str(size(IM2,1)) 'x' num2str(size(IM2,2))])
            
       
            RD=abs(IM1(:)-IM2(:))./IM1(:)*100;
                 G=find(~isinf(RD(:)) & ~isnan(RD(:)));
            text(0.5, 0.95,['%RD mean: ' num2str(nanmean(abs(RD(G)))) ', std: ' num2str(nanstd(abs(RD(G))))])
            
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
        
        
        function noise_coeff=calculateNoiseCoefficientsMatrix(noisecov)
            noise_coeff=zeros(size(noisecov));
            for itemp = 1:size(noisecov,1)
                for jtemp = 1:size(noisecov,1)
                    noise_coeff(itemp,jtemp) = noisecov(itemp,jtemp)/sqrt(noisecov(itemp,itemp)*noisecov(jtemp,jtemp));
                end
            end
        end
        
        
                function NK =resize2DMatrixbyInterpolationVoxelSpace(theK,x1,x2)
    
        [nX, nY]=size(theK);
       [newxg,newyg]=ndgrid(x1,x2);
       %old grid
       [oldxg,oldyg]=ndgrid([1:nX],[1:nY]);
       %scattered interpolation
        if(isreal(theK))
         F=scatteredInterpolant([oldxg(:) oldyg(:)],double(theK(:)));
        %composition of the channels
        NK=F(newxg,newyg)
        else
        FR=scatteredInterpolant([oldxg(:) oldyg(:)],double(real(theK(:))));
        FI=scatteredInterpolant([oldxg(:) oldyg(:)],double(imag(theK(:))));
        %composition of the channels
        NK=FR(newxg,newyg) +1i*FI(newxg,newyg);
        end    
        
                end
        
    end
end

