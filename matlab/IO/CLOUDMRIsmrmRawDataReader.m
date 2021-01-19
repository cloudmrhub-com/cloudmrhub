classdef CLOUDMRIsmrmRawDataReader<CLOUDMRRD
    
    properties
        Dset %H5 dataset
        Hdr %H5
        ICE_FACTOR=3200;
        DATA
        computePF=0
    end
    
    methods
        function this = CLOUDMRSiemensRawDataReader(f)
            %data are collected and if data is a structure the possible
            %(image) (noise) and refscan will have is own hdr
                   
            if nargin>0
                this.setFilename(f);
            end
            
            
        end
        
        
        function readDataset(this)
            if ~isempty(this.getFilename())
                this.Dset = ismrmrd.Dataset(this.getFilename(), 'dataset');
                fprintf(1,'\nread Dataset for file %s\n\n',this.getFilename())  ;
                
            else
                fprintf(1,'\nno filename set in %s\n\n',class(this));
                
            end
        end
        
        function readHdr(this)
            if isempty(this.Hdr)
                if ~isempty(this.Dset)
                    this.Hdr = ismrmrd.xml.deserialize(this.Dset.readxml);
                    fprintf(1,'\nread xmxl dor file %s\n\n',this.getFilename());
                    
                else
                    this.readDataset();
                    if ~isempty(this.Dset)
                        this.readHdr();
                    end
                    
                end
                
            end
        end
        
        
        function readKSpace(this)
            %%TODO   if cartesian
            
            this.readHdr();
            
            %% Encoding and reconstruction information
            % Matrix size
            
            enc=this.getEncodedSpaceInfo();
            enc_Nx = enc.Nx;
            enc_Ny = enc.Ny;
            enc_Nz = enc.Nz;
            
            rec=this.getReconSpaceInfo();
            
            rec_Nx = rec.Nx;
            rec_Ny = rec.Ny;
            rec_Nz = rec.Nz;
            
            
            nSlices = this.getNSlice();
            
            nCoils = this.getNCoils();
            
            nReps = this.getNRepetition();
            
            nContrasts = this.getNContrast();
            
            nAvg = this.getNAverage();
            
            %% Read all the data
            % Reading can be done one acquisition (or chunk) at a time,
            % but this is much faster for data sets that fit into RAM.
            D = this.Dset.readAcquisition();
            
            % Note: can select a single acquisition or header from the block, e.g.
            % acq = D.select(5);
            % hdr = D.head.select(5);
            % or you can work with them all together
            
            %% Ignore noise scans
            % TODO add a pre-whitening example
            % Find the first non-noise scan
            % This is how to check if a flag is set in the acquisition header
            isNoise = D.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
            %             firstScan = find(isNoise==0,1,'first');
            firstScan = 1;
            
            if firstScan > 1
                noise = D.select(1:firstScan-1);
            else
                noise = [];
            end
            meas  = D.select(firstScan:D.getNumber);
            clear D;
            
            %% Reconstruct images
            % Since the entire file is in memory we can use random access
            % Loop over repetitions, contrasts, slices
            nimages = 0;
            
            %move kspace in the phase direction
            
            E=this.getEncodedSpaceInfo();
            
            K = zeros(nAvg,nContrasts,nReps,enc_Nx, enc_Ny, nSlices, nCoils);
            Kmask = zeros(nAvg,nContrasts,nReps,enc_Nx, enc_Ny, nSlices, nCoils);
            
            PF=this.isPF();
            xx=this.getPF();
            for avg= 1:nAvg
                display(['average' num2str(avg)]);
                for rep = 1:nReps
                    display(['repetition' num2str(rep)]);
                    for contrast = 1:nContrasts
                        display(['contrast' num2str(contrast)]);
                        for slice = 1:nSlices
                            display(['slice' num2str(slice)]);
                            
                            %                        nimages = nimages + 1;
                            
                            % Initialize the K-space storage array
                            
                            % Select the appropriate measurements from the data
                            %                             acqs = find(  (meas.head.idx.average==(avg-1)) & (meas.head.idx.contrast==(contrast-1)) ...
                            %                                 & (meas.head.idx.repetition==(rep-1)) ...
                            %                                 & (meas.head.idx.slice==(slice-1)) );
                            
                            
                            acqs = find(  (meas.head.idx.average==(avg-1)) & (meas.head.idx.contrast==(contrast-1)) ...
                                & (meas.head.idx.repetition==(rep-1)) ...
                                & (meas.head.idx.slice==(slice-1)) );
                            
                            for p = 1:length(acqs)
                                ky = meas.head.idx.kspace_encode_step_1(acqs(p)) + 1; %matlab
                                kz = meas.head.idx.kspace_encode_step_2(acqs(p)) + 1; %matlab
                                
                                
                                K(avg,contrast,rep,:,ky,kz,:) = meas.data{acqs(p)};
                                Kmask(avg,contrast,rep,:,ky,kz,:) = ones(size(meas.data{acqs(p)}));
                            end
                            
                            
                            if(PF && this.computePF)
                                npf=this.getPF();                                
                                KK=reshape(this.resamplePF(reshape(K(avg,contrast,rep,:,this.getPF(),slice,:),enc_Nx, numel(npf),1,nCoils)),1,1,1,enc_Nx, enc_Ny, 1, nCoils);
                                K(avg,contrast,rep,:,:,slice,:)=KK;
                            end
                            
                        end
                        
                    end
                end
            end
            %resample the KSPACE
            
            %frequency oversampling
                        if(rec_Nx~=enc_Nx)
            K=this.filterOS(K);
                        end
            this.DATA=K;
            
            
        end
        
        
        function E=getEncodingLimits(this)
            this.readHdr();
            E.k1.max=this.Hdr.encoding.encodingLimits.kspace_encoding_step_1.maximum +1;
            E.k1.min=this.Hdr.encoding.encodingLimits.kspace_encoding_step_1.minimum +1;
            E.k1.center=this.Hdr.encoding.encodingLimits.kspace_encoding_step_1.center +1;
            
            E.k2.max=this.Hdr.encoding.encodingLimits.kspace_encoding_step_2.maximum +1;
            E.k2.min=this.Hdr.encoding.encodingLimits.kspace_encoding_step_2.minimum +1;
            E.k2.center=this.Hdr.encoding.encodingLimits.kspace_encoding_step_2.center +1;
            
            %             E.acceleration=
        end
        
        
        function O=resamplePF(this,K)
            if(this.isPF)
                %old=[1 size(K,1) size(K,1);1 size(K,2) size(K,2)];
                E=this.getEncodedSpaceInfo();
                %new=[1  size(K,1) E.Nx;1 size(K,2) E.Ny];
                RE=this.getEncodingLimits();
                
                %                 T= (RE.k1.max-RE.k1.center)+1;
                T=E.Ny-RE.k1.max;
                nx= this.getPF()+T;
                O=zeros(E.Nx,E.Ny,1,this.getNCoils); %posso farlo meglio
                O(:,nx,:,:)=K;
                
                
                [~, f] = pocs( permute(squeeze(O(:,:,1,:)),[3 1 2]), 200 );
                
                O=permute(reshape(f,this.getNCoils,E.Nx,E.Ny,1),[2 3 4 1]);
                
                
                
                %O(:,:,1,c)=resempleIT(squeeze(K(:,:,1,c)),old,new);
                
                
                
            end
        end
        
        
        
        function O=isPF(this)
            E=this.getEncodedSpaceInfo();
            
            Er=this.getEncodingLimits();
            
            if(E.Ny ~= Er.k1.max)
                O=1;
            else
                O=0;
            end
        end
        
        
        function E=getParallelInformation(this)
            this.readHdr();
            E.parallel=this.Hdr.encoding.parallelImaging;
            
        end
        
        
        
        function E=getPF(this)
            this.readHdr();
            E=[this.Hdr.encoding.encodingLimits.kspace_encoding_step_1.minimum : this.Hdr.encoding.encodingLimits.kspace_encoding_step_1.maximum ]+1;
            
        end
        
        
        
        
        
        function enc=getEncodedSpaceInfo(this)
            this.readHdr();
            
            %% Encoding and reconstruction information
            % Matrix size
            enc.Nx = this.Hdr.encoding.encodedSpace.matrixSize.x;
            enc.Ny = this.Hdr.encoding.encodedSpace.matrixSize.y;
            enc.Nz = this.Hdr.encoding.encodedSpace.matrixSize.z;
            
            % Field of View
            enc.FOVx = this.Hdr.encoding.encodedSpace.fieldOfView_mm.x;
            enc.FOVy = this.Hdr.encoding.encodedSpace.fieldOfView_mm.y;
            enc.FOVz = this.Hdr.encoding.encodedSpace.fieldOfView_mm.z;
        end
        
        
        function rec=getReconSpaceInfo(this)
            this.readHdr();
            %% Encoding and reconstruction information
            % Matrix size
            rec.Nx = this.Hdr.encoding.reconSpace.matrixSize.x;
            rec.Ny = this.Hdr.encoding.reconSpace.matrixSize.y;
            rec.Nz = this.Hdr.encoding.reconSpace.matrixSize.z;
            
            % Field of View
            rec.FOVx = this.Hdr.encoding.reconSpace.fieldOfView_mm.x;
            rec.FOVy = this.Hdr.encoding.reconSpace.fieldOfView_mm.y;
            rec.FOVz = this.Hdr.encoding.reconSpace.fieldOfView_mm.z;
            
        end
        
        function o=getTrajectory(this)
            this.readHdr();
            try
                o = this.Hdr.encoding.trajectory;
            catch
                fprintf(1,'any problems reading the slice number')
            end
        end
        
        
        
        function o=getNCoils(this)
            %ismrm
            this.readHdr();
            try
                o = this.Hdr.acquisitionSystemInformation.receiverChannels;
            catch
                fprintf(1,'problems reading the coils number')
                o = 1;
            end
        end
        
        function o=getNRepetition(this)
            %ismrm
            this.readHdr();
            try
                o = this.Hdr.encoding.encodingLimits.repetition.maximum + 1;
            catch
                fprintf(1,'any problems reading the repetition number')
                o = 1;
            end
        end
        
        
        
        function o=getNAverage(this)
            
            this.readHdr();
            try
                o = this.Hdr.encoding.encodingLimits.average.maximum + 1;
            catch
                fprintf(1,'any problems reading the averages number')
                o = 1;
            end
        end
        
        function o=getNContrast(this)
            
            this.readHdr();
            try
                o = this.Hdr.encoding.encodingLimits.contrast.maximum + 1;
            catch
                fprintf(1,'any problems reading the constrast number')
                o = 1;
            end
        end
        
        function o=getNSlice(this)
            
              this.readHdr();
            try
                o = this.Hdr.encoding.encodingLimits.slice.maximum + 1;
            catch
                fprintf(1,'any problems reading the slice number')
                o = 1;
            end
        end
        
    function o=getDset(this)
            o=this.Dset; %H5 dataset
        end
        
        function o=getHdr(this)
            o=this.Hdr; %H5 dataset
        end
        
        
        
      
        
        
        
        function o=readImageKSpace(this)
        
            if(isempty(this.DATA))
                this.readKSpace();
            end
                o=this.DATA;
            
            
            
        end
        
        
        function o=readNoiseKSpace(this)
            display('pay attention no kspace dat is readable in H5');
            o=this.readImageKSpace();
            
        end
        
        
        
        function o=filterOS(this,K)
          
            
            %K = zeros(nAvg,nContrasts,nReps,enc_Nx, enc_Ny, nSlices, nCoils);
          o=K(:,:,:,1:2:end,:,:,:);
            
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
    
    methods(Static)
        
         
        
       
        
    end
    
end
