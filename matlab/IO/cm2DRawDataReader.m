classdef cm2DRawDataReader<cmOutput
    %cm2DRawDataReader is the basic data structure for the KSpace data read
    %    
    %  
    %   Eros Montin eros.montin@gmail.com
    % 2021/01
    %properties(Access=private) 7D kspace
    
    properties
        FileName
        Reader % a kaspace reader class
        RawDataImageKSpace
        RawDataNoiseKSpace
        RawDataRefsKSpace
        AverageRawDataImageKSpaceFlag=1
        AverageRawDataNoiseKSpaceFlag=0;
    end
    
    methods
        function this = cm2DRawDataReader(f)
            if nargin>0
                this.setFilename(f);
            end
        end
        function o=getFilename(this)
            o=this.FileName;
        end
        function setFilename(this,f)
            this.FileName=f;
        end
        
        function readFile(this)
            [pt,n,e]=fileparts(this.getFilename);
            switch lower(e)
                case '.h5'
                    this.Reader=CLOUDMRIsmrmRawDataReader(this.getFilename());
                case '.dat'
                    this.Reader=cm2DSiemensRawDataReader(this.getFilename());
            end
            
        end
        
        function setRawDataImageKSpace(this,f)
            this.RawDataImageKSpace=f;
        end
        
        function setRawDataNoiseKSpace(this,f)
            this.RawDataNoiseKSpace=f;
        end
        
        function setAverageRawDataImageKSpaceFlag(this,f)
            this.AverageRawDataImageKSpaceFlag=f;
        end
        
        function o=getAverageRawDataImageKSpaceFlag(this)
            o= this.AverageRawDataImageKSpaceFlag;
        end
        
                function setAverageRawDataNoiseKSpaceFlag(this,f)
            this.AverageRawDataImageKSpaceFlag=f;
        end
        
        function o=getAverageRawDataNoiseKSpaceFlag(this)
            o= this.AverageRawDataImageKSpaceFlag;
        end
        
        
     
        
        
        
        function o=getRawDataImageKSpace(this)
            
            if isempty(this.RawDataImageKSpace)
                            if (isempty(this.Reader))
                                this.readFile();
                            end            
                this.RawDataImageKSpace=this.Reader.readImageKSpace();
            end
            
            
            o=this.RawDataImageKSpace;
                   
            
            
        end
        
        
        
        function o=setcm2DRawDataReaderImageKSpaceFrom2DSlice(this,K,av,c,r)
            %it must be 4D freq,phase slice and coils
            [NF NP SL NC]=size(K);
           % nAvg,nContrasts,nReps,enc_Nx, enc_Ny, nSlices, nCoils
            NK=zeros(1,1,1,NF,NP,SL,NC);
            for s=1:SL
            NK(av,c,r,:,:,s,:)=K;
            end
            if isempty(this.RawDataImageKSpace)
            this.setRawDataImageKSpace(NK);
            end
        end
        
        
        function o=getRawDataNoiseKSpace(this)
            
            
           
            
            if isempty(this.RawDataNoiseKSpace)
                 if (isempty(this.Reader))
                this.readFile();
                end
                this.RawDataNoiseKSpace=this.Reader.readRawDataNoiseKSpace();
            end
            
            
            
            o=this.RawDataNoiseKSpace;
            
            
            
        end
        
        
        function O=getKSpaceDimensionsName(this)
            %display('1: Average','2: contrast','3: repetition,4: Frequency Encode','5: Phase Encode','6: Slice','7: Coils,');
            O={'1: Average','2: contrast','3: repetition','4: Frequency Encode','5: Phase Encode','6: Slice','7: Coils'};
        end
        
        function o=getKSpaceNoiseSlice(this,a,c,r,s)
            
          
            if (strcmp(a,'avg') || (a==0) || this.AverageRawDataNoiseKSpaceFlag())
               X=this.averageKS(this.getRawDataNoiseKSpace());
               a=0;
            else
                X=this.getRawDataNoiseKSpace();
            end
            
            o=this.reduceKSpace2DSlice(X,a,c,r,s);
            end
        
        
        
        
        
        
        
        function o=getRawDataRawDataImageKSpaceSlice(this,a,c,r,s)
            %this function get me the slice of the signal
            %average,contrast,repetition,slice)
            if ( this.getAverageRawDataImageKSpaceFlag())
               X=this.averageKS( this.getRawDataImageKSpace());
               a=0;
            else
                X= this.getRawDataImageKSpace();
            end
            
           
            
            
            
            o=this.reduceKSpace2DSlice(X,a,c,r,s);
        end
        
        function o=getNumberImageSlices(this)
            
            o=this.getNSlices(this.getRawDataImageKSpace());
        end
        
        
                function o=getNumberNoiseSlices(this)
            
            o=this.getNSlices(this.getRawDataNoiseKSpace());
        end
        
        
        function o=getNumberRepetition(this)
            
            o=this.getNRepetition(this.getRawDataImageKSpace());
        end
        
        
    
    
            function O= reduceKSpace2DSlice(this,K,a,c,r,s)
            %average,contrast,repetition,slice
            %if average =0 (we make average) else
            %we provide that particular average image the output is 3D (freq,fase,coils)
            
            S=cm2DRawDataReader.getKspace3Dsize(K);
            C=cm2DRawDataReader.getNCoil(K);
            %the average can be 1 if the kspace is averaged
            if ( this.getAverageRawDataImageKSpaceFlag())
                K=cm2DRawDataReader.averageKS(K);
                O=reshape(K(1,c,r,:,:,s,:),[S(1) S(2) C]);
            else
                O=reshape(K(a,c,r,:,:,s,:),[S(1) S(2) C]);
            end
            
            end
end
            
    
    methods (Static)
        
        
        function O=getNSlices(K)
            O= size(K,6);
        end
        
        
        function O=averageKS(K)
            O= mean(K,1);
        end
        
        function O=getNRepetition(K)
            O= size(K,3);
        end
        
        
        function O=getNAverages(K)
            O= size(K,1);
        end
        

                function O=getNContrast(K)
            O= size(K,2);
                end

        
                

        
        
        
        function O=getKspace3Dsize(K)
            
            O= [size(K,4) size(K,5) size(K,6)];
        end
        
        
        function O=getNCoil(K)
            
            O= size(K,7);
        end
        
    end
    
end

