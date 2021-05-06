classdef reconstructRawDataFilter<handle
    %need to have mroptimum code
    properties(Access=protected)
        SignalFilename %filename of the signal
        NoiseFilename %filename of the noise
        Reconstructor %the reconstruction you want
        NoiseReader
        SignalReader
    end
    
    methods
        function setSignalFilename(this,fn)
            this.SignalFilename=fn;
            
            A=reconstructRawDataFilter('signal');
            A.setFilename(this.getSignalFilename);              
            this.SignalReader=A.getOutput();
        end
        
        
                function o=getSignalFilename(this)
            o=this.SignalFilename;            
                end
        
                        function o=getNoiseFilename(this)
                                o=this.NoiselFilename;
                        end
        
        function setNoiseFilename(this,fn)
            this.NoiseFilename=fn;
            A=reconstructRawDataFilter('noise');
            A.setFilename(this.getSignalFilename);              
            this.NoiseReader=A.getOutput();   
        end
        
        function setReconstructor(this,recon)
            this.Reconstructor=recons;
        end
        
        function getOutput(this)
            
        end
    end
    
    
    methods(Static)
        function jop=getJop(recon)
            jop=mro2DReconGetDefaultOptionsForType(recon);
        end
    end
    
end