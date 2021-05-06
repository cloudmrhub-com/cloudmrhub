classdef readRawDataFilter<handle
    
    properties(Access=protected)
        Filename %filename of the signal
        isSignal
    end
    properties
        Reader
    end
    
    methods
        function this=readRawDataFilter(n)
            %declare if it's noise or signal
            switch(n)
                case {'noise'}
                    this.isSignal=0;
                case {'signal','data'}
                    this.isSignal=1;
            end
            
        end
        
        function setFilename(this,fn)
            
            this.Filename=fn;
            this.Reader=cm2DRawDataReader();
            this.Reader.setIsSignalFile(this.isSignal);
            this.Reader.setFilename(this.getFilename());              
        end

                function o=getFilename(this)
                o=this.Filename;
                end
        function o=getOutput(this)
            o=this.Reader;
        end
        
        
        function o=getFileinfo(this)
            J=this.Reader;
            o=J.getFileinfo();
        end
        
        
        

    end
    
end