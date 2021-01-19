classdef cm2DSignalToNoiseRatioPseudoMultipleReplicasWien<cm2DSignalToNoiseRatioPseudoMultipleReplicas
    %main class of array combining methods, the constructor is ovewritten by
    %the class constructor
    
    properties
    BoxSize=3;
            
    end
    
    
    
    
    
    methods
        function this = cm2DSignalToNoiseRatioPseudoMultipleReplicasWien(Recon,x)
            %the class expects a 3D matrix composed by a tile of 2D images
            %or nothing
            this.Type='CR';
            
                if nargin>0
                 this.setReconstructor(Recon)
            end
            
            if nargin>1
                this.add2DStackOfImages(x)
            end
            
        end
        
             
        function readConf(this,js)
            
            this.NumberOfPseudoReplicas=js.NR;
            try this.BoxSize=js.boxsize;end;
        end
        
        
        
        function O=getParams(this)
            
            O.Type=this.Type;
            O.NR=this.NumberOfPseudoReplicas;
            O.boxsize=this.BoxSize;
            
        end
        
        
        
        
        
        function o=getBoxSize(this)
            o=this.BoxSize;
        end
        

         function o=setBoxSize(this,b)
            this.BoxSize=b;
        end

        function noise=getPseudoMultipleReplicaNoiseImage(this)
            ref=this.getPseudoMultipleReplicaSignalImage();
            M= this.getImageArrayMEAN();
            noise= this.getWienNoiseImage(ref-M,this.BoxSize);
        end
        
        
        
            
        
        
        
       
        
                function o=getPseudoMultipleReplicaSignalImage(this)
            o= this.Reconstructor.getOutput();
        end
        
        
        
    end
    
    
    methods (Static)
         function NOISE=getWienNoiseImage(noiseonly,box)
            %noiseonly is the reconstructed noiseonly image 2D 1 slice, box the size
            %of the box
            
            [NC,NR]=size(noiseonly);
            
            kx=box;
            ky=box; 
            NOISE=nan(NC,NR);
        


             PADDEDNOISE = padarray(noiseonly,[kx ky],NaN);
                for ic=1:NC
                    for ir=1:NR
                        pic=kx+ic+[-kx:kx];
                        pir=ky+ir+[-ky:ky];
                        try
                            NOISE(ic,ir)=nanstd(reshape(PADDEDNOISE(pic,pir),[],1));
                        catch
                            
                        end
                    end
                end
                
        end
    end
    
    

    
end



