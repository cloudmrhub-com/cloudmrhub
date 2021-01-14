classdef cm2DSignalToNoiseRatioMultipleReplicas<cm2DSignalToNoiseRatio
    %main class of array combining methods, the constructor is ovewritten by
    %the class constructor
    
    properties
        imageArray
        MEAN
        STD
    end
    
    
    
    
    
    methods
        function this = cm2DSignalToNoiseRatioMultipleReplicas(x)
            %the class expects a 3D matrix composed by a tile of 2D images
            %or nothing
            this.Type='MR';
           
            if nargin>0
             this.add2DStackOfImages(x)                
            end
        end
        
        
        function add2DImage(this,IM)
            this.imageArray=cat(3,this.imageArray,  IM);
        end
        

         function add2DStackOfImages(this,x)
             for t=1:size(x,3)
                    this.add2DImage(x(:,:,t));
                end
        end

        
       function o=getImageArray(this)
            o=this.imageArray;
       end
        
        
        function setImageArray(this,set2dimages)
            
            this.imageArray=set2dimages;
            this.logIt('','ok')
        end

        
        function calculate(this)
                            
                this.MEAN=this.getImageArrayMEAN();
                this.STD=this.getImageArraySTD();
                this.SNR=this.MEAN./this.STD;
                        
        end
        
        
        function O=getSNR(this)
            if isempty(this.SNR)
                this.calculate();
            end
            
            O=this.SNR;
        end
        
        
        function O=getImageArrayMEAN(this)
            if isempty(this.MEAN)
                this.MEAN=nanmean(abs(this.getImageArray()),3);
            end
            
            O=this.MEAN;
        end
        
        
        function O=getImageArraySTD(this)
            if isempty(this.STD)
                this.STD=nanstd(abs(this.getImageArray()),0,3);
            end
            
            O=this.STD;
        end
        
        
        function plotImageArray(this,p)
            if nargin~=2
            p=0.5;
            end
            figure();
        im=this.getImageArray();for t=1:100;imshow(im(:,:,t),[]);title(['Replicas number:' num2str(t)]);pause(p); end
        end

        
        
    end

    
end



