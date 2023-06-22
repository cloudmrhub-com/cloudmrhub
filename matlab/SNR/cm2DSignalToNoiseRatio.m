classdef cm2DSignalToNoiseRatio<cmOutput
    %main class of array combining methods, the constructor is ovewritten by
    %the class constructor
    properties
        imageSize
        SNR
        Type
        SubType
    end
    methods
        function this = cm2DSignalToNoiseRatio()
            %the class expects a 3D matrix composed by a tile of 2D images
            %or nothing
            this.logIt('SNR Calculation started instantiated','start')
           
        end
         function O=getOutput(this)
             O=this.getSNR();
         end
    end
end



