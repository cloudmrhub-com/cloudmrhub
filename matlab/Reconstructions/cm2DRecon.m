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
            %is a third element is given that's the prewhitened signal
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
    
end



