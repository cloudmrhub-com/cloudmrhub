classdef cm2DReconGRAPPA<cm2DReconWithSensitivityAutocalibrated
    %this is an initial class derived by the code developped at cbi.
    %sensitivity.
    %http://mriquestions.com/senseasset.html
    %last update 2020 Dec
    %update 2020 Jan
    %---------------------------------------------------------------------
    %kspace is at the fully sampled size and with zeros on the accelerated
    %positions
    
    properties(Access=private)
        GrappaKernel ; % must be in [odd even] numbers format
        SosRecon=0; %default is adapt
      end
    
    methods
        
             function this=cm2DReconGRAPPA()
            this.setAccelerationFrequency (1);
            this.setHasSensitivity(0);
        end        
        function o=getAutocalibration(this)
            o=[this.getAutocalibrationFrequency() this.getAutocalibrationPhase()];
        end
        
        
                function setAcceleration(this,a)
                   this.setAccelerationPhase(a)
                   this.setAccelerationFrequency(1)
                end
                
                                function o=getAcceleration(this)
                   o=this.getsetAccelerationPhase();
                                end

                
                function setAutocalibration(this,acl)
                %acl array length 2
            this.setAutocalibrationFrequency(acl(1)) 
            this.setAutocalibrationPhase(acl(2));
        end
        

        function o=getGrappaKernel(this)
            
            o=this.GrappaKernel;
        end
        
        function o=setGrappaKernel(this,GK)
            
                %acl array length 2
                
                if(isequal(mod(GK,2),[1 0])) %firt must be odd the second even
            this.GrappaKernel=GK;
                else
                    display('Grappa Kernel must be [odd, even]')
                    this.errorMessage()
                end
        end

        
             function o=getSosRecon(this)
            o=this.SosRecon;
        end
        
        function setSosRecon(this,SR)
                
            this.SosRecon=SR;
        end
        
        
        %override
      function o=getAccelerationFrequency(this)
          this.logIt('GRAPPA 2D can be accelerated only in the phase dimension','warning')  
          o=1;
        end
        %override
        function setAccelerationFrequency(this,R1)
            this.logIt('GRAPPA 2D can be accelerated only in the phase dimension','warning')  
            
        end
        
        
        function R=getR(this)
            this.logIt('get R for grappa','?');
            try
            SS=this.getSignalSize();
             R = this.getAccelerationPhase();
             np=SS(2);
             
    if mod(np,R) ~= 0
        tempR = R;
        while (mod(np,tempR) ~= 0)
            tempR = tempR-1;
        end
        R = tempR;
        ss=['Acceleration R reduced to ' num2str(R) ', so the number of lines can be exactly divided by R'];
        disp(ss);
        this.logIt(ss,'ok');
    end
            catch
                this.logIt('can not get the R for grappa','error');
            end
        end

                    
        
        
       function ACS=getSignalPrewhitenedAutocalibrationsArea(this)
           ACS=this.getAutocalibrationsLinesKSpace(this.getPrewhitenedSignal(),this.getAutocalibrationFrequency, this.getAutocalibrationPhase);
         end
        
        function MRimage=getOutput(this)
            %Riccardo Lattanzi riccardo.lattanzi@nyulangone.org
            %GRAPPA_img_recon December 2020
            
              [SS]=this.getSignalSize();
             nf=SS(1);
             np=SS(2);
             nc = this.getSignalNCoils();
             R=this.getR();
             grappa_kernel=this.getGrappaKernel();
             data_acs=this.getSignalPrewhitenedAutocalibrationsArea();
            pw_signalrawdata=this.getPrewhitenedSignal(); 
            % decimate data
            k_temp=this.getUndersampledKSpace(pw_signalrawdata,1,R); % nf x np/R x nc

             Rn=eye(nc);
             sos_coil_combination=this.getSosRecon(); %sos or adapt 0 or 1
             fullKS=this.getGRAPPAKspace(k_temp,data_acs,grappa_kernel,R);
            
  
             %recon     
             if sos_coil_combination
    % sum of squares 
    recon=cm2DReconRSS();
else
    % adaptive combine reconstruction
    recon=cm2DReconAdapt();
    %MRimage= adapt_array_2d(img_matrix,Rn,1);
             end
    recon.setPrewhitenedSignal(fullKS);
    MRimage = recon.getOutput();
    MRimage=this.demimicAccelerated2DImage(MRimage);
    %T=GRAPPA_img_recon(k_temp,eye(nc),data_acs,grappa_kernel,R,sos_coil_combination);
           
             
        end
        
        
 
        
        
        
        
    end
      methods (Static)
            
          function krecon=getGRAPPAKspace(rawdata,acs,ksize,R)
              
 ksf=ksize(1); % MUST BE AN ODD NUMBER
 ksp=ksize(2); % MUST BE AN EVEN NUMBER

[nf,np_acc,nc] = size(rawdata);
np = np_acc*R;
[nfacs,npacs,nc]=size(acs);

%% CALIBRATION STEP %%
% calculate GRAPPA reconstruction weights %

% initialize source points (Nc*Nkernelpoints X Nrepetition)
src = zeros(nc*ksp*ksf,(npacs-(ksp-1)*R)*(nfacs-(ksf-1)));
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% There ksf*ksp*nc points for each repetition of the kernel within the acs 
% lines. (ksp-1) and (ksf-1) because the kernel has to stay within the acs 
% lines, multiplied by R because when the kernel hit the edge, the space 
% between kernel lines include the skipped kspace lines
%
% Example for R = 2 along phase and kernel size 5x4 (ksfxksp):
%
% -O-O-O-O-O- 1
% - - - - - - 2
% -O-O-O-O-O- 3
% - - -X- - - 4
% -O-O-O-O-O- 5
% - - - - - - 6
% -O-O-O-O-O- 7
% - - - - - - 8
% - - - - - - 9
% - - - - - - 10
% ... etc ...
%
% The circles are the source points, and the X are the target points.
% When the kernel slides downward and the current line 7 hits the bottom of
% the acs lines, the (ksp-1)*R = 3*2 = 6 lines will not be used to collect
% additional source points
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% initialize target points (Nrepetitions X Nc*Nmissingkspacelines)
trg = zeros(nc*(R-1),(npacs-(ksp-1)*R)*(nfacs-(ksf-1)));
% All source points in all coils are used to fit one target point in each individual coil

dx = (ksf-1)/2;
dy = (ksp/2-1)*R;
src_count = 0;                         
% reordering data
for irow = dx+1:nfacs-dx % target point is at the center of the kernel
    for icol = 1:npacs-(ksp-1)*R % span all the acs lines 
        src_count = src_count + 1;
        % source points (nc*ksp*ksf)
        src(:,src_count) = reshape(acs(irow-dx:irow+dx,icol:R:icol+(ksp-1)*R,:),[nc*ksp*ksf 1]);
        % target points (nc*(R-1))
        trg(:,src_count) = reshape(acs(irow,icol+1+dy:icol+dy+R-1,:),[nc*(R-1) 1]);
    end
end
% Find GRAPPA weights by solving w*S_source = S_target
% weights (nc*(R-1) X nc*ksp*ksf)
ws = (trg)*pinv(src);

%% SYNTHESIS STEP %%
% reconstruct missing k-space points

% The borders of k-space are zero-filled to accomodate the GRAPPA convolution operator
krecon = zeros(nf+2*dx,np+2*dy+1,nc); 

% Available undersampled data are copied into zero-padded matrix

krecon(dx+1:end-dx,dy+1:R:np+dy,:) = rawdata; 

for irow = dx+1:nf+dx
    for icol= 1:R:np
        src = reshape(krecon(irow-dx:irow+dx,icol:R:icol+(ksp-1)*R,:),[nc*ksf*ksp 1]);
        % Apply weights to the source points
        krecon(irow,icol+dy+1:icol+dy+R-1,:)=reshape(ws*src,[(R-1) nc]);
    end
end
% remove the extra borders of k-spae
krecon = krecon(dx+1:nf+dx,dy+1:np+dy,:);

%% IMAGE RECONSTRUCTION STEP %%
% reconstruct coil images and combine them
% 
% dim = [1,2];
% tmp = ifftshift(krecon);
% 
% for idim = dim
%     tmp = ifft(tmp,[],idim);
% end
% % Generate aliased multi-channel images
% img_matrix = fftshift(tmp);
% 
% % Since iFFT scales data by 1/sqrt(N), we need to scale the image back by the same quantity
% img_matrix = img_matrix*sqrt(nf*np);
% 
% if 0
% % Since only nf*np/Rtot samples are actually sampled, we need to scale the image by the square root of the
% % total acceleration factor to maintain the same scaling as for the noise (fully sampled)
% img_matrix = img_matrix*sqrt(Rtot);
%end
          end
          
            function TEST=therecon()
                    
                %instantiate the class
                TEST=cm2DReconGRAPPA();
                [K,nc]=TEST.getMRoptimumTestData();
                
                R1=2; % no for grappa 2d
                R2=4;
                ACL1=4;
                ACL2=2;
                GK=[5 4];
                KS=undersampleGRAPPA2D(K,R1,R2,ACL1,ACL2);
                TEST.setSignalKSpace(KS);
                TEST.setNoiseCovariance(nc);
                %please set noise before set the sens
                %TEST.setCoilSensitivityMatrix(KS)
                
                TEST.setAccelerationPhase(R2);
                TEST.setAutocalibration([ACL1 ACL2])
                TEST.setGrappaKernel(GK);
                end
        end
end


