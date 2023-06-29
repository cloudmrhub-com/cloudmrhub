function [KOUT, SHRINK,MASK, info,ACL]=undersample2DData(K,frequencyacceleration,phaseacceleration,frequencyautocalibration,phaseautocalibration)
%K is a kspace 2d data (frequency,phase,coil)
% frequency accelration, phaseacceleration,
% set ACS lines for mSense simulation (fully sampled central k-space
% region) % output image is 100% with 




[nX, nY, nCoils]=size(K);

if frequencyacceleration>=nX
    KOUT=NaN;
    return;
end

if phaseacceleration>=nY
    KOUT=NaN;
    return;
end


ACL=[];
ACShw = phaseautocalibration/2 ; 
ACShh = frequencyautocalibration/2 ; 

Ysamp_u = [1:phaseacceleration:nY] ; % undersampling
if (~isnan(phaseautocalibration))
    if(phaseautocalibration>0)
        Ysamp_ACS = [floor(nY/2)-ACShw+1 : floor(nY/2)+ACShw] ; % GRAPPA autocalibration lines
        Ysamp = union(Ysamp_u, Ysamp_ACS) ; % actually sampled lines
    else
        Ysamp_ACS=[];
        Ysamp=Ysamp_u;
    end
else
    Ysamp_ACS=[];
    Ysamp=Ysamp_u;
end


nYsamp = length(Ysamp) ; % number of actually sampled
nYsamp_u=length(Ysamp_u);

Xsamp_u = [1:frequencyacceleration:nX] ; % undersampling
if (~isnan(frequencyautocalibration))
    
    if (frequencyautocalibration>0)
        
        Xsamp_ACS = [floor(nX/2)-ACShh+1 : floor(nX/2)+ACShh] ; % GRAPPA autocalibration lines
        Xsamp = union(Xsamp_u, Xsamp_ACS) ; % actually sampled lines
    else
        Xsamp_ACS=[];
        Xsamp=Xsamp_u;
    end
    
else
    Xsamp_ACS=[];
    Xsamp=Xsamp_u;
    
end

nXsamp = length(Xsamp) ; % number of actually sampled
nXsamp_u=length(Xsamp_u);




% Ysamp indexes the actually sampled lines to the encoded k-space line number. 
% For example, if there were just regular factor 2 undersampling 
% (with no ACS lines), Ysamp would have length 128 and be [1 3 5 ... 255].
% With ACS lines, the elements of Ysamp are separated by 2 near the k-space
% edges, and by 1 in the central ACS region.


KOUT =zeros(size(K));
KOUT(Xsamp,Ysamp,:) = K(Xsamp,Ysamp,:);

if nargout>1
    SHRINK=zeros(nXsamp_u,nYsamp_u,nCoils);
    for t=1:nXsamp_u
        for y=1:nYsamp_u
            SHRINK(t,y,:)=K(Xsamp_u(t),Ysamp_u(y),:);
        end
    end
    
end

MASK= zeros(size(K));


MASK(Xsamp,Ysamp,:)=1;
info.NX=nXsamp;
info.NY=nYsamp;
info.NC=nCoils;
info.ACL.X=Xsamp_ACS;
info.ACL.Y=Ysamp_ACS;





   
if nargout>4
    ACL=K(Xsamp_ACS,Ysamp_ACS,:);
end

        