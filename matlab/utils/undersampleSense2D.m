function [KOUT, SHRINK,MASK, info]=undersampleSense2D(K,frequencyacceleration,phaseacceleration)
%K is a kspace 2d data (frequency,phase,coil)
% phaseacceleration
% set ACS lines for mSense simulation (fully sampled central k-space
% region)




[nX, nY, nCoils]=size(K);


%check autocalibration > size K

    autocalibration=0;



if (~exist('frequencyacceleration','var'))
    frequencyacceleration=1;
end
if (~exist('phaseacceleration','var'))
    phaseacceleration=1;
end



if autocalibration>=nY
    KOUT=NaN;
    return;
end


if mod(autocalibration,2)>0
    KOUT=NaN;
    return;
end

ACShw = autocalibration/2 ; 

Ysamp_u = [1:phaseacceleration:nY] ; % undersampling 
Ysamp_ACS = [floor(nY/2-ACShw+1) : floor(nY/2+ACShw)] ; % GRAPPA autocalibration lines
Ysamp = union(Ysamp_u, Ysamp_ACS) ; % actually sampled lines
nYsamp = length(Ysamp) ; % number of actually sampled

Xsamp_u = [1:frequencyacceleration:nX] ; % undersampling 
Xsamp_ACS = [floor(nX/2-ACShw+1) : floor(nX/2+ACShw)] ; % GRAPPA autocalibration lines
Xsamp = union(Xsamp_u, Xsamp_ACS) ; % actually sampled lines
nXsamp = length(Xsamp) ; % number of actually sampled





% Ysamp indexes the actually sampled lines to the encoded k-space line number. 
% For example, if there were just regular factor 2 undersampling 
% (with no ACS lines), Ysamp would have length 128 and be [1 3 5 ... 255].
% With ACS lines, the elements of Ysamp are separated by 2 near the k-space
% edges, and by 1 in the central ACS region.


KOUT = zeros(size(K));
KOUT(Xsamp,Ysamp,:) = K(Xsamp,Ysamp,:);

if nargout>1
    SHRINK=zeros(nXsamp,nYsamp,nCoils);
    for t=1:nXsamp
        for y=1:nYsamp
            SHRINK(t,y,:)=K(Xsamp(t),Ysamp(y),:);
        end
    end
    
end

MASK= zeros(size(K));


MASK(Xsamp,Ysamp,:)=1;
info.NX=nXsamp;
info.NY=nYsamp;
info.NC=nCoils;


