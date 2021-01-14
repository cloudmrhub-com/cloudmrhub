function [KOUT, SHRINK,MASK, info]=undersampleGRAPPA2D(K,frequencyacceleration,phaseacceleration,frequencyautocalibration,phaseautocalibration)
%K is a kspace 2d data (frequency,phase,coil)
% phaseacceleration
% set ACS lines for mSense simulation (fully sampled central k-space
% region) % output image is 100% with 

[nX, nY, nCoils]=size(K);


%check autocalibration > size K
if (~exist('phaseautocalibration','var'))
        phaseautocalibration=nY;
 end
 
% 
if (~exist('frequencyautocalibration','var'))
    frequencyautocalibration=nX;
end


if (~exist('frequencyacceleration','var'))
    frequencyacceleration=1;
end
if (~exist('phaseacceleration','var'))
    phaseacceleration=1;
end



if phaseautocalibration>nY
    KOUT=NaN;
    return;
end

if frequencyautocalibration>nX
    KOUT=NaN;
    return;
end




[KOUT, SHRINK,MASK, info]=undersample2DData(K,frequencyacceleration,phaseacceleration,frequencyautocalibration,phaseautocalibration);
