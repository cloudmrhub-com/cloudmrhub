function [KOUT, SHRINK,MASK, info]=undersampleSense2D(K,frequencyacceleration,phaseacceleration)
%K is a kspace 2d data (frequency,phase,coil)
% phaseacceleration
% set ACS lines for mSense simulation (fully sampled central k-space
% region)
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

[KOUT, SHRINK,MASK, info]=undersample2DData(K,frequencyacceleration,phaseacceleration,NaN,NaN);