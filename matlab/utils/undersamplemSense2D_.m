[nX, nY, nCoils]=size(K);


%check autocalibration > size K
if (~exist('autocalibration','var'))
    autocalibration=NaN;
end


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


[KOUT, SHRINK,MASK, info]=undersample2DData(K,frequencyacceleration,phaseacceleration,NaN,autocalibration);