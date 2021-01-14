function [out] = testSNR___(o,signal,noise,sensitivities)


ACMs=CLOUDMRgetclassfromOptions(o);

ACMs.setConf(o);
ACMs.setSignalKSpace(signal);
 
ACMs.setNoiseKSpace(noise);

if(exist('sensitivities','var'))
   ACMs.setSourceCoilSensitivityMap(sensitivities)
end

out=ACMs.getSNR();

end

