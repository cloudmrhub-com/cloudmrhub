function [O]=spread2DRawDataKSpaceNoiseInChannels(NRD,noisefile)
%cloudmrrd rd file, set noisefile to true  if kspace of the noise is on a
%separate file, false if it is on the signal kspace (multiraid)

if noisefile
    K=NRD.getRawDataImageKSpace();
else
    K=NRD.getRawDataNoiseKSpace();
end

O=[];

NCLOILS=NRD.getNCoil(K);
NSLICES=NRD.getNSlices(K);
NREP=NRD.getNRepetition(K);
NAV=NRD.getNAverages(K);
NCNT=NRD.getNContrast(K);

for na=1:NAV
    for nco=1:NCNT
        for nrep=1:NREP
            for nsl=1:NSLICES
                O=cat(2,O,squeeze(K(na,nco,nrep,:,:,nsl,:)));
            end
        end
    end

end
