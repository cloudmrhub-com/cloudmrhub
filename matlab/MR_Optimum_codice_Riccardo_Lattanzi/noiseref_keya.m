data_file = [];
which_slice = 1;           % select which slice to use in rawdata with multiple slices
K_ICE_AMPL_SCALE_FACTOR = 3200; % ICE multiplies data by this quantity

if isempty(data_file)
    [data_file,data_path] = uigetfile('*.dat','Select the data raw file');
    disp('***  loading MR signal from:');
    disp(['   ' fullfile(data_path,data_file)]);
    disp('---');
    data_file = fullfile(data_path,data_file);
end

twix_obj = mapVBVD(data_file);

noisedata = twix_obj.noise.unsorted;
kdata = twix_obj.image();
clear image_obj
totalslices = size(kdata,5);
if which_slice > totalslices,
    which_slice = 1;
    disp('*** WARNING: whichslice was out of range so it was set to 1');
end
signalrawdata = permute(squeeze(kdata(:,:,:,:,which_slice,:,:,:,:,:,:,:,:,:,:,:)),[1,3,2]); % [nfreq nphase ncoil]
signalrawdata = signalrawdata*K_ICE_AMPL_SCALE_FACTOR;

noiserawdata = permute(squeeze(noisedata(:,:,:,:,which_slice,:,:,:,:,:,:,:,:,:,:,:)),[1,3,2]); % [nfreq nphase ncoil]
noiserawdata = noiserawdata*K_ICE_AMPL_SCALE_FACTOR;

clear noise_obj
clear kdata