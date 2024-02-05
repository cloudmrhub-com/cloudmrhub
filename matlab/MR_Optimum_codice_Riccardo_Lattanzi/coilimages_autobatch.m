
% READ ALL SIGNAL DATA FROM A SELECTED DIRECTORY AND PLOTS COIL IMAGES FOR
% EACH OF THEM
%
% Riccardo Lattanzi: 30th October 2014

which_slice = 1;
signallims = [1.2 1.8];
fixcaxisrange = 0;
min_signal = signallims(1);
max_signal = signallims(2);

[filename, pathname] = uigetfile('*.dat','Pick a signal file');
fullFileName = fullfile(pathname, filename);
disp('-----------------------------------')
disp(['Reading file: ' fullFileName])
disp('...')
image_obj = mapVBVD(fullFileName,'doAverage');
kdata = image_obj.image();
clear image_obj
totalslices = size(kdata,5);
if which_slice > totalslices,
    which_slice = 1;
    disp('*** WARNING: whichslice was out of range so it was set to 1');
end
if size(kdata,2) > 1,
    signalrawdata = permute(squeeze(kdata(:,:,:,:,which_slice,:,:,:,:,:,:,:,:,:,:,:)),[1,3,2]); % [nfreq nphase ncoil]
else % single coil case
    signalrawdata = squeeze(kdata(:,:,:,:,which_slice,:,:,:,:,:,:,:,:,:,:,:)); % [nfreq nphase]
end
signalrawdata = signalrawdata*3200;
nrow = size(signalrawdata,1);
ncol = size(signalrawdata,2);
nchan = size(signalrawdata,3);

disp('***  start reconstructing individual coil images...');
img_matrix = MRifft(signalrawdata,[1,2]);
img_matrix = img_matrix*sqrt(nrow*ncol);

for icoil = 1:nchan
    currentcoilimg = img_matrix(:,:,icoil);
    disp(['coil #' num2str(icoil) ': min = ' num2str(abs(min(currentcoilimg(:)))) ', max = ' num2str(abs(max(currentcoilimg(:))))]);
    if fixcaxisrange
        % NOTE: need this check to avoid setting caxis limits outside a coil range
        if min_signal < abs(min(currentcoilimg(:)))
            min_signal = abs(min(currentcoilimg(:)));
        end
        if max_signal > abs(max(currentcoilimg(:)))
            max_signal = abs(max(currentcoilimg(:)));
        end
        if max_signal < min_signal
            max_signal = min_signal + 1;
        end
    end
end
signallims = [min_signal max_signal];
if (nchan/4 < 1)
    figure;
    set(gcf,'name',['File: ' filename 'Plot Limits = [' num2str(min_signal) ',' num2str(max_signal) ']']);
    for ichan = 1:nchan
        subplot(1,nchan,ichan);
        imshow(abs(img_matrix(:,:,ichan)),signallims);
        title(['Chan ', num2str(ichan)]);
    end
    colormap(used_colormap);
else
    figure;
    set(gcf,'name',['File: ' filename 'Plot Limits = [' num2str(min_signal) ',' num2str(max_signal) ']']);
    for ichan = 1:nchan
        if mod(nchan,4) == 0
            subplot(4,nchan/4,ichan);
        else
            subplot(4,floor(nchan/4)+1,ichan);
        end
        imshow(abs(img_matrix(:,:,ichan)),signallims);
        title(['Chan ', num2str(ichan)]);
        
    end
    colormap(jet);
end
