
% READ ALL NOISE DATA FROM A SELECTED DIRECTORY AND GENERATES NOISE
% COVARIANCE AND NOISE COEFFICIENTS PLOTS FOR EACH OF THEM
%
% Riccardo Lattanzi: 30th October 2014

plot_single_cases = 0; % 1 --> plot noise covariance and noise coefficients on separate figures for all data files

which_slice = 1;
use_clims = 1; % 1 --> GENERATE THE FIGURE WITH ALL THE PLOTS WITH THE SAME COLORBAR
noisecovlims = [7.4121e-14 5.8958e-11];
noiseimglims = [0 0.8];
% NOTE: if use_clims = 1, the noise coefficients are plotted from 0 to 1
noise_BW = 0.792;

dirname = uigetdir; % select the directory
filePattern = fullfile(dirname, '*.dat'); % read all the DAT files inside
NoiseFiles = dir(filePattern);
disp(['The following ' num2str(length(NoiseFiles)) ' noise files will be analized'])
for ifile = 1:length(NoiseFiles)
    tempname = NoiseFiles(ifile);
    disp(['*** ' tempname.name ' ***']);
end

min_noisecov = noisecovlims(1);
max_noisecov = noisecovlims(2);
min_noiseimg = noiseimglims(1);
max_noiseimg = noiseimglims(2);

for ifile = 1:length(NoiseFiles)
    tempname = NoiseFiles(ifile);
    fullFileName = fullfile(dirname, tempname.name);
    disp('-----------------------------------')
    disp(['Current file: ' tempname.name])
    %     fprintf(1, 'Now reading %s\n', fullFileName);
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
        noiserawdata = permute(squeeze(kdata(:,:,:,:,which_slice,:,:,:,:,:,:,:,:,:,:,:)),[1,3,2]); % [nfreq nphase ncoil]
    else % single coil case
        noiserawdata = squeeze(kdata(:,:,:,:,which_slice,:,:,:,:,:,:,:,:,:,:,:)); % [nfreq nphase]
    end
    noiserawdata = noiserawdata*3200;
    nrow = size(noiserawdata,1);
    ncol = size(noiserawdata,2);
    
    noiseimg = MRifft(noiserawdata,[1,2]);
    % reconstruct individual coils' images and apply FFT scale factor
    % iFFT scales data by 1/N, so I need to scale back, as the noise covariance
    % matrix is calculated in k-space (not Fourier transformed)
    noiseimg = mean(noiseimg,3)*sqrt(nrow*ncol);
    if min_noiseimg < abs(min(noiseimg(:)))
        min_noiseimg = abs(min(noiseimg(:)));
    end
    if max_noiseimg > abs(max(noiseimg(:)))
        max_noiseimg = abs(max(noiseimg(:)));
    end
    if max_noiseimg < min_noiseimg
        max_noiseimg = min_noiseimg + 1;
    end
    allnoiseimg{ifile} = noiseimg;
    disp('***  start computing noise statistics...');
    [noisecov,noise_bandwidth] = calc_noise_cov(noiserawdata,0); % 0 --> without noise BW correction
    noisecov = noisecov/noise_BW;
    % check that clims are within actual range of values
    if min_noisecov < abs(min(noisecov(:)))
        min_noisecov = abs(min(noisecov(:)));
    end
    if max_noisecov > abs(max(noisecov(:)))
        max_noisecov = abs(max(noisecov(:)));
    end
    if max_noisecov < min_noisecov
        max_noisecov = min_noisecov + 1;
    end
    allnoisecov{ifile} = noisecov;
    noise_coeff = zeros(size(noisecov));
    for itemp = 1:size(noisecov,1)
        for jtemp = 1:size(noisecov,1)
            noise_coeff(itemp,jtemp) = noisecov(itemp,jtemp)/sqrt(noisecov(itemp,itemp)*noisecov(jtemp,jtemp));
        end
    end
    allnoisecoeff{ifile} = noise_coeff;
    disp('***  ...noise statistics computed');
    if plot_single_cases
        disp('***  ...plotting magnitude of noise covariance matrix');
        figure;
        imshow(abs(noisecov.'),[]); colormap(jet)
        colorbar
        title('NOISE COVARIANCE MATRIX');
        set(gcf,'name',tempname.name);
        disp('***  ...plotting magnitude of noise coefficient matrix');
        figure;
        imshow(abs(noise_coeff.'),[]); colormap(jet)
        colorbar
        title('NOISE COEFFICIENTS MATRIX');
        set(gcf,'name',tempname.name);
    end
end

% PLOT ALL TOGETHER
noisecovlims = [min_noisecov max_noisecov];
noiseimglims = [min_noiseimg max_noiseimg];
plot_rows = floor(length(NoiseFiles)/3);
if floor(length(NoiseFiles)/3) ~= ceil(length(NoiseFiles)/3)
    plot_rows = plot_rows + 1;
end
if use_clims
    figure;
    for iplot = 1:length(NoiseFiles)
        subplot(plot_rows,3,iplot);
        imshow(abs(cell2mat(allnoisecov(iplot))),noisecovlims); colormap(jet)
        tempname = NoiseFiles(iplot);
        tempname = tempname.name;
        MIDindex = strfind(tempname, 'MID');
        if isempty(MIDindex)
            title(tempname,'FontSize',14);
        else
            title(tempname(MIDindex:MIDindex+5),'FontSize',14);
        end
        set(gcf,'name',['NOISE COVARIANCE - Plot Limits = [' num2str(min_noisecov) ',' num2str(max_noisecov) ']'])
    end
    figure
    for iplot = 1:length(NoiseFiles)
        subplot(plot_rows,3,iplot);
        imshow(abs(cell2mat(allnoisecoeff(iplot))),[0 1]); colormap(jet)
        tempname = NoiseFiles(iplot);
        tempname = tempname.name;
        MIDindex = strfind(tempname, 'MID');
        if isempty(MIDindex)
            title(tempname,'FontSize',14);
        else
            title(tempname(MIDindex:MIDindex+5),'FontSize',14);
        end
        set(gcf,'name','NOISE COEFFICIENTS')
    end
    figure;
    for iplot = 1:length(NoiseFiles)
        subplot(plot_rows,3,iplot);
        imshow(abs(cell2mat(allnoiseimg(iplot))),noiseimglims); colormap(gray)
        tempname = NoiseFiles(iplot);
        tempname = tempname.name;
        MIDindex = strfind(tempname, 'MID');
        if isempty(MIDindex)
            title(tempname,'FontSize',14);
        else
            title(tempname(MIDindex:MIDindex+5),'FontSize',14);
        end
        set(gcf,'name',['NOISE IMAGES - Plot Limits = [' num2str(min_noiseimg) ',' num2str(max_noiseimg) ']'])
    end
else
    figure;
    for iplot = 1:length(NoiseFiles)
        subplot(plot_rows,3,iplot);
        imshow(abs(cell2mat(allnoisecov(iplot))),[]); colormap(jet)
        tempname = NoiseFiles(iplot);
        tempname = tempname.name;
        MIDindex = strfind(tempname, 'MID');
        if isempty(MIDindex)
            title(tempname,'FontSize',14);
        else
            title(tempname(MIDindex:MIDindex+5),'FontSize',14);
        end
        set(gcf,'name','NOISE COVARIANCE')
        colorbar
    end
    figure
    for iplot = 1:length(NoiseFiles)
        subplot(plot_rows,3,iplot);
        imshow(abs(cell2mat(allnoisecoeff(iplot))),[]); colormap(jet)
        tempname = NoiseFiles(iplot);
        tempname = tempname.name;
        MIDindex = strfind(tempname, 'MID');
        if isempty(MIDindex)
            title(tempname,'FontSize',14);
        else
            title(tempname(MIDindex:MIDindex+5),'FontSize',14);
        end
        set(gcf,'name','NOISE COEFFICIENTS')
        colorbar
    end
    figure;
    for iplot = 1:length(NoiseFiles)
        subplot(plot_rows,3,iplot);
        imshow(abs(cell2mat(allnoiseimg(iplot))),[]); colormap(gray)
        tempname = NoiseFiles(iplot);
        tempname = tempname.name;
        MIDindex = strfind(tempname, 'MID');
        if isempty(MIDindex)
            title(tempname,'FontSize',14);
        else
            title(tempname(MIDindex:MIDindex+5),'FontSize',14);
        end
        set(gcf,'name','NOISE IMAGES')
        colorbar
    end
end




