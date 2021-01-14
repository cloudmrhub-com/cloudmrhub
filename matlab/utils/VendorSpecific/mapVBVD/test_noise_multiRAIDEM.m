
singleraidsignalfile = '/home/montie01/Downloads/Signal_Data_Multislice.dat';
singleraidnoisefile = '/home/montie01/Downloads/Noise_Data_Multislice.dat';
multiraidfile = '/home/montie01/Downloads/Signal_Data_Multislice_MULTIRAID.dat';

% singleraidnoisefile = '../data/test_multiRAID_data/meas_MID00027_FID188181_Multislice_no_RF_single.dat';
% multiraidfile = '../data/test_multiRAID_data/meas_MID00024_FID188178_Multislice.dat';


plot_coil_sensitivities = 1;

% READ noise-only single RAID file
image_obj = mapVBVD(singleraidnoisefile,'removeOS');
kdata = image_obj.image();
singleR_noise_dwelltime = image_obj.hdr.Meas.RealDwellTime(2);
disp(['Dwell time single RAID noise = ' num2str(singleR_noise_dwelltime)])
nchan = size(kdata,2);
nslice = size(kdata,5);
singleR_noisedata = squeeze(permute(kdata,[1 3 2 5 4]));
% singleR_noisedata = singleR_noisedata(:,:,:,1); % Slice #1


%multislice
sl=1;
NOISE=squeeze(singleR_noisedata(:,:,:,sl));
for s =2:nslice
    NOISE=cat(2,NOISE, squeeze(singleR_noisedata(:,:,:,sl)));
end

[singleR_noisecov,singleR_noise_bandwidth] = calc_noise_cov(NOISE,1);
clear image_obj kdata

% READ noise prescan from multi RAID file
image_obj = mapVBVD(multiraidfile,'removeOS');
if (length(image_obj) > 1)
    if (length(fieldnames(image_obj{end})) == 0)
        image_obj = image_obj{1};
        prescan_obj = [];
    else
        prescan_obj = image_obj{end-1};
        image_obj = image_obj{end};
    end
else
    prescan_obj = [];
end

kdata_dwelltime = image_obj.hdr.Meas.RealDwellTime(2);
disp(['Dwell time multi RAID signal = ' num2str(kdata_dwelltime)])
% knoise_dwelltime = prescan_obj.hdr.Meas.RealDwellTime(2);
knoise_dwelltime = 5000;
correction_factor = sqrt(knoise_dwelltime/kdata_dwelltime);

% READ signal from multi RAID file
kdata = image_obj.image();
multiR_signaldata = squeeze(permute(kdata,[1 3 2 5 4]));
multiR_signaldata = multiR_signaldata(:,:,:,1);
clear image_obj kdata

multiR_nrow = size(multiR_signaldata,1);
multiR_ncol = size(multiR_signaldata,2);
disp(['Multi RAID signal dims: ' num2str(multiR_nrow) ' x ' num2str(multiR_ncol)])

knoise = prescan_obj.noise();

%qui entrano in gioco gli averages


display(['multiraid noise has ' num2str(size(knoise,6)) 'averages'])

%tu qui togli gli averages
% knoise = knoise(:,:,:,1,1,1);

multiraidNoise=permute(knoise,[1 3 2,4,5,6]);

ave=1;
NOISE=squeeze(multiraidNoise(:,:,:,1,1,ave));
for s =2:size(multiraidNoise,6)
    NOISE=cat(2,NOISE, squeeze(multiraidNoise(:,:,:,1,1,ave)));
end


knoise=NOISE;

knoise = knoise*correction_factor;

multiR_noisedata = knoise;

[multiR_noisecov,multiR_noise_bandwidth] = calc_noise_cov(multiR_noisedata,1);

disp('------------------------')
disp(['    Single RAID noise BW = ' num2str(singleR_noise_bandwidth)])
disp('------------------------')
disp(['    Multi RAID noise BW = ' num2str(multiR_noise_bandwidth)])
disp('------------------------')
disp(['    Single RAID noise std = ' num2str(std(real(singleR_noisedata(:))))])
disp('------------------------')
disp(['    Multi RAID noise std = ' num2str(std(real(multiR_noisedata(:))))])

figure;
subplot(3,2,1)
imagesc(abs(singleR_noisecov)); colorbar; axis square; title('Single RAID')
subplot(3,2,2)
imagesc(abs(multiR_noisecov)); colorbar; axis square; title('Multi RAID')
subplot(3,2,3)
imagesc(abs(singleR_noisecov./multiR_noisecov)); colorbar; axis square; title('Abs(Single/Multi) RAID')
subplot(3,2,4)
imagesc(abs(singleR_noisecov-multiR_noisecov)); colorbar; axis square; title('Abs(Single - Multi) RAID')
subplot(3,2,5)
hist(real(singleR_noisedata(:))); colorbar
subplot(3,2,6)
hist(real(multiR_noisedata(:))); colorbar

% ---------------------------
% ESTIMATE SNR FOR COMPARISON
% ---------------------------

% READ signal single RAID file
image_obj = mapVBVD(singleraidsignalfile);
kdata = image_obj.image();
singleR_dwelltime = image_obj.hdr.Meas.RealDwellTime(2);
disp(['Dwell time single RAID signal = ' num2str(singleR_dwelltime)])
singleR_signaldata = squeeze(permute(kdata,[1 3 2 5 4]));
singleR_signaldata = singleR_signaldata(:,:,:,1); % slice #1
clear image_obj kdata

singleR_nrow = size(singleR_signaldata,1);
singleR_ncol = size(singleR_signaldata,2);
disp(['Single RAID signal dims: ' num2str(singleR_nrow) ' x ' num2str(singleR_ncol)])

% ---------------
% SINGLE RAID SNR
% ---------------
singleR_img_matrix = MRifft(singleR_signaldata,[1,2]);
% reconstruct individual coils' images and apply FFT scale factor
% iFFT scales data by 1/N, so I need to scale back, as the noise covariance
% matrix is calculated in k-space (not Fourier transformed)
singleR_img_matrix = singleR_img_matrix*sqrt(singleR_nrow*singleR_ncol);
reference_image = sqrt(sum(abs(singleR_img_matrix).^2,3));
singleR_coilsens_set = singleR_img_matrix./repmat(reference_image,[1 1 nchan]);
disp('***  ...coil sensitivities computed...');

if plot_coil_sensitivities
    if (nchan/4 < 1)
        figure;
        set(gcf,'name','Coil Sensitivities Multi RAID');
        for ichan = 1:nchan
            subplot(1,nchan,ichan);
            %                             % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
            %                             imshow(abs(coilsens_set((size(coilsens_set,1)/4+1):end-size(coilsens_set,1)/4,:,ichan)),[]);
            %                             title(['Chan ', num2str(ch_id_track(ichan))]);
            imshow(abs(singleR_coilsens_set(:,:,ichan)),[]);
            title(['Chan ', num2str(ichan)]);
        end
        colormap(used_colormap);
    else
        figure;
        set(gcf,'name','Coil Sensitivities Multi RAID');
        for ichan = 1:nchan
            if mod(nchan,4) == 0
                subplot(4,nchan/4,ichan);
            else
                subplot(4,floor(nchan/4)+1,ichan);
            end
            %                             % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
            %                             imshow(abs(singleR_coilsens_set((size(singleR_coilsens_set,1)/4+1):end-size(singleR_coilsens_set,1)/4,:,ichan)),[]);
            %                             title(['Chan ', num2str(ch_id_track(ichan))]);
            imshow(abs(singleR_coilsens_set(:,:,ichan)),[]);
            title(['Chan ', num2str(ichan)]);
        end
        colormap(parula);
    end
end
singleR_snr_opt = zeros(singleR_nrow,singleR_ncol);
for irow = 1:singleR_nrow
    for icol = 1:singleR_ncol
        s_matrix = squeeze(singleR_coilsens_set(irow,icol,:));
        singleR_snr_opt(irow,icol) = sqrt(2)*abs((s_matrix')*inv(singleR_noisecov)*squeeze(singleR_img_matrix(irow,icol,:)))/...
            sqrt((s_matrix')*inv(singleR_noisecov)*s_matrix);
    end
end
singleR_snr_map = abs(singleR_snr_opt);
singleR_snr_map = fliplr(singleR_snr_map);


% ---------------
% MULTI RAID SNR
% ---------------
multiR_img_matrix = MRifft(multiR_signaldata,[1,2]);
% reconstruct individual coils' images and apply FFT scale factor
% iFFT scales data by 1/N, so I need to scale back, as the noise covariance
% matrix is calculated in k-space (not Fourier transformed)
multiR_img_matrix = multiR_img_matrix*sqrt(multiR_nrow*multiR_ncol);
reference_image = sqrt(sum(abs(multiR_img_matrix).^2,3));
multiR_coilsens_set = multiR_img_matrix./repmat(reference_image,[1 1 nchan]);
disp('***  ...coil sensitivities computed...');

if plot_coil_sensitivities
    if (nchan/4 < 1)
        figure;
        set(gcf,'name','Coil Sensitivities Single RAID');
        for ichan = 1:nchan
            subplot(1,nchan,ichan);
            %                             % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
            %                             imshow(abs(coilsens_set((size(coilsens_set,1)/4+1):end-size(coilsens_set,1)/4,:,ichan)),[]);
            %                             title(['Chan ', num2str(ch_id_track(ichan))]);
            imshow(abs(multiR_coilsens_set(:,:,ichan)),[]);
            title(['Chan ', num2str(ichan)]);
        end
        colormap(used_colormap);
    else
        figure;
        set(gcf,'name','Coil Sensitivities Single RAID');
        for ichan = 1:nchan
            if mod(nchan,4) == 0
                subplot(4,nchan/4,ichan);
            else
                subplot(4,floor(nchan/4)+1,ichan);
            end
            %                             % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
            %                             imshow(abs(multiR_coilsens_set((size(multiR_coilsens_set,1)/4+1):end-size(multiR_coilsens_set,1)/4,:,ichan)),[]);
            %                             title(['Chan ', num2str(ch_id_track(ichan))]);
            imshow(abs(multiR_coilsens_set(:,:,ichan)),[]);
            title(['Chan ', num2str(ichan)]);
        end
        colormap(parula);
    end
end

multiR_snr_opt = zeros(multiR_nrow,multiR_ncol);
for irow = 1:multiR_nrow
    for icol = 1:multiR_ncol
        s_matrix = squeeze(multiR_coilsens_set(irow,icol,:));
        multiR_snr_opt(irow,icol) = sqrt(2)*abs((s_matrix')*inv(multiR_noisecov)*squeeze(multiR_img_matrix(irow,icol,:)))/...
            sqrt((s_matrix')*inv(multiR_noisecov)*s_matrix);
    end
end
multiR_snr_map = abs(multiR_snr_opt);
multiR_snr_map = fliplr(multiR_snr_map);




% snr_ratio = singleR_snr_map./multiR_snr_map;
% snr_ratio_max = max(snr_ratio(:));
% snr_ratio_mean = mean(snr_ratio(:));
% disp('------------------------')
% disp(['   Max SNR Ratio = ' num2str(100*(snr_ratio_max-1)) ' %'])
% disp('------------------------')
% disp(['   Mean SNR Ratio = ' num2str(100*(snr_ratio_mean-1)) ' %'])
% 
% 
% 
% figure;
% subplot(3,2,1)
% imagesc(abs(singleR_snr_map)); colorbar; axis image; title('Single RAID')
% subplot(3,2,2)
% imagesc(abs(multiR_snr_map)); colorbar; axis image; title('Multi RAID')
% subplot(3,2,3)
% imagesc(abs(singleR_snr_map./multiR_snr_map)); colorbar; axis image; title('Abs(Single/Multi) RAID')
% subplot(3,2,4)
% imagesc(abs(singleR_snr_map-multiR_snr_map)); colorbar; axis image; title('Abs(Single - Multi) RAID')
% subplot(3,2,5)
% hist(abs(singleR_snr_map(:))); colorbar
% subplot(3,2,6)
% hist(abs(multiR_snr_map(:))); colorbar
% 
% 
% 
% 



















if 0
    
    figure;
    set(gcf,'name','IMAGE_OBJ --> IMAGE (k-space)')
    for icoil = 1:nchan
        subplot(4,4,icoil)
        hist(real(kdata(:,:,icoil)));
        title(['coil #' num2str(icoil)]);
    end
    figure;
    set(gcf,'name','PRESCAN_OBJ --> IMAGE (k-space)')
    for icoil = 1:nchan
        subplot(4,4,icoil)
        hist(real(knoise(:,:,icoil)));
        title(['coil #' num2str(icoil)]);
    end
    
    figure;
    set(gcf,'name','IMAGE_OBJ --> IMAGE (image domain)')
    for icoil = 1:nchan
        subplot(4,4,icoil)
        mydata = kdata(:,:,icoil);
        kdata_img = MRifft(mydata,[1,2]);
        hist(real(kdata_img));
        title(['coil #' num2str(icoil)]);
    end
    figure;
    set(gcf,'name','PRESCAN_OBJ --> IMAGE (image domain)')
    for icoil = 1:nchan
        subplot(4,4,icoil)
        mydata2 = knoise(:,:,icoil);
        knoise_img = MRifft(mydata2,[1,2]);
        hist(real(knoise_img));
        title(['coil #' num2str(icoil)]);
    end
    
end

if 0
    % adjust dwell time
    testdata = kdata(:,:,1);
    testnoise = knoise(:,:,1);
    figure;
    subplot(1,2,1)
    hist(real(testdata)); colorbar
    subplot(1,2,2)
    hist(real(testnoise)); colorbar
end

if 0
kdata_noisecov = zeros(nchan);
for iCh = 1:nchan,
    for jCh = 1:nchan,
        kdata_noisecov(iCh,jCh)=sum(sum(kdata(:,:,iCh).*conj(kdata(:,:,jCh))))/(size(kdata,1)*size(kdata,2));
    end;
end;
kdata_noisecoeff = zeros(nchan);
for iCh = 1:nchan,
    for jCh = 1:nchan,
        kdata_noisecoeff(iCh,jCh)= kdata_noisecov(iCh,jCh)/sqrt(kdata_noisecov(iCh,iCh)*kdata_noisecov(jCh,jCh));
    end;
end;

knoise_noisecov = zeros(nchan);
for iCh = 1:nchan,
    for jCh = 1:nchan,
        knoise_noisecov(iCh,jCh)=sum(sum(knoise(:,:,iCh).*conj(knoise(:,:,jCh))))/(size(knoise,1)*size(knoise,2));
    end;
end;
knoise_noisecoeff = zeros(nchan);
for iCh = 1:nchan,
    for jCh = 1:nchan,
        knoise_noisecoeff(iCh,jCh)= knoise_noisecov(iCh,jCh)/sqrt(knoise_noisecov(iCh,iCh)*knoise_noisecov(jCh,jCh));
    end;
end;

knoise2_noisecov = zeros(nchan);
for iCh = 1:nchan,
    for jCh = 1:nchan,
        knoise2_noisecov(iCh,jCh)=sum(sum(knoise2(:,:,iCh).*conj(knoise2(:,:,jCh))))/(size(knoise2,1)*size(knoise2,2));
    end;
end;
knoise2_noisecoeff = zeros(nchan);
for iCh = 1:nchan,
    for jCh = 1:nchan,
        knoise2_noisecoeff(iCh,jCh)= knoise2_noisecov(iCh,jCh)/sqrt(knoise2_noisecov(iCh,iCh)*knoise2_noisecov(jCh,jCh));
    end;
end;

knoisesum_noisecov = zeros(nchan);
for iCh = 1:nchan,
    for jCh = 1:nchan,
        knoisesum_noisecov(iCh,jCh)=sum(sum(knoisesum(:,:,iCh).*conj(knoisesum(:,:,jCh))))/(size(knoisesum,1)*size(knoisesum,2));
    end;
end;
knoisesum_noisecoeff = zeros(nchan);
for iCh = 1:nchan,
    for jCh = 1:nchan,
        knoisesum_noisecoeff(iCh,jCh)= knoisesum_noisecov(iCh,jCh)/sqrt(knoisesum_noisecov(iCh,iCh)*knoisesum_noisecov(jCh,jCh));
    end;
end;
end

if 0
    figure;
    set(gcf,'name','Noise Covariance and Noise Coefficient')
    
    subplot(2,4,1);
    imagesc(abs(kdata_noisecov))
    title('Noise Cov Img')
    subplot(2,4,2);
    imagesc(abs(knoise_noisecov))
    title('Noise Cov Noi')
    subplot(2,4,3);
    imagesc(abs(knoise2_noisecov))
    title('Noise Cov Noi 2')
    subplot(2,4,4);
    imagesc(abs(knoisesum_noisecov))
    title('Noise Cov Noi Sum')
    
    
    subplot(2,4,5);
    imagesc(abs(kdata_noisecoeff))
    title('Noise Coef Img')
    subplot(2,4,6);
    imagesc(abs(knoise_noisecoeff))
    title('Noise Coef Noi')
    subplot(2,4,7);
    imagesc(abs(knoise2_noisecoeff))
    title('Noise Coef Noi 2')
    subplot(2,4,8);
    imagesc(abs(knoisesum_noisecoeff))
    title('Noise Coef Noi Sum')
end


