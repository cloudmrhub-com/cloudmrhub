% -------------------------------------------------------------------------
% 
%  SNR Analysis Toolbox Version 4.0
%
% ----------------------------------
% (snr_toolbox_batch.m)
%
% reconstructs images in SNR units from MR raw data
%
% Riccardo Lattanzi <riccardo.lattanzi@nyumc.org>
% September 13, 2013
%
% Notes:
%
% - channel id shifted if there is a channel_id = 0 (RL 10-28-2010)
% - 2D acceleration for SENSE (RL 01-21-2011)
% - routine to read VD11 data, e.g. from Skyra (RL 07-13-2012)
% - RSS SNR without noise covariance matrix (RL 07-13-2012)
% - No need to specify whether is VB or VD data (RL 09-13-2013)
% - Extract noise info from signal rawdata for noise covariance calculation (RL 09-13-2013)
%
% KNOWN ISSUES:
% - multislice option not implemented for VB13-17 raw data because it is not clear how dimensions are
%                     ordered (user has to squeeze the matrix before processing the raw data if more than 1 slice)
% - normally consecutive numbers are assigned to the used receive channels in the raw data,  
%                     but for some 7T datasets the raw data matrix has all 32 receive channels
%                     with zero-entries for non-used channels 
% - currently assumes that if there is a "0" channel_id it is always at the end
% - to avoid issues in loops, when there is a channel ID = 0, all channels are renumbered by adding 1 (in VB13-17 data)
% - the code assume that body_coil data is associated with channel ID = 0, but need to verify this
%
% -------------------------------------------------------------------------

prismadata = 1;

% ---------------------
% RECONSTRUCTION METHOD
% ---------------------
% available: 
%    'rss'     = root sum of square reconstruction
%    'rss_psi' = root sum of square reconstruction with full noise covariance matrix
%    'opt'     = optimal Roemer reconstruction
%    'sense'   = SENSE parallel imaging reconstruction
%    'grappa'  = GRAPPA parallel imaging reconstruction (not yet working)

recon_method = 'opt';

% --------------------
%  PSEUDO-MR SETTINGS
% --------------------
use_pseudo_mr = 0; % 1 --> calculate SNR map using pseudo multiple replicas
nreplicas = 128;   % number of replicas syntethically generated

% -------------------
%  GENERAL SETTINGS
% -------------------
which_slice = 1;           % select which slice to use in rawdata with multiple slices
use_noise_ref = 0;
calc_noise_BW = 0;
noise_BW_correction = 0;   % 1 --> correct noise statistics for the system noise bandwidth
compute_noise_coeff = 1;   % 1 --> compute the noise correlation matrix from the noise covariance matrix
split_noise = 0;           % 1 --> model real and imag noise separately
compute_coil_snr = 0;      % 1 --> compute an SNR map for each coil
mag_correction = 0;        % 1 --> apply magnitude bias correction
noise_mag_biased = 0;      % 1 --> SNR magnitude correction in both numerator and denominator, 0 --> only numerator
normalize_sens = 1;        % 1 --> normalize coil sensitivities
simple_sens = 1;           % 1 --> coil sensitivities calculated as the coil images normalized by the RSS reconstruction (FAST)
% 0 --> coil sensitivities manipulated to remove noise with the method used in Siemens scanners (SLOW)

decimate_data = 1;         % 1 --> decimate fully sampled data to simulate acceleration
compute_g_factor = 1;      % 1 --> compute geometry factor for parallel imaging

save_log_file = 0;         % save a copy if the command window in a txt file
user_label = 'RL';         % user signature
warningsoffflag = 1;       % 1 --> set off some MATLAB warnings

% ---------------
%  PLOT SETTINGS
% ---------------
plot_coil_images = 1;
plot_noise_coeff = 1;
plot_noisecov = 1;
plot_reconstructed_image = 1;
plot_g_factor = 0;
plot_coil_sensitivities = 1;
plot_snr_map = 1;
used_colormap = 'jet';
show_colorbar = 1;
% clims = [0 400];
clims = [];
usedatamask = 1; % mask data for better g-factor display

if plot_g_factor && ~compute_g_factor
    plot_g_factor = 0;
end

% -------------------------------------
%  PARALLEL IMAGING SIMULATION SETTINGS
% -------------------------------------
acc = 1;                   % acceleration factor (iPAT) ** only integer values **
% acc = [2 2];                 % 2D acceleration factor [Freq Phase] ** only integer values **

% if (~decimate_data) && (length(acc) > 1)
%     disp('------------------------------------------------------------------');
%     disp('*** ERROR: 2D acceleration available only when decimating data ***');
%     disp('------------------------------------------------------------------');
%     keyboard
% end

switch recon_method
    case 'rss'
        decimate_data = 0;
        acc = 1;
    case 'rss_psi'
        decimate_data = 0;
        acc = 1;
    case 'opt'
        decimate_data = 0;
        acc = 1;
    case 'grappa'
        if acc == 1
            recon_method = 'opt';
            decimate_data = 0;
        end
end

if compute_g_factor && (acc == 1)
    compute_g_factor = 0;
    plot_g_factor = 0;
end

sens_file = []; % file with coil sensitivities. If [] the user will be asked to browse the disks and select them

mat_rawdata = 0;
grappa_kernel = [4,3];
nacsy=32;                  % autocalibration lines used in GRAPPA
nacsx=32;

% ----------------------
%  INITIALIZATION
% ----------------------
K_ICE_AMPL_SCALE_FACTOR = 3200; % ICE multiplies data by this quantity

if warningsoffflag,
    warning('off','MATLAB:nearlySingularMatrix')
    warning('off','MATLAB:divideByZero')
end

homedir = './';
datadir = [homedir 'data/'];
logdir      = [homedir 'logfiles/'];
addpath(logdir);
if save_log_file
    logfilename = [logdir 'snrtool_log_user_' user_label '_' datestr(now,30) '.txt'];
    diary(logfilename);
    disp('----------------------------------------------------------------');
    disp(['a copy of the command window will be stored in: ./' logfilename]);
    disp('----------------------------------------------------------------');
end

% NOISE and DATA files. If [] the user will be asked to browse the disks and select them

noise_file = [datadir 'RawData_Nova24ChCoil/meas_MID407_SNR24_Tra_BW200_FID17932.dat']; %[datadir 'RawData_Nova24ChCoil/meas_MID404_SNR24_Noi_BW200_FID17929.dat'];
data_file = [datadir 'RawData_Nova24ChCoil/meas_MID407_SNR24_Tra_BW200_FID17932.dat'];
% noise_file = [datadir 'bore_liner/meas_MID114_GRE_BC_SNR_Tra_RF000v_FID28590.dat'];
% data_file = [datadir 'bore_liner/meas_MID111_GRE_BC_SNR_Tra_RF625v_FID28587.dat'];
% noise_file = [datadir 'SNR_QED_Sag/meas_MID231_gre_SNR_noi_50mmSl_RF_000v_FID11325.dat'];
% data_file = [datadir 'SNR_QED_Sag/meas_MID233_gre_SNR_sag_50mmSl_RF_180v_FID11327.dat'];
% noise_file = [datadir 'bore_liner/meas_MID771_GRE_24Ch_SNR_Tra_RF000v_FID98609.dat'];
% data_file = [datadir 'bore_liner/meas_MID770_GRE_24Ch_SNR_Tra_RF680v_FID98608.dat'];
% noise_file = [datadir 'radial_dGEMRIC/meas_MID412_tse_MP_noise_FID173799.dat'];
% data_file = [datadir 'radial_dGEMRIC/meas_MID405_T1_p2_SAG_B1corr_DAN_FID173792.dat'];

% noise_file = [datadir 'Greg/2011_5_4_PARKS_SNR/meas_MID803_SNR_Tra_BW300_NOISE_FID25017.dat'];
% data_file = [datadir 'Greg/2011_5_4_PARKS_SNR/meas_MID800_SNR_Tra_BW300_AX_FID25014.dat'];
% data_file = [datadir 'Greg/2011_5_4_PARKS_SNR/meas_MID801_SNR_Tra_BW300_SAG_FID25015.dat'];
% data_file = [datadir 'Greg/2011_5_4_PARKS_SNR/meas_MID802_SNR_Tra_BW300_COR_FID25016.dat'];

% noise_file = [datadir '3T_Skyra/Shoulder16ch_small_SNRRawData/meas_MID254_SNR_Sag_BW300_FID13612.dat'];

noise_file = [];
data_file = [];

% ----------------------------
%  LOAD FILES AND READ RAWDATA
% ----------------------------


if ~use_noise_ref
    
    
    % -- read noise --
    if isempty(noise_file)
        [noise_file,noise_path] = uigetfile('*.dat','Select the noise raw file');
        disp('***  loading noise data from:');
        disp(['   ' fullfile(noise_path,noise_file)]);
        disp('---');
        noise_file = fullfile(noise_path,noise_file);
    end
    
%     [image_obj noise_obj phasecor_obj refscan_obj refscanPC_obj RTfeedback_obj phasestab_obj] = mapVBVD(noise_file,'doAverage');
   image_obj = mapVBVD(noise_file,'doAverage');
   
   if prismadata
       image_obj = image_obj{2};
   end
   
   kdata = image_obj.image();
    clear image_obj
    totalslices = size(kdata,5);
    if which_slice > totalslices,
        which_slice = 1;
        disp('*** WARNING: whichslice was out of range so it was set to 1');
    end
    noiserawdata = permute(squeeze(kdata(:,:,:,:,which_slice,:,:,:,:,:,:,:,:,:,:,:)),[1,3,2]); % [nfreq nphase ncoil]
    noiserawdata = noiserawdata*K_ICE_AMPL_SCALE_FACTOR;
    
    % -- read data --
    if isempty(data_file)
        [data_file,data_path] = uigetfile('*.dat','Select the data raw file');
        disp('***  loading MR signal from:');
        disp(['   ' fullfile(data_path,data_file)]);
        disp('---');
        data_file = fullfile(data_path,data_file);
    end
    
    %     [image_obj noise_obj phasecor_obj refscan_obj refscanPC_obj RTfeedback_obj phasestab_obj] = mapVBVD(data_file,'doAverage');
    image_obj = mapVBVD(data_file,'doAverage');
    
    if prismadata
        image_obj = image_obj{2};
    end
    
    kdata = image_obj.image();
    clear image_obj
    totalslices = size(kdata,5);
    if which_slice > totalslices,
        which_slice = 1;
        disp('*** WARNING: whichslice was out of range so it was set to 1');
    end
    signalrawdata = permute(squeeze(kdata(:,:,:,:,which_slice,:,:,:,:,:,:,:,:,:,:,:)),[1,3,2]); % [nfreq nphase ncoil]
    signalrawdata = signalrawdata*K_ICE_AMPL_SCALE_FACTOR;
    clear kdata
    
else
    
    % -- read signal --
    if isempty(data_file)
        [data_file,data_path] = uigetfile('*.dat','Select the data raw file');
        disp('***  loading MR signal from:');
        disp(['   ' fullfile(data_path,data_file)]);
        disp('---');
        data_file = fullfile(data_path,data_file);
    end
    
    
    % [image_obj noise_obj phasecor_obj refscan_obj] = mapVBVD(noise_file);
    % [image_obj noise_obj] = mapVBVD(noise_file);
%     [image_obj noise_obj phasecor_obj refscan_obj refscanPC_obj RTfeedback_obj phasestab_obj] = mapVBVD(data_file,'doAverage');
    twix_obj = mapVBVD(data_file);
    
    noise_obj = twix_obj.noise.unsorted;   
    kdata = twix_obj.image();
    clear image_obj
    totalslices = size(kdata,5);
    if which_slice > totalslices,
        which_slice = 1;
        disp('*** WARNING: whichslice was out of range so it was set to 1');
    end
    signalrawdata = permute(squeeze(kdata(:,:,:,:,which_slice,:,:,:,:,:,:,:,:,:,:,:)),[1,3,2]); % [nfreq nphase ncoil]
    signalrawdata = signalrawdata*K_ICE_AMPL_SCALE_FACTOR;

    
    %     if ~isempty(refscan_obj.dataSize) % It is a GRAPPA dataset
    %         kref = refscan_obj;
    %         acs = permute(squeeze(kref(:,:,:,:,which_slice)),[1,3,2]); % [nfreq nphase ncoil]
    %         clear refscan_obj
    %     end
    
    kdata = noise_obj();
%     if isempty(noise_obj.dataSize) && use_noise_ref
%         use_noise_ref = 0;
%         disp('*** WARNING: noise reference is not available so use_noise_ref was set to 1');
%     end
    noiserawdata = permute(squeeze(kdata(:,:,:,:,which_slice,:,:,:,:,:,:,:,:,:,:,:)),[1,3,2]); % [nfreq nphase ncoil]
    noiserawdata = noiserawdata*K_ICE_AMPL_SCALE_FACTOR;

    clear noise_obj
    clear kdata
end


% if decimate_data && (length(ch_id_track) == 1)
%     disp('**WARNING** cannot reconstruct undersampled data with one coil,');
%     decimate_data = 0;
%     acc = 1;
%     disp('the flag "decimate_data" has been set to zero and "acc" to one')
% end
if decimate_data && (length(acc) > 1)
    originalrawdata = signalrawdata;
    signalrawdata = signalrawdata(1:acc(1):end,1:acc(2):end,:);
else
    originalrawdata = signalrawdata;
    signalrawdata = signalrawdata(:,1:acc:end,:);
end
nrow = size(signalrawdata,1);
ncol = size(signalrawdata,2);
nchan = size(signalrawdata,3);

fullrow = size(originalrawdata,1);
fullcol = size(originalrawdata,2);

% ------------------------------------------ %
%  COMPUTE SNR MAPS FOR SINGLE CHANNEL DATA  %
% ------------------------------------------ %

if nchan == 1 % body coil or single channel case
    if ch_id_track == 0
        disp('    *** COMPUTING SNR FOR BODYCOIL DATA ***');
    else
        disp(['    *** COMPUTING SNR FOR SINGLE CHANNEL DATA (channel #' num2str(ch_id_track) ')'])
    end
    
    mag_correction = 0;
    img_matrix = MRifft(signalrawdata,[1,2]);
    img_matrix = img_matrix*sqrt(nrow*ncol);

    
    if plot_reconstructed_image
        figure;
        set(gcf,'name','Single-Channel Reconstructed Image');
        % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
        imshow(abs(img_matrix((size(img_matrix,1)/4+1):end-size(img_matrix,1)/4,:)),[]);
        colormap(used_colormap);
        colorbar
        if ch_id_track == 0
            title('Body Coil');
        else
            title(['Chan ', num2str(ch_id_track)]);
        end
    end
    if noise_BW_correction,
        noise_bandwidth = mrir_noise_bandwidth(noiserawdata,1);
        if ( noise_bandwidth < 0.6 ),
            warning('noise bandwidth is too low; this data is unlikely to be pure noise');
        end;
        %         noisestd = std(noiserawdata(:))/noise_bandwidth;
        noisestd = sqrt((noiserawdata(:)')*noiserawdata(:)/length(noiserawdata(:)))/noise_bandwidth;
    else
        %         noisestd = std(noiserawdata(:));
        noisestd = sqrt((noiserawdata(:)')*noiserawdata(:)/length(noiserawdata(:)));
        noise_bandwidth = 1;
    end
    if use_pseudo_mr
        [snr_map,g_null] = calc_snr_pseudomr(signalrawdata,noisestd,nreplicas,1,recon_method,0,1,simple_sens);
    else
        snr_map = sqrt(2)*abs(img_matrix/noisestd);
    end
    snr_map = snr_map((size(snr_map,1)/4+1):end-size(snr_map,1)/4,:);

% ------------------------------------------ %
%  COMPUTE SNR MAPS FOR MULTI-CHANNEL DATA   %
% ------------------------------------------ %

else
    
    %% Calculate Noise Covariance Matrix
    %  (from noise only prescan)
    % ***************************************** %
    
    disp('***  start computing noise statistics...');
%     noisecovjp = mrir_array_stats_matrix(noiserawdata, 'cov',noise_BW_correction); % Jonathan's code (equivalent result)

%     noiserawdata = MRifft(noiserawdata,[1,2]);
%     noiserawdata = noiserawdata*sqrt(256);
    
    [noisecov,noise_bandwidth] = calc_noise_cov(noiserawdata,noise_BW_correction);
    if split_noise
        % scale the noise covariance by sqrt(2) to use a model where real and imag
        % noise are independent with the same sigma
        noisecov = noisecov / sqrt(2);
    end
    if compute_noise_coeff
        noise_coeff = zeros(size(noisecov));
        for itemp = 1:size(noisecov,1)
            for jtemp = 1:size(noisecov,1)
                noise_coeff(itemp,jtemp) = noisecov(itemp,jtemp)/sqrt(noisecov(itemp,itemp)*noisecov(jtemp,jtemp));
            end
        end
    end
    disp('***  ...noise statistics computed');
    if plot_noisecov
        disp('***  ...plotting magnitude of noise covariance matrix');
        figure;
%         pcolor(abs(noisecov));
        imshow(abs(noisecov.'),[]); colormap(jet)
        colorbar
        title('NOISE COVARIANCE MATRIX');
    end
    if plot_noise_coeff
        disp('***  ...plotting magnitude of noise coefficient matrix');
        figure;
        imshow(abs(noise_coeff.'),[]); colormap(jet)
        colorbar
        title('NOISE COEFFICIENTS MATRIX');
    end
    
    %% Generate Images and SNR Maps for Individual Coils
    % Using Fourier transform and the diagonal elements of the noise
    % covariance matrix
    % ******************************************* %
    
    disp('***  start reconstructing individual coil images...');
    img_matrix = MRifft(signalrawdata,[1,2]);
    % reconstruct individual coils' images and apply FFT scale factor
    % iFFT scales data by 1/N, so I need to scale back, as the noise covariance
    % matrix is calculated in k-space (not Fourier transformed)
    img_matrix = img_matrix*sqrt(nrow*ncol);

    disp('***  ...individual coil images reconstructed');
    % Plot images for individual coils
    if plot_coil_images
        if (nchan/4 < 1)
            figure;
            set(gcf,'name','Individual Coil Images');
            for ichan = 1:nchan
                subplot(1,nchan,ichan);
%                 % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
%                 imshow(abs(img_matrix((size(img_matrix,1)/4+1):end-size(img_matrix,1)/4,:,ichan)),[]);
%                 if ~vd11_flag
%                     title(['Chan ', num2str(ch_id_track(ichan))]);
%                 end
                imshow(abs(img_matrix(:,:,ichan)),[]);
                title(['Chan ', num2str(ichan)]);
            end
            colormap(used_colormap);
        else
            figure;
            set(gcf,'name','Individual Coil Images');
            for ichan = 1:nchan
                if mod(nchan,4) == 0
                    subplot(4,nchan/4,ichan);
                else
                    subplot(4,floor(nchan/4)+1,ichan);
                end
%                 % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
%                 imshow(abs(img_matrix((size(img_matrix,1)/4+1):end-size(img_matrix,1)/4,:,ichan)),[]);
%                 title(['Chan ', num2str(ch_id_track(ichan))]);
                imshow(abs(img_matrix(:,:,ichan)),[]);
                title(['Chan ', num2str(ichan)]);
                
            end
            colormap(used_colormap);
        end
    end
    % Compute SNR and plot SNR maps for individual coils
    if compute_coil_snr
        coil_snr_data = zeros(nrow/2,ncol,nchan);
        if (nchan/4 < 1)
            figure;
            set(gcf,'name','Individual Coil SNR Maps');
            for ichan = 1:nchan
                noisesamples = noiserawdata(:,:,ichan);
%                 noisestd = std(noisesamples(:))/noise_bandwidth;
%                 noisestd = sqrt((noisesamples(:)')*noisesamples(:)/length(noisesamples(:)))/noise_bandwidth;
%                 coil_snr_map = sqrt(2)*abs(img_matrix((size(img_matrix,1)/4+1):end-size(img_matrix,1)/4,:,ichan)/noisecov(ichan,ichan));
                coil_snr_map = sqrt(2)*abs(img_matrix((size(img_matrix,1)/4+1):end-size(img_matrix,1)/4,:,ichan)/sqrt(noisecov(ichan,ichan)));
                coil_snr_data(:,:,ichan) = coil_snr_map;
                
                subplot(1,nchan,ichan);
                % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                imshow(coil_snr_map,[]);
                title(['Chan ', num2str(ch_id_track(ichan))]);
            end
            colormap(used_colormap);
        else
            figure;
            set(gcf,'name','Individual Coil SNR Maps');
            for ichan = 1:nchan
                noisesamples = noiserawdata(:,:,ichan);
%                 noisestd = std(noisesamples(:))/noise_bandwidth;
%                 noisestd = sqrt((noisesamples(:)')*noisesamples(:)/length(noisesamples(:)))/noise_bandwidth;
%                 coil_snr_map = sqrt(2)*abs(img_matrix((size(img_matrix,1)/4+1):end-size(img_matrix,1)/4,:,ichan)/noisecov(ichan,ichan));
                coil_snr_map = sqrt(2)*abs(img_matrix((size(img_matrix,1)/4+1):end-size(img_matrix,1)/4,:,ichan)/sqrt(noisecov(ichan,ichan)));
                coil_snr_data(:,:,ichan) = coil_snr_map;
                if mod(nchan,4) == 0
                    subplot(4,nchan/4,ichan);
                else
                    subplot(4,floor(nchan/4)+1,ichan);
                end
                % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                imshow(coil_snr_map,[]);
                title(['Chan ', num2str(ch_id_track(ichan))]);
            end
            colormap(used_colormap);
        end
    end
    
    % ***************************************** %
    %%% IMAGE RECONSTRUCTION IN SNR UNITS
    % ***************************************** %
    switch recon_method
    
        %% Root Sum of Squares Reconstruction WITHOUT Noise Covariance Matrix
        case 'rss'
            disp('***  start RSS reconstruction...');
            snr_rss = zeros(size(signalrawdata,1),size(signalrawdata,2));

            
            snr_rss_test = zeros(size(signalrawdata,1),size(signalrawdata,2));
            
            if use_pseudo_mr
                [snr_rss,g_rss] = calc_snr_pseudomr(signalrawdata,noisecov,nreplicas,nchan,recon_method,0,1,simple_sens);
            else
                for irow = 1:size(signalrawdata,1)
                    for icol = 1:size(signalrawdata,2)
                        
                        
                        
                        S = squeeze(img_matrix(irow,icol,:));

%                         signalmag = sqrt(2)*abs( S' * S );
% 
%                         noisepower_test = abs(S' * S);
%                         
% %                         noisestd = sqrt((noiserawdata(:)')*noiserawdata(:)/length(noiserawdata(:)))/noise_bandwidth;
% 
%                         snr_rss(irow,icol) = signalmag / (std(noiserawdata(:))*sqrt(noisepower_test));
                        
                        
        signalmag = sqrt(2)*abs( S' * S );
        
        noisepower = abs(S' * noisecov * S);
        
        
        snr_rss(irow,icol) = signalmag / sqrt(noisepower);
                        
                        
%                         snr_rss(irow,icol) = sqrt(2*(squeeze(img_matrix(irow,icol,:))')*inv(diag(diag(noisecov)))*squeeze(img_matrix(irow,icol,:)));
                        %                 snr_rss(irow,icol) = sqrt(2*sqrt(nrow*ncol)*(squeeze(img_matrix(irow,icol,:))')*squeeze(img_matrix(irow,icol,:)));
                    end
                end
            end
            snr_map = real(snr_rss); % it should be real, this is just to remove residual imaginary components that would give an error when plotting the SNR map
            snr_map = fliplr(snr_map);
            % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
            snr_map = snr_map((size(snr_map,1)/4+1):end-size(snr_map,1)/4,:);
            
            
            snr_map_test = real(snr_rss_test); % it should be real, this is just to remove residual imaginary components that would give an error when plotting the SNR map
            snr_map_test = fliplr(snr_map_test);
            % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
            snr_map_test = snr_map_test((size(snr_map_test,1)/4+1):end-size(snr_map_test,1)/4,:);
            
            
            
            if mag_correction
                snr_corrected = mrir_array_combine_noise_correction(snr_map, nchan, 'rss', noise_mag_biased);
            end
            disp('***  ...RSS reconstruction done');
            if plot_reconstructed_image
                rss_image = sqrt(sum(abs(img_matrix).^2,3));
                figure;
                set(gcf,'name','RSS Reconstructed Image');
                % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                imshow(abs(rss_image((size(rss_image,1)/4+1):end-size(rss_image,1)/4,:)),[]);
                colormap(used_colormap);
                title('RSS Image');
            end

        %% Root Sum of Squares Reconstruction WITH Noise Covariance Matrix
        case 'rss_psi'
            disp('***  start RSS reconstruction with noise covariance matrix...');
            snr_rss_psi = zeros(size(signalrawdata,1),size(signalrawdata,2));
            if use_pseudo_mr
                [snr_rss_psi,g_rss_psi] = calc_snr_pseudomr(signalrawdata,noisecov,nreplicas,nchan,recon_method,0,1,simple_sens);
            else
                for irow = 1:size(signalrawdata,1)
                    for icol = 1:size(signalrawdata,2)
                        
                        snr_rss_psi(irow,icol) = sqrt(2*(squeeze(img_matrix(irow,icol,:))')*inv(noisecov)*squeeze(img_matrix(irow,icol,:)));
                        %                 snr_rss(irow,icol) = sqrt(2*sqrt(nrow*ncol)*(squeeze(img_matrix(irow,icol,:))')*squeeze(img_matrix(irow,icol,:)));
                    end
                end
            end
            snr_map = real(snr_rss_psi); % it should be real, this is just to remove residual imaginary components that would give an error when plotting the SNR map
            snr_map = fliplr(snr_map);
            % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
            snr_map = snr_map((size(snr_map,1)/4+1):end-size(snr_map,1)/4,:);
            if mag_correction
                snr_corrected = mrir_array_combine_noise_correction(snr_map, nchan, 'rss', noise_mag_biased);
            end
            disp('***  ...RSS reconstruction with noise covariance matrix done');
            if plot_reconstructed_image
                rss_image = sqrt(sum(abs(img_matrix).^2,3));
                figure;
                set(gcf,'name','RSS + Noise Covariance Matrix Reconstructed Image');
                % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                imshow(abs(rss_image((size(rss_image,1)/4+1):end-size(rss_image,1)/4,:)),[]);
                colormap(used_colormap);
                title('RSS + Noise Covariance Matrix Image');
            end
            
        %% Roemer Optimal Reconstruction
        case 'opt'
            disp('***  start Roemer-optimal (OPT) reconstruction...');
            if use_pseudo_mr
                [snr_opt,g_opt,opt_image] = calc_snr_pseudomr(signalrawdata,noisecov,nreplicas,nchan,recon_method,0,1,simple_sens);
            else
                if simple_sens
                    reference_image = sqrt(sum(abs(img_matrix).^2,3));
                    coilsens_set = img_matrix./repmat(reference_image,[1 1 nchan]);
                else
                    [recon,coilsens_set]=adapt_array_2d(img_matrix);
                end
                disp('***  ...coil sensitivities computed...');
                
                if plot_coil_sensitivities
                    if (nchan/4 < 1)
                        figure;
                        set(gcf,'name','Coil Sensitivities');
                        for ichan = 1:nchan
                            subplot(1,nchan,ichan);
%                             % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
%                             imshow(abs(coilsens_set((size(coilsens_set,1)/4+1):end-size(coilsens_set,1)/4,:,ichan)),[]);
%                             title(['Chan ', num2str(ch_id_track(ichan))]);
                            imshow(abs(coilsens_set(:,:,ichan)),[]);
                            title(['Chan ', num2str(ichan)]);
                            
                        end
                        colormap(used_colormap);
                    else
                        figure;
                        set(gcf,'name','Coil Sensitivities');
                        for ichan = 1:nchan
                            if mod(nchan,4) == 0
                                subplot(4,nchan/4,ichan);
                            else
                                subplot(4,floor(nchan/4)+1,ichan);
                            end
%                             % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
%                             imshow(abs(coilsens_set((size(coilsens_set,1)/4+1):end-size(coilsens_set,1)/4,:,ichan)),[]);
%                             title(['Chan ', num2str(ch_id_track(ichan))]);
                            imshow(abs(coilsens_set(:,:,ichan)),[]);
                            title(['Chan ', num2str(ichan)]);
                         end
                        colormap(used_colormap);
                    end
                end
                
                snr_opt = zeros(size(signalrawdata,1),size(signalrawdata,2));
                if plot_reconstructed_image
                    opt_image = snr_opt;
                end
                for irow = 1:size(signalrawdata,1)
                    for icol = 1:size(signalrawdata,2)
                        s_matrix = squeeze(coilsens_set(irow,icol,:));
                        snr_opt(irow,icol) = sqrt(2)*abs((s_matrix')*inv(noisecov)*squeeze(img_matrix(irow,icol,:)))/...
                            sqrt((s_matrix')*inv(noisecov)*s_matrix);
                        if plot_reconstructed_image
                            opt_image(irow,icol) = abs((s_matrix')*inv(noisecov)*squeeze(img_matrix(irow,icol,:)));
                        end
                    end
                end
            end
%             snr_map = abs(snr_opt);
            snr_map = real(snr_opt);
            snr_map = fliplr(snr_map);
%             % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
%             snr_map = snr_map((size(snr_map,1)/4+1):end-size(snr_map,1)/4,:);
            if mag_correction
                snr_corrected = mrir_array_combine_noise_correction(snr_map, nchan, 'opt', noise_mag_biased);
            end
            disp('***  ...Roemer-optimal (OPT) reconstruction done');
            if plot_reconstructed_image
                figure;
                set(gcf,'name','OPT Reconstructed Image');
%                 % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
%                 imshow(abs(opt_image((size(opt_image,1)/4+1):end-size(opt_image,1)/4,:)),[]);
                imshow(abs(opt_image),[]);
                colormap(used_colormap);
                title('OPT Image');
            end
            
            %% SENSE Reconstruction
        case 'sense'
            disp('***  start SENSE reconstruction...');
            if use_pseudo_mr
                [snr_sense,g_sense] = calc_snr_pseudomr(signalrawdata,noisecov,nreplicas,nchan,recon_method,compute_g_factor,acc,simple_sens);
            else
                % --** COMPUTE COIL SENSITIVITIES **--
                if decimate_data
                    if (length(acc)>1)
                        if (mod(acc(1),2)==0)
                            img_matrix=ifftshift(img_matrix,1);
                        end
                        if (mod(acc(2),2)==0)
                            img_matrix=ifftshift(img_matrix,2);
                        end
                        
                    else
                        if (mod(acc,2)==0)
                            img_matrix=ifftshift(img_matrix,2);
                        end
                    end
                    % Coil sensitivities from original rawdata
                    sens_matrix = MRifft(originalrawdata,[1,2]);
                    sens_matrix = sens_matrix*sqrt(size(originalrawdata,1)*size(originalrawdata,2));
                    if simple_sens
                        reference_image = sqrt(sum(abs(sens_matrix).^2,3));
                        coilsens_set = sens_matrix./repmat(reference_image,[1 1 nchan]);
                    else
                        [recon,coilsens_set] = adapt_array_2d(sens_matrix);
                    end
                    disp('***  ...coil sensitivities computed...');
                else
                    if mod(acc,2)==0,
                        img_matrix=ifftshift(img_matrix,2);
                    end
                    % Coil sensitivities from file
                    if isempty(sens_file)
                        [sens_file,sens_path] = uigetfile('*.dat*','Select raw file with coil sensitivities');
                        disp('***  loading sensitivity data from:');
                        disp(['   ' fullfile(sens_path,sens_file)]);
                        disp('---');
                        sens_file = fullfile(sens_path,sens_file);
                    end
                    [raw, noise, ref, phasecor, pc_order, centerlines, lastheader, rest, phasestab, ch_id_track] = ...
                        read_meas_vb13(sens_file);
                    sensrawdata = permute(raw,[1 3 2]);
                    sensrawdata = K_ICE_AMPL_SCALE_FACTOR*sensrawdata;
                    coilsens_set = MRifft(sensrawdata,[1,2]);
                    coilsens_set = coilsens_set*sqrt(size(sensrawdata,1)*size(sensrawdata,2));
                    if simple_sens
                        reference_image = sqrt(sum(abs(coilsens_set).^2,3));
                        coilsens_set = coilsens_set./repmat(reference_image,[1 1 nchan]);
                    else
                        [recon,coilsens_set] = adapt_array_2d(coilsens_set);
                    end
                    disp('***  ...coil sensitivities computed...');
                end
                %             [s_r,g]=sense_std_1d(permute(img_matrix,[2 1 3]),permute(coilsens_set,[2 1 3]),'pinv');
                % --** PLOT COIL SENSITIVITIES **--
                if plot_coil_sensitivities
                    if (nchan/4 < 1)
                        figure;
                        set(gcf,'name','Coil Sensitivities');
                        for ichan = 1:nchan
                            subplot(1,nchan,ichan);
                            % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                            imshow(abs(coilsens_set((size(coilsens_set,1)/4+1):end-size(coilsens_set,1)/4,:,ichan)),[]);
                            title(['Chan ', num2str(ch_id_track(ichan))]);
                        end
                        colormap(used_colormap);
                    else
                        figure;
                        set(gcf,'name','Coil Sensitivities');
                        for ichan = 1:nchan
                            if mod(nchan,4) == 0
                                subplot(4,nchan/4,ichan);
                            else
                                subplot(4,floor(nchan/4)+1,ichan);
                            end
                            % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                            imshow(abs(coilsens_set((size(coilsens_set,1)/4+1):end-size(coilsens_set,1)/4,:,ichan)),[]);
                            title(['Chan ', num2str(ch_id_track(ichan))]);
                        end
                        colormap(used_colormap);
                    end
                end
                
                % --** SENSE RECONSTRUCTION **--
                snr_sense = zeros([fullrow fullcol]);
                
                if plot_reconstructed_image
                    sense_image = snr_sense;
                end
                if compute_g_factor
                    g_sense = snr_sense;
                end
                
                for irow = 1:nrow
                    if (length(acc)>1)
                        freq_set = floor(irow:fullrow/acc(1):(fullrow + 0.5));
                    else
                        freq_set = irow;
                    end
                    for icol = 1:ncol
                        if (length(acc)>1)
                            phase_set = floor(icol:fullcol/acc(2):(fullcol+0.5));
                        else
                            phase_set = floor(icol:fullcol/acc:(fullcol+0.5));
                        end
                        if acc == 1
                            s_matrix = squeeze(coilsens_set(freq_set,phase_set,:));
                        else
                            s_matrix = squeeze(coilsens_set(freq_set,phase_set,:));
                            if length(acc)>1
                                s_matrix = reshape(s_matrix,[length(freq_set)*length(phase_set) nchan]);
                            end
                            s_matrix = s_matrix.';
                        end
                        u_matrix = inv((s_matrix')*inv(noisecov)*s_matrix)*(s_matrix')*inv(noisecov);
                        %                     u_matrix = inv(s_matrix*inv(noisecov)*(s_matrix'))*s_matrix*inv(noisecov);
                        snr_sense(freq_set,phase_set) = reshape(sqrt(2)*(u_matrix)*squeeze(img_matrix(irow,icol,:))./diag(sqrt((u_matrix)*noisecov*(u_matrix'))),[length(freq_set) length(phase_set)]);
                        if plot_reconstructed_image
                            sense_image(freq_set,phase_set) = reshape(u_matrix*squeeze(img_matrix(irow,icol,:)),[length(freq_set) length(phase_set)]);
                        end
                        if compute_g_factor
                            if (length(acc)>1)
                                g_sense(freq_set,phase_set) = reshape( sqrt(acc(1)*acc(2)*diag(inv((s_matrix')*inv(noisecov)*s_matrix)).*diag((s_matrix')*inv(noisecov)*s_matrix)), [length(freq_set) length(phase_set)]);
                            else
                                g_sense(freq_set,phase_set) = reshape( sqrt(acc*diag(inv((s_matrix')*inv(noisecov)*s_matrix)).*diag((s_matrix')*inv(noisecov)*s_matrix)), [length(freq_set) length(phase_set)]);
                            end
                        end
                    end
                end
            end
            %             snr_sense = snr_sense(nrow/4 + 1:end - nrow/4,:);
            snr_map = abs(snr_sense);
            snr_map = fliplr(snr_map);
            % Crop FOV (i.e. remove the effect of the 2-fold frequency
            % oversampling)
            snr_map = snr_map((size(snr_map,1)/4+1):end-size(snr_map,1)/4,:);
            if mag_correction
                snr_corrected = mrir_array_combine_noise_correction(snr_map, nchan, 'sense', noise_mag_biased);
            end
            disp('***  ...SENSE reconstruction done');
            if plot_reconstructed_image
                figure;
                set(gcf,'name','SENSE Reconstructed Image');
                % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                imshow(abs(sense_image((size(sense_image,1)/4+1):end-size(sense_image,1)/4,:)),[]);
                title('SENSE Image');
                colormap(used_colormap);
            end
            if plot_g_factor && compute_g_factor
                figure;
                set(gcf,'name','SENSE g-factor');
                % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                imshow(abs(g_sense((size(g_sense,1)/4+1):end-size(g_sense,1)/4,:)),[]);
                if (length(acc)>1)
                    title(['SENSE g-factor  ACC = ' num2str(acc(1)) ' x ' num2str(acc(2))] );
                else
                    title(['SENSE g-factor  ACC = 1 x ' num2str(acc)] );
                end
                colormap(used_colormap);
                colorbar
                
                figure;
                set(gcf,'name','SENSE 1/g-factor');
                inv_g = 1./g_sense((size(g_sense,1)/4+1):end-size(g_sense,1)/4,:);
                inv_g(isnan(inv_g)) = 0;
                inv_g(isinf(inv_g)) = 0;
                % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                imshow(abs(inv_g),[]);
                if (length(acc)>1)
                    title(['SENSE 1/g-factor  ACC = ' num2str(acc(1)) ' x ' num2str(acc(2))] );
                else
                    title(['SENSE 1/g-factor  ACC = 1 x ' num2str(acc)] );
                end
                colormap(used_colormap);
                colorbar
                
            end
    %% GRAPPA Reconstruction
        case 'grappa'
            disp('***  start GRAPPA reconstruction...');
            if use_pseudo_mr
                [snr_grappa,g_grappa] = calc_snr_pseudomr(signalrawdata,noisecov,nreplicas,nchan,recon_method,compute_g_factor,acc,simple_sens);
            else
            signalrawdata = permute(signalrawdata,[2 1 3]); % PE x FE x Coil
            if decimate_data
                kernel_npoints = grappa_kernel(1)*grappa_kernel(2);
                kernel_y = grappa_kernel(1);
                kernel_x = grappa_kernel(2);
                fullrows = size(originalrawdata,1);
                fullcols = size(originalrawdata,2);
                ctr_kspace_x = (fullrows/2)+1;
                ctr_kspace_y = (fullcols/2)+1;
                
                % autocalibration lines are taken from the fully-sampled data
%                 data_acs=signalrawdata(nrow/2-nacsy/2+1:nrow/2+nacsy/2,ncol/2-nacsx/2+1:ncol/2+nacsx/2,:);
                data_acs = originalrawdata(fullrows/2-nacsx/2+1:fullrows/2+nacsx/2,fullcols/2-nacsy/2+1:fullcols/2+nacsy/2,:);
                data_acs = permute(data_acs,[2 1 3]);
%                 % reconstruction of the fully sampled data
%                 [s_full,cmap_full,p_comb_full]=adapt_array_2d(fftshift(fft2(fftshift(signalrawdata))),noisecov,1);
                
                % k-space grappa reconstruction of the accelerated data (for each coil)
                [sk_acc,ws] = grappa1_2d(signalrawdata,data_acs,acc,grappa_kernel,1);
                disp('*** ...GRAPPA weights calculated...');
                
                if plot_reconstructed_image
                    % adaptive combination of coil images in a composite image
%                     [grappa_image,cmap,p_comb] = adapt_array_2d(fftshift(fft2(fftshift(sk_acc))),noisecov,1);
                    [grappa_image,cmap,p_comb] = adapt_array_2d(MRifft(sk_acc,[1 2]),noisecov,1);
                    grappa_image = grappa_image*sqrt(size(sk_acc,1)*size(sk_acc,2));
                    grappa_image = permute(grappa_image,[2 1 3]);
                    figure;
                    set(gcf,'name','GRAPPA K-space Reconstructed Image');
                    % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                    imshow(abs(grappa_image((size(grappa_image,1)/4+1):end-size(grappa_image,1)/4,:)),[]);
                    title('GRAPPA Image (k-space reconstruction)');
                    colormap(used_colormap);
                    if usedatamask
                        datamask = zeros(size(grappa_image));
                        datamask(abs(grappa_image) > 0.001) = 1;
                        figure;
                        imshow(datamask((size(datamask,1)/4+1):end-size(datamask,1)/4,:),[]);
                    end
                end
                
                % --- begin grappa reconstruction in image domain --- %
                sens_matrix = permute(MRifft(sk_acc,[1 2]),[2 1 3])*sqrt(size(sk_acc,1)*size(sk_acc,2));
                acc_imgs = zeros(fullcols,fullrows,nchan);
                W_matrix_full = zeros(fullcols,fullrows,nchan,nchan);
                
                % zero-fill the raw data to match final image dimensions
                data_acc_zerofilled = zeros(fullcols,fullrows,nchan);
                data_acc_zerofilled(1:acc:fullcols,:,:) = signalrawdata;
                for icoil = 1:nchan
                    %                     % create aliased coil image
                    %                     data_acc_img(:,:,icoil) = fftshift(fft2(fftshift(data_acc_zerofilled(:,:,icoil))));
                    
                    % the weights to be applied to all coils in order to reconstruct the un-aliased image of the i-th coil
                    ws_zerofilled = zeros(fullcols,fullrows,nchan);
                    W_single = zeros(fullcols,fullrows);
                    for jcoil = 1:nchan
                        % ws dimensions are [acc-1]x[nchan*kernel_x*kernel_y]x[nchan]
                        
                        %                         ws_temp = reshape(squeeze(ws(:,kernel_npoints*(jcoil-1)+1:kernel_npoints*jcoil,icoil)),[kernel_y,kernel_x]);
                        ws_temp = reshape(squeeze(ws(:,jcoil:nchan:end,icoil)),[kernel_y,kernel_x]); % @ works only for acc = 2;
                        
                        ws_temp = rot90(rot90(ws_temp));
                        
                        % center and zero-pad the convolution kernel
                        if mod(kernel_y,2) == 0,
                            neg_rows = kernel_y/2;
                            pos_rows = neg_rows;
                        else
                            neg_rows = ceil(kernel_y/2);
                            pos_rows = floor(kernel_y/2);
                        end
                        if mod(kernel_x,2) == 0,
                            neg_cols = kernel_x/2;
                            pos_cols = neg_cols;
                        else
                            neg_cols = ceil(kernel_x/2);
                            pos_cols = floor(kernel_x/2);
                        end
                        for iky = 1:kernel_y
                            if iky <= neg_rows
                                for ikx = 1:kernel_x
                                    if ikx <= neg_cols
                                        ws_zerofilled(ctr_kspace_y-(neg_rows*acc-(iky-1)*acc),ctr_kspace_x-(neg_cols-(ikx-1)),jcoil)=ws_temp(iky,ikx);
                                    else
                                        ws_zerofilled(ctr_kspace_y-(neg_rows*acc-(iky-1)*acc),ctr_kspace_x+(ikx-neg_cols-1),jcoil)=ws_temp(iky,ikx);
                                    end
                                end
                            else
                                for ikx = 1:kernel_x
                                    if ikx <= neg_cols
                                        ws_zerofilled(ctr_kspace_y+(iky-neg_rows-1)*acc,ctr_kspace_x-(neg_cols-(ikx-1)),jcoil)=ws_temp(iky,ikx);
                                    else
                                        ws_zerofilled(ctr_kspace_y+(iky-neg_rows-1)*acc,ctr_kspace_x+(ikx-neg_cols-1),jcoil)=ws_temp(iky,ikx);
                                    end
                                end
                            end
                        end
                        % set the center of the kernel equal to 1
                        if jcoil == icoil
                            ws_zerofilled(ctr_kspace_y-1,ctr_kspace_x-1,jcoil) = 1;
                        end
                        
                        %                       % create and store the weights in image space
                        %                         W_matrix_full(:,:,icoil,jcoil) = fftshift(fft2(fftshift(ws_zerofilled(:,:,jcoil))));
                        W_matrix_full(:,:,icoil,jcoil) = MRifft(ws_zerofilled(:,:,jcoil),[1 2])*sqrt(size(ws_zerofilled,1)*size(ws_zerofilled,2));
                        %                         W_single = W_single + sqrt(fullcols*fullrows)*MRifft(ws_zerofilled(:,:,jcoil),[1 2]).*(sens_matrix(:,:,jcoil)'); % OLD VERSION
                        acc_imgs(:,:,icoil) = acc_imgs(:,:,icoil) + sqrt(size(ws_zerofilled,1)*size(ws_zerofilled,2))*MRifft(ws_zerofilled(:,:,jcoil),[1 2]).*MRifft(data_acc_zerofilled(:,:,jcoil),[1 2])*sqrt(size(data_acc_zerofilled,1)*size(data_acc_zerofilled,2));
                    end
                    %                     W_single = sum(MRifft(ws_zerofilled,[1 2]).*(sens_matrix.'),3);
                    %                     acc_imgs(:,:,icoil) = W_single.*MRifft(data_acc_zerofilled(:,:,icoil),[1 2])*sqrt(fullcols*fullrows); %OLD VERSION
                    
                    %                     acc_imgs(:,:,icoil) = sum(MRifft(ws_zerofilled,[1 2]).*MRifft(data_acc_zerofilled,[1 2]),3);
                end
                
                snr_den = zeros(fullrows,fullcols);
                snr_num = zeros(fullrows,fullcols);
                
                if compute_g_factor
                    gfactor_den = zeros(fullrows,fullcols);
                    gfactor_num = zeros(fullrows,fullcols);
                    g_factor_coil = zeros(fullrows,fullcols,nchan);
                    %                     g_factor_coil_temp = zeros(fullrows*fullcols,nchan);
                end
                
                %                 img_sos = sqrt(sum(abs(acc_imgs).^2,3));
                if plot_reconstructed_image
                    [grappa_image_imdom,cmap,p_comb] = adapt_array_2d(acc_imgs,noisecov,1);
                    p_comb = permute(p_comb,[2 3 1]);
                    grappa_image_imdom = permute(grappa_image_imdom,[2 1 3]);
                    figure;
                    set(gcf,'name','GRAPPA Image-Domain Reconstructed Image');
                    % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                    imshow(abs(grappa_image_imdom((size(grappa_image_imdom,1)/4+1):end-size(grappa_image_imdom,1)/4,:)),[]);
                    title('GRAPPA Image (image-domain reconstruction)');
                    colormap(used_colormap);
                end
                % -- Compute SNR map -- %
                
                %                 for icoil = 1:nchan
                %                         wimg_temp = W_matrix_full(:,:,icoil,:);
                %                         wimg_temp = reshape(wimg_temp,[fullcols*fullrows nchan]);
                %                         g_factor_coil_temp(:,icoil) = sqrt(abs(wimg_temp*noisecov*(wimg_temp')))/sqrt(noisecov(icoil,icoil));
                %                 end
                %                 g_factor_coil_temp = reshape(g_factor_coil_temp,[fullcols fullrows nchan]);
                
                
                for icol = 1:fullcols
                    for irow = 1:fullrows
                        %                         snr_den_temp = (p_comb(:,irow,icol).')*squeeze(W_matrix_full(irow,icol,:))*noisecov*( (p_comb(:,irow,icol).')*squeeze(W_matrix_full(irow,icol,:)) )';
                        %                         snr_den(irow,icol) = snr_den_temp(1,1);
                        %                         snr_num_temp = sqrt(2)*(p_comb(:,irow,icol).')*squeeze((W_matrix_full(irow,icol,:).*data_acc_img(irow,icol,:)));
                        %                         snr_num(irow,icol) = snr_num_temp(1,1);
                        
                        g_factor_coil(irow,icol,:) = sqrt(diag(abs(squeeze(W_matrix_full(icol,irow,:,:))*noisecov*(squeeze(W_matrix_full(icol,irow,:,:))'))))./sqrt(real(diag(noisecov)));
                        
                        gfactor_num(irow,icol) = sqrt(abs((squeeze(p_comb(icol,irow,:)).')*squeeze(W_matrix_full(icol,irow,:,:))*noisecov*( (squeeze(p_comb(icol,irow,:)).')*squeeze(W_matrix_full(icol,irow,:,:)) )'));
                        gfactor_den(irow,icol) = sqrt(abs((squeeze(p_comb(icol,irow,:)).')*eye(nchan)*noisecov*( (squeeze(p_comb(icol,irow,:)).')*eye(nchan) )'));
                        
                        
                    end
                end
                g_grappa = gfactor_num./gfactor_den;
                %                 snr_grappa = snr_num./sqrt(abs(snr_den));
                %                 snr_grappa = snr_grappa(nrow/4 + 1:end - nrow/4,:);
                if plot_g_factor && compute_g_factor
                    figure;
                    set(gcf,'name','GRAPPA g-factor');
                    if usedatamask
                        g_grappa = g_grappa.*datamask;
                    end
                    % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                    imshow(abs(g_grappa((size(g_grappa,1)/4+1):end-size(g_grappa,1)/4,:)),[]);
                    if (length(acc)>1)
                        title(['GRAPPA g-factor  ACC = ' num2str(acc(1)) ' x ' num2str(acc(2))] );
                    else
                        title(['GRAPPA g-factor  ACC = 1 x ' num2str(acc)] );
                    end
                    colormap(used_colormap);
                    colorbar
                    
                    figure;
                    set(gcf,'name','GRAPPA 1/g-factor');
                    inv_g = 1./g_grappa((size(g_grappa,1)/4+1):end-size(g_grappa,1)/4,:);
                    inv_g(isnan(inv_g)) = 0;
                    inv_g(isinf(inv_g)) = 0;
                    % Crop FOV (i.e. remove the effect of the 2-fold frequency oversampling)
                    imshow(abs(inv_g),[]);
                    if (length(acc)>1)
                        title(['GRAPPA 1/g-factor  ACC = ' num2str(acc(1)) ' x ' num2str(acc(2))] );
                    else
                        title(['GRAPPA 1/g-factor  ACC = 1 x ' num2str(acc)] );
                    end
                    colormap(used_colormap);
                    colorbar
                end
                
            end
                
                %                 figure;
                %                 subplot(1,2,1),imshow(abs(s_full),[]);
                %                 title('acc=1 (no acc)')
                %                 subplot(1,2,2),imshow(abs(s_acc),[]);
                %                 title(strcat('acc=',int2str(acc),' (GRAPPA)'));
%             else
                %                 data_acc=signalrawdata;
                %                 data_acs=signalrawdata(nrow/2-nacsy/2+1:nrow/2+nacsy/2,ncol/2-nacsx/2+1:ncol/2+nacsx/2,:);
                %                 % grappa reconstruction of the accelerated data
                %                 [sk_acc,ws] = grappa1_2d(data_acc,data_acs,acc,[4,3],1);
                %                 [s_acc,cmap,p_comb] = adapt_array_2d(fftshift(fft2(fftshift(sk_acc))),noisecov,1);
            end
            disp('***  ...GRAPPA reconstruction done');
    end
end

%% Plotting SNR Maps
if plot_snr_map
    iptsetpref('ImshowInitialMagnification','fit');
    max_snr = max(snr_map(:));
    min_snr = min(snr_map(:));
    mean_snr = mean(snr_map(:));
    figure;
    %         imagesc(snr_map,clims);
    h_im = imshow(snr_map,clims);
    axis image
    colormap(used_colormap);
    if show_colorbar
        colorbar;
    end
    set(gcf,'name','SNR map');
    title_label = char(['SNR [Max = ' num2str(max_snr,'%5.1d') ', Min = ' num2str(min_snr,'%5.1d'), ',Mean = ' num2str(mean_snr,'%5.1d') ']']);
    title(title_label)
    imgui;
    
%     if (recon_method == 'rss')
%     max_snr_test = max(snr_map_test(:));
%     min_snr_test = min(snr_map_test(:));
%     mean_snr_test = mean(snr_map_test(:));
%     figure;
%     %         imagesc(snr_map,clims);
%     h_im = imshow(snr_map_test,clims);
%     axis image
%     colormap(used_colormap);
%     if show_colorbar
%         colorbar;
%     end
%     set(gcf,'name','SNR map TEST');
%     title_label = char(['SNR [Max = ' num2str(max_snr_test,'%5.1d') ', Min = ' num2str(min_snr_test,'%5.1d'), ',Mean = ' num2str(mean_snr_test,'%5.1d') ']']);
%     title(title_label)
%     imgui;
%     end
    
    
    
    if mag_correction
        max_snr_cor = max(snr_corrected(:));
        min_snr_cor = min(snr_corrected(:));
        mean_snr_cor = mean(snr_corrected(:));
        figure;
        %         imagesc(snr_corrected,clims);
        h_im_cor = imshow(snr_corrected,clims);
        axis image
        colormap(used_colormap);
        if show_colorbar
            colorbar;
        end
        set(gcf,'name','SNR map (magnitude corrected)');
        title_label_cor = char(['SNR (mag corrected) [Max = ' num2str(max_snr_cor,'%5.1d') ', Min = ' num2str(min_snr_cor,'%5.1d'), ',Mean = ' num2str(mean_snr_cor,'%5.1d') ']']);
        title(title_label_cor)
        imgui;
    end
end

if save_log_file
    diary off
end

% END snr_toolbox_batch.m
