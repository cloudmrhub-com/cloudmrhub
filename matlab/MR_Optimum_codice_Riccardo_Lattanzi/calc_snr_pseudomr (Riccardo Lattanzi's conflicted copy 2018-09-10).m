function [snr,g,img] = calc_snr_pseudomr(signalrawdata,noisecov,nreplicas,nchan,recon_method,compute_g,acc,simple_sens)

disp(['     Pseudo multiple replicas beginning (' num2str(nreplicas) ' replicas) ...']);
nrow = size(signalrawdata,1);
ncol = size(signalrawdata,2);

if compute_g == 0
    g = [];
else
    g = zeros(size(signalrawdata));
end

if nchan == 1
    image_stack = zeros([size(signalrawdata) nreplicas]);
    img_matrix = MRifft(signalrawdata,[1,2]);
    img_matrix = img_matrix*sqrt(nrow*ncol);
    for ireplica = 1:nreplicas
        gaussian_whitenoise = sqrt(noisecov)*randn(size(signalrawdata)); % it's the noise bandwidth in the single-channel case
        noisy_kspace = signalrawdata + gaussian_whitenoise;
        image_stack(:,:,ireplica) = MRifft(noisy_kspace,[1,2])*sqrt(nrow*ncol);
    end
    image_noise = std(real(image_stack),0,3);
    %     image_signal = mean(real(image_stack),3);
    %     snr = image_signal./image_noise;
    snr = real(img_matrix)./image_noise;
else
    switch recon_method
        case 'opt'
            [V,D] = eig(noisecov);
            corr_noise_factor = V*sqrt(D)*inv(V);
            
            img_matrix = MRifft(signalrawdata,[1,2]);
            img_matrix = img_matrix*sqrt(nrow*ncol);

            if simple_sens
                reference_image = sqrt(sum(abs(img_matrix).^2,3));
                coilsens_set = img_matrix./repmat(reference_image,[1 1 nchan]);
            else
                [recon,coilsens_set]=adapt_array_2d(img_matrix);
            end
            disp('***  ...coil sensitivities computed...');
            
            image_stack = zeros([nrow ncol nreplicas]);
            for ireplica = 1:nreplicas
                %         gaussian_whitenoise = (rand(size(signalrawdata))-.5);
                %         gaussian_whitenoise = randn(size(signalrawdata));
                gaussian_whitenoise = sqrt(0.5)*(randn(size(signalrawdata)) +1i*randn(size(signalrawdata)));
                gaussian_whitenoise = reshape(gaussian_whitenoise,[nrow*ncol nchan]);
                gaussian_whitenoise = corr_noise_factor*(gaussian_whitenoise.');
                gaussian_whitenoise = reshape((gaussian_whitenoise.'),[nrow ncol nchan]);
                noisy_kspace = signalrawdata + gaussian_whitenoise;
                noisy_img = MRifft(noisy_kspace,[1,2]);
                noisy_img = noisy_img*sqrt(nrow*ncol);
                
                for irow = 1:size(signalrawdata,1)
                    for icol = 1:size(signalrawdata,2)
                        s_matrix = squeeze(coilsens_set(irow,icol,:));
                        image_stack(irow,icol,ireplica) = abs((s_matrix')*inv(noisecov)*squeeze(noisy_img(irow,icol,:)));
                        if ireplica == 1
                            img(irow,icol) = abs((s_matrix')*inv(noisecov)*squeeze(img_matrix(irow,icol,:)));
                        end
                    end
                end
                disp(['    Replica #' num2str(ireplica) ' done']);
            end
            image_noise = std(real(image_stack),0,3);
            %     image_signal = mean(real(image_stack),3);
            %     snr = image_signal./image_noise;
            snr = real(img)./image_noise;
    end
end
