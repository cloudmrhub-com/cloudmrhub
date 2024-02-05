function [noisecov,noise_bandwidth] = calc_noise_cov(noise,bw_correction)



nchan = size(noise,3);
% noisecov = zeros(nchan);
% for jj = 1:nchan,
%     tmp1 = noise(:,:,jj);
%     eta_j = tmp1(:);
%     for kk = 1:nchan,
%         tmp2 = noise(:,:,kk);
%         eta_k = tmp2(:);
%         noisecov(jj,kk) = 0.5 * [var(eta_j + eta_k) + i*var(eta_j + i*eta_k) ...
%             - (1+i)*( var(eta_j) + var(eta_k) )];
%     end;
% end;


% Siemens method
noisecovbis = zeros(nchan);
for iCh = 1:nchan,
    for jCh = 1:nchan,
        noisecovbis(iCh,jCh)=sum(sum(noise(:,:,iCh).*conj(noise(:,:,jCh))))/(size(noise,1)*size(noise,2));
    end;
end;
noisecov = noisecovbis;

% norm = sqrt(diag(diag(noisecovbis))); % extract diagonal elements on a diag-matrix
% noisecovbis2 = noisecovbis + eye(size(norm))*max(max(abs(norm)))/1000; % apply subtle changes to correlation matrix (numerical stability)
% norm = sqrt(diag(diag(noisecovbis2))); % extract diagonal elements
% noisecovbisnorm = inv(norm)*noisecovbis2*inv(norm); %norm matrix
% cor_norm_abs = abs(noisecovbisnorm);  % absolute values disp('Magnitude Values of Correlation Matrix = ')



if bw_correction
    noise_bandwidth = mrir_noise_bandwidth(noise,1);
    if ( noise_bandwidth < 0.74 ) | ( noise_bandwidth > 0.81 ),
        warning('noise bandwidth is too low or too high; this data is unlikely to be pure noise');
    end;
    noisecov = noisecov/noise_bandwidth;
else
    noise_bandwidth = 1;
end
