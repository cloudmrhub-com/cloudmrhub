function [noise_bandwidth, noise_bandwidth_chan, N_power_spectrum_avg] = mrir_noise_bandwidth(noise, varargin)
%MRIR_NOISE_BANDWIDTH  calculate noise spectrum and return effective noise bandwidth
%
% bandwidth = mrir_noise_bandwidth(noise)
%
% 
% references:
%
%  Kellman P, McVeigh ER.
%  Image reconstruction in SNR units: a general method for SNR measurement.
%  Magn Reson Med. 2005 Dec;54(6):1439-47.
%  PMID: 1626157  
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  DISPLAY = 0;
  if ( nargin >= 2 ),
    DISPLAY = varargin{1};
  end;


  %==--------------------------------------------------------------------==%
  
  dims = size(noise); % [freq phase channel]
  power_spectrum = abs(fft(noise, [], 1)).^2;
  power_spectrum_chan = squeeze(mean(power_spectrum, 2));
  power_spectrum_norm = mean(power_spectrum_chan([1:dims(1)/4 3*dims(1)/4:end],:),1); % for each channel takes the mean within the central region only
  norm_power_spectrum_chan = power_spectrum_chan./repmat(power_spectrum_norm, [dims(1), 1]);
  % individual coils' noise equivalent bandwidth
  noise_bandwidth_chan = sum(norm_power_spectrum_chan, 1)./dims(1);
  
  % global noise equivalent bandwidth
  N_power_spectrum_avg_normalized = mean(norm_power_spectrum_chan,2);
  noise_bandwidth = mean(noise_bandwidth_chan,2);
  
  
  if ( DISPLAY ),
      figure('name', mfilename);
      rh = rectangle('position', [0.25*dims(1), 0.001, 0.5*dims(1), 1.098]);
      set(rh, 'FaceColor', 0.9*[1 1 1], 'LineStyle', 'none');
      hold on;
      plot(fftshift(N_power_spectrum_avg_normalized));
      set(gca, 'XTick', [1, dims(1)*[0.25 0.50 0.75 1.0]], 'XTickLabel', ...
          [-0.5, -0.25 0, +0.25 +0.5]);
      set(gca, 'XLim', [1, dims(1)], 'YLim', [0, 1.1]);
      set(gca, 'YGrid', 'on');
      xlabel('frequency (normalized)');
      ylabel('power (DC normalized)');
      title(sprintf('average of normalized noise power spectrum, BW=%2.3f', noise_bandwidth));
      box on;
  end;
  
  
  return;
  


%%%  N = 10000;
%%%
%%%  noise_true = complex(randn(N,1),randn(N,1)) * 3.0;
%%%  noise_freq = fft(noise_true);
%%%
%%%  T = (tukeywin(N, 0.5) + 0.4) / 1.4;
%%%
%%%  noise_filt = abs((noise_freq) .* T) .* exp(i*angle(noise_freq));
%%%
%%%  BW = sum(T.^2)/N
%%%
%%%  noise_meas = (ifft(noise_filt));
%%%
%%%
%%%  var(real(noise_meas)) / var(real(noise_true))

