function [snr_map_corrected, varargout] = mrir_array_combine_noise_correction(snr_map, Ncha_true, varargin)
%MRIR_ARRAY_COMBINE_NOISE_CORRECTION
%
% [snr_map_corrected, snr_correction] = mrir_array_combine_noise_correction(snr_map, Ncha, [METHOD]);
%
%  where METHOD is either 'rss', 'opt', or 'sense'
%
%
% [snr_map_corrected, snr_correction] = mrir_array_combine_noise_correction(snr_map, Ncha, [METHOD], [FLAG__BIASED_SNR]);
% 
% FLAG__BIASED_SNR = 1 corrects numerator *and* denominator,
% FLAG__BIASED_SNR = 0 corrects numerator *only*
%
%
% See also MRIR_ARRAY_COMBINE_NOISE_CORRECTION_LOOKUP.
  
% Henkelman RM.
% Measurement of signal intensities in the presence of noise in MR images.
% Med Phys. 1985 Mar-Apr;12(2):232-3.
% PMID: 4000083
%
% Constantinides CD, Atalar E, McVeigh ER.
% Signal-to-noise measurements in magnitude images from NMR phased arrays.
% Magn Reson Med. 1997 Nov;38(5):852-7.
% PMID: 9358462
%
% Kellman P, McVeigh ER.
% Image reconstruction in SNR units: a general method for SNR measurement.
% Magn Reson Med. 2005 Dec;54(6):1439-47.
% PMID: 16261576

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jun/20
% $Id: mrir_array_combine_noise_correction.m,v 1.2 2009/08/16 00:56:51 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  METHOD_COMBO = 'rss';
  if ( nargin >= 3 ),
    if ( ~isempty(varargin{1}) ),
      METHOD_COMBO = varargin{1};
    end;
  end;

  % SNR map can either have bias in numerator and denominator (FLAG__BIASED_SNR=1),
  % or bias in just the numerator (FLAG__BIASED_SNR=0);
  FLAG__BIASED_SNR = 1;
  if ( nargin >= 4 ),
    FLAG__BIASED_SNR = varargin{2};
  end;
  
  DEBUG = 0;
  global VERBOSE; if ( isempty(VERBOSE) ), VERBOSE = 0; end;


  %==--------------------------------------------------------------------==%

  switch lower(METHOD_COMBO),
   case {'rss', 'cov'},
    Ncha = Ncha_true;
   case {'opt', 'sense'},
    Ncha = 1;
   otherwise,
    error('unknown coil channel combination method: "%s"', METHOD_COMBO);
  end;


  % heuristic limit for SNR correction lookup: any measured SNRs above this
  % value are assumed to not need correction
  SNR_MAX_LIMIT = 1.5 * Ncha;

  snr_max = min( [max(snr_map(:)), SNR_MAX_LIMIT] );
  snr_min = max(snr_map(:));


  % because we only know the analytic expression for the mapping from true
  % SNR to measured SNR, here we build a range of values of the true SNR
  % (a.k.a., An/sigma), compute the expected measured SNR (a.k.a.,
  % Mnbar/sigma_Mn), then use this as a lookup table (with interpolation) to
  % find the true SNR values corresponding to the input measured values.

  Nevaluations = 100;

  % since correction is always a positive subtraction, "snr_max" guaranteed
  % to be large enough
  snr_true = linspace(0, snr_max, Nevaluations);

  %meas/meas true/true meas/true  meas/true  meas/true    
  [snr_meas, snr_true, std_ratio, avg_ratio, snr_bias] = mrir_array_combine_noise_correction_lookup(snr_true, Ncha);

  
  if ( FLAG__BIASED_SNR ),
    snr_uncorr = snr_meas;
  else,
    if ( VERBOSE ),
      disp(sprintf('==> [%s]: correcting magnitude bias in SNR numerator only...', mfilename));
    end;
    snr_uncorr = snr_bias;
  end;

  

  [snr_uncorr, ind_unique] = unique(snr_uncorr);
  snr_true = snr_true(ind_unique);

  correction_subtractive = zeros(size(snr_map));

  mask_corrected = zeros(size(snr_map));
  mask_upperbound = zeros(size(snr_map));
  mask_lowerbound = zeros(size(snr_map));
  mask_error = zeros(size(snr_map));

  snr_lowerbound = min(snr_uncorr);
  snr_upperbound = max(snr_uncorr);


  for pixel = 1:numel(snr_map),
    if ( snr_map(pixel) < snr_lowerbound ),
      correction_subtractive(pixel) = snr_map(pixel);
      mask_lowerbound(pixel) = 1;
    elseif ( snr_map(pixel) > snr_upperbound ),
      correction_subtractive(pixel) = 0;
      mask_upperbound(pixel) = 1;
    else,
      % interpolate to find subtraction:             x                      y              xi
      correction_subtractive(pixel) = interp1(snr_uncorr, (snr_uncorr - snr_true), snr_map(pixel), 'linear');
      mask_corrected(pixel) = 1;
    end;

    if ( isnan(correction_subtractive(pixel)) ),
      correction_subtractive(pixel) = 0;
      mask_corrected(pixel) = 0;
      mask_error(pixel) = 1;
    end;

  end;

  snr_map_corrected = snr_map - correction_subtractive;

  correction_divisive = snr_map ./ snr_map_corrected;

  index_corrected  = find(mask_corrected);
  index_upperbound = find(mask_upperbound);
  index_lowerbound = find(mask_lowerbound);
  index_error      = find(mask_error);


  if ( nargout > 1 ),
    varargout{1} = correction_divisive;
    varargout{2} = correction_subtractive;
  end;

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_combine_noise_correction.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
