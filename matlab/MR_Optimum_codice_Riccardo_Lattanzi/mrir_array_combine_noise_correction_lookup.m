function [Ameas_per_Smeas, Atrue_per_Strue, Smeas_per_Strue, Ameas_per_Atrue, Ameas_per_Strue] ...
    = mrir_array_combine_noise_correction_lookup(snr_true, Ncha, varargin)
%MRIR_ARRAY_COMBINE_NOISE_CORRECTION_LOOKUP
%
% [snr_meas, snr_true] = mrir_array_combine_noise_correction_lookup(snr_true, Ncha)
%
%
% example:
%
%   mrir_array_combine_noise_correction_lookup(0, 1);


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

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jul/01
% $Id: mrir_array_combine_noise_correction_lookup.m,v 1.1 2009/02/01 22:02:07 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  VERBOSE = 0;

  USE__MATLAB_SYM = 1;
  if ( nargin >= 3 ),
    USE__MATLAB_SYM = varargin{1};
  end;

  if ( ~USE__MATLAB_SYM && ~exist('chgm', 'file') ),
    dstr = '"chgm.m" not available -- using MATLAB symbolic toolbox';
    disp(sprintf('<!> [%s]: %s', mfilename, dstr));
    USE__MATLAB_SYM = 1;
  end;

  if ( USE__MATLAB_SYM && ~license('test','symbolic_toolbox') ),
    dstr = 'MATLAB symbolic toolbox not available -- using "chgm.m"';
    disp(sprintf('<!> [%s]: %s', mfilename, dstr));
    USE__MATLAB_SYM = 0;
  end;

  if ( ~license('test','symbolic_toolbox') && ~exist('chgm', 'file') ),
    error('no methods available for special function evaluation');
  end;

  if ( VERBOSE ),
  if ( USE__MATLAB_SYM ),
    dstr = 'evaluating with "hypergeom"';
  else,
    dstr = 'evaluating with "chgm"';
  end;
  disp(sprintf('==> [%s]: %s', mfilename, dstr));
  end;


  %==--------------------------------------------------------------------==%

  Nevaluations = numel(snr_true);

  % user inputs the "true" SNR, a.k.a., "A_n / sigma"
  % (make a copy whose contents will be altered)
  Atrue_per_Strue = snr_true;

  % leading constant term for Mn_bar, the measured average magnitude
  % intensity
  %%% (Constantinides et al., 1997, eq. 2)
  alpha = prod([ 1 : 2 : (2*Ncha-1) ]) / (  2^(Ncha-1) * factorial(Ncha-1)  );

  %==--------------------------------------------------------------------==%

  t0 = clock;

  % measured average magnitude intensity in terms of the confluent
  % hypergeometric function
  %%% (Constantinides et al., 1997, eq. 2)


  for entry = 1:Nevaluations,

    if ( USE__MATLAB_SYM ),
      % evaluate confluent hypergeometric function with "hypergeom", from
      % the MATLAB symbolic toolbox --- slower, but returns all numerical
      % values
      hg = hypergeom( -0.5, Ncha, -0.5 * Atrue_per_Strue(entry).^2 );
    else,
      % evaluate confluent hypergeometric function with "chgm", translated
      % from fortran by the 'f2matlab' project --- faster, but returns NaNs
      % for some larger values
      [a, b, x, hg] = chgm( -0.5, Ncha, -0.5 * Atrue_per_Strue(entry).^2 );
    end;

    % "measured" SNR, a.k.a., "M_n / sigma"
    Ameas_per_Strue(entry) = alpha * sqrt(pi/2) * hg;
  end;

  t1 = clock;
  runtime_seconds = etime(t1,t0);

  if ( VERBOSE ),
    TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                   fix(runtime_seconds/60/60), ...
                   rem(fix(runtime_seconds/60), 60), ...
                   rem(fix(runtime_seconds), 60));

    dstr = sprintf('total evaluation time = %s', TIME);
    disp(sprintf('<t> [%s]: %s', mfilename, dstr));
  end;


  %==--------------------------------------------------------------------==%

  invalid_entries = find(~isfinite(Ameas_per_Strue));

  Atrue_per_Strue(invalid_entries) = [];   % user-supplied
  Ameas_per_Strue(invalid_entries) = [];   % calculated


  % compute the correction factor for the noise standard deviation, a.k.a.,
  % "sigma_M_n / sigma"
  %%% (Henkelman, 1985, eq. 7; Constantinides et al., 1997, eq. 3)
  Smeas_per_Strue = sqrt(2*Ncha + Atrue_per_Strue.^2 - Ameas_per_Strue.^2);

  Ameas_per_Atrue = Ameas_per_Strue ./ Atrue_per_Strue;

  Ameas_per_Smeas = Ameas_per_Strue ./ Smeas_per_Strue;




  if ( Nevaluations == 1 && isempty(invalid_entries) ),
    fprintf(['\n [A_%d = %2.2f sigma]:', ...
             '    M_%d = %2.2f sigma, sigma_M_%d = %0.3f sigma\n\n'], ...
             Ncha, Atrue_per_Strue, ...
             Ncha, Ameas_per_Strue, ...
             Ncha, Smeas_per_Strue);
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonnyreb/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_combine_noise_correction_lookup.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
