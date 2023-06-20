function Rn =calucaltenoisecovariance(noise,bandwidth)
%freq,phase,coils
%update 09/24/2020 Riccardo Lattanzi

nchan = size(noise,3);
if (nargin<2)
    bandwidth=1;
end
  Rn = zeros(nchan);
  noise_samples = reshape(noise,[size(noise,1)*size(noise,2) nchan]);
  noise_samples = noise_samples.';
  Rn = 1/(2*size(noise_samples,2))*(noise_samples*noise_samples');
  Rn = Rn/bandwidth;

end
