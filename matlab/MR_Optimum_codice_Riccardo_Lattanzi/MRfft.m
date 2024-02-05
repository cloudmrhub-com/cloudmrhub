function output = MRfft(input, dim)

% function output = MRfft(input, dim)
%
% MR fft convention which takes
% input:
% with the center pixel in the middle
%
% output:
% DC Fourier coefficient at the center

% Ernie Yeh ENY 2/14/02


tmp = ifftshift(input);

for idim = dim,
    tmp = fft(tmp,[],idim);
end

output = fftshift(tmp);