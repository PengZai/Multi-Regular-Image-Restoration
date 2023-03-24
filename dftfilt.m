function g = dftfilt(f, H)
%DFTFILT Performs frequency domain filtering
%G = DFTFILT(F, H) filters F in the frequency domain using the
%filter transfer function H. The output, G, is the filtered image,
%which has the same size as F.

%Obtain the FFT of the padded input
F = fft2(f, size(H, 1), size(H, 2));
%Perform filtering
g = real(ifft2(H.*F));
%Crop to orginal size
g = g(1:size(f,1), 1:size(f,2));
