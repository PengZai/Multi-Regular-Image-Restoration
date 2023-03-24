clc;
clear all;
close all;


%orignal_img = imread('.\dataSet\cameraman.tif');
orignal_img = imread('.\Ximage\bga3_0_055.jpg');
orignal_img = imresize(orignal_img, [256, 256]);


if(length(size(orignal_img))==3)
    orignal_img=im2double(rgb2gray(orignal_img));
else
    orignal_img=im2double(orignal_img);
end

[n,m] = size(orignal_img);

sigma = 0.5;
middle = n/2 + 1;
h = zeros(size(orignal_img));
for i=-4:4
   for j=-4:4
      h(i+middle,j+middle)= (1/(1+i*i+j*j));
   end
end

degrade_img =imfilter(orignal_img,h,'circular','conv')+ sigma*randn(m,n);




%cut off distance (0~0.5)*PQ(1)
%D0 = 0.14*m;

%频谱能量估计，确定分区阈值D0
orignal_img_fft2 = abs(fft2(orignal_img));
orignal_img_energy = sum(sum(orignal_img_fft2));
for i=0.14:0.001:0.5
    estimation_D = i*m;
    estimation_LP_filter = lpfilter('gaussian', n, m, estimation_D);
    estimation_fft2 = abs(fft2(dftfilt(orignal_img, estimation_LP_filter)));
    estimation_energy = sum(sum(estimation_fft2));
    beta = estimation_energy/orignal_img_energy;
    if beta >= 0.6
        D0 = estimation_D;
        break
    end
end

LP_filter = lpfilter('gaussian', n, m, D0);
HP_filter = hpfilter('gaussian', n, m, D0);

show_img_fft = log(1 + abs(fftshift(fft2(orignal_img))));
%show_img_fft = log(1 + abs(fft2(orignal_img)));
show_LP_filter = abs(fftshift(LP_filter));
show_HP_filter = abs(fftshift(HP_filter));



smooth_img = dftfilt(orignal_img, LP_filter);
smooth_degrade = dftfilt(degrade_img, LP_filter);
detail_img = dftfilt(orignal_img, HP_filter);
detail_degrade = dftfilt(degrade_img, HP_filter);



f = smooth_img + detail_img;


figure, imshow(orignal_img, []); title('orignal');

figure, imshow(smooth_img, []); title('smooth_img');
figure, imshow(detail_img, []); title('detail_img');



figure, imshow(degrade_img, []); title('degrade');
figure, imshow(smooth_degrade, []); title('smooth_degrade');
figure, imshow(detail_degrade, []); title('detail_degrade');


