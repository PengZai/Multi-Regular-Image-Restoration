clc;clear all;close all;

%1 各向异性， 2 各向同性
tv_method = 1;

%图像读取和预处理
%led
orignal_img = imread('./Ximage/led1_0_08.jpg');
canny_threshold=[0.015,0.022];

%bga
% orignal_img = imread('./Ximage/bga3_0_055.jpg');
% canny_threshold=[0.01,0.015];

if(length(size(orignal_img))==3)
    orignal_img=im2double(rgb2gray(orignal_img));
else
    orignal_img=im2double(orignal_img);
end
orignal_img = imresize(orignal_img, [256, 256]);
[n, m] = size(orignal_img);

orignal_img_canny=edge(orignal_img,'canny', canny_threshold);







%创建正则化模型的模糊矩阵H
load kernels.mat
H=k{7};

H_FFT=psf2otf(H,[n,m]);

L1_H = ifft2(H_FFT);
TV_H = H_FFT;

% definde the function handles that compute 
% the blur and the conjugate blur.  R和RT都会使图像变得模糊
%R_once = @(x) real(ifft2(fft2(h_once).*fft2(x)));
R = @(x) real(ifft2(fft2(L1_H).*fft2(x)));
RT = @(x) real(ifft2(conj(fft2(L1_H)).*fft2(x)));

% define the function handles that compute 
% the products by W (inverse DWT) and W' (DWT)
wav = daubcqf(2);
W = @(x) midwt(x,wav,3);
WT = @(x) mdwt(x,wav,3);


%Finally define the function handles that compute 
% the products by A = RW  and A' =W'*R' 
%AT做正变换、A做逆变换，  x' = A(AT(x)), 但是会更加模糊
A = @(x) R(W(x));
AT = @(x) WT(RT(x));


%noise parameter
sigma=0;

% generate noisy blurred observations
y=imfilter(orignal_img,H,'circular','conv')+ sigma*randn(m,n);


%频谱能量估计，确定分区阈值D0
y_fft2 = abs(fft2(y));
y_energy = sum(sum(y_fft2));
for i=0.1:0.01:0.5
    estimation_D = i*m;
    estimation_LP_filter = lpfilter('gaussian', n, m, estimation_D);
    estimation_fft2 = abs(fft2(dftfilt(y, estimation_LP_filter)));
    pause(0.1)
    estimation_energy = sum(sum(estimation_fft2));
    beta = estimation_energy/y_energy;
    if beta>0.8
        D0 = estimation_D;
        break
    end
end

%图像分区
LP_filter = lpfilter('gaussian', n, m, D0);
HP_filter = hpfilter('gaussian', n, m, D0);
smooth_img = dftfilt(y, LP_filter);
detail_img = dftfilt(y, HP_filter);


% regularization parameter
%L1 parameter
lambda_la1 = 7e-4;

% set stop tolerance
tolA = 1e-5;

%TV parameter
lambda_tv=4e-4;



%Run_D_ADMM
[x_admm,iter]=My_D_ADMM_C(smooth_img, TV_H, lambda_tv,tv_method,1e-5);
%Run GPSR_Basic
x_l1 = W(My_GPSR_Basic(detail_img, A, AT, lambda_la1, tolA));
objective = x_admm + x_l1;
objective_canny=edge(objective,'canny', canny_threshold);

%visualisation
figure(1), imshow(orignal_img, []); title('orignal');
figure(2); imshow(orignal_img_canny, []); title('orignal_img canny')
figure(3); imshow(y, []);  title('degrade');

%figure(4), imshow(smooth_img,[]);  title('smooth_img');
%figure(5), imshow(detail_img,[]); title('detail_img');

%figure(6); imshow(x_admm, []);  title('x_admm');
%figure(7); imshow(x_l1, []);  title('L1');

figure(8); imshow(objective, []);  title('objective');
figure(9); imshow(objective_canny, []); title('objective canny');







