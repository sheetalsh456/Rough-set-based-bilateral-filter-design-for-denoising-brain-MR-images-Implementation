original_image = brain_mri;
original_image = rgb2gray(original_image);
%imshow(I);
[row,col] = size(original_image);
noise = 15*randn(row, col);
noise = uint8(noise);
noisy_image = original_image + noise;
I = double(original_image)/255.0;
H = double(noisy_image)/255.0;
N = 3;
granule_size = 2;

disp('Trilateral Filter');
TF = trilateral(noisy_image, N, granule_size);
%imshow(TF);
MSE = sum(sum((I-TF).^2))/(row*col);
PSNR = 10*log10(1.0/MSE);
fprintf('MSE: %f ', MSE);
fprintf('\nPSNR: %f dB \n\n', PSNR);

disp('Bilateral Filter');
BF = bilateral(H, 2, [3, 0.2]);
%imshow(BF);
MSE = sum(sum((I-BF).^2))/(row*col);
PSNR = 10*log10(1.0/MSE);
fprintf('MSE: %f ', MSE);
fprintf('\nPSNR: %f dB \n', PSNR);