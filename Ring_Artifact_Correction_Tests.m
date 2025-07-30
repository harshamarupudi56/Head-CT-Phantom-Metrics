close all
clear all
clc
restoredefaultpath
read_data                  = 0 ;
fit_data                   = 0 ;
mfilename = '/home/sriharsha.marupudi/TIGRE-master/MATLAB';
addpath(genpath(mfilename));
current_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/' ;
addpath(genpath(current_file));
mfilename = 'home/sriharsha.marupudi/Desktop/PCD/1112023/Gel_Phantom/readingtools/';
addpath(genpath(mfilename));
use_par = 1; %use parralelisation
addpath '/home/sriharsha.marupudi/Desktop/PCD/1112023'
results_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/Gel_Phantom/results/';
mfilename =  '/home/sriharsha.marupudi/MATLAB Add-Ons/Collections/TV-L1 Image Denoising Algorithm/';
addpath(genpath(mfilename));

%% Input data
%-1- filenames
filename_phantom    = '/gpfs_projects/sriharsha.marupudi/PCD/Harsha/2024/July/Acquisitions_20240716/ACR_Phantom_20240716/SP/';
filename_openbeam    = '/gpfs_projects/sriharsha.marupudi/PCD/Harsha/2024/July/Acquisitions_20240716/STC/';

cd (current_file)
initiate_Gel_Phantom_20240716

%%
if (fit_data)
    cd (STC_file)
    load('fit_data_ACR_SP25_20240716.mat','A_low','B_low','t');
    cd (current_file)

    % A_low = A_low(roix,roiy);
    % B_low = B_low(roix,roiy);

    cd data/
    save('fit_data_ACR_SP25_20240716.mat','A_low','B_low','t'); 
    cd ../
else
    cd data/
    load('fit_data_ACR_SP25_20240716.mat','A_low','B_low','t'); 
    cd ../
end

%% Read the images
AC = 0;
if(read_data)
    calibration = 0;
    cropping = 1; 
    % Read phantom:
    cd(filename_phantom)
    n_frames_per = 4;
    read_data_1E_average;
    cd(current_file)
    Vp_low = V_low;
    cd data/
    save('data_phantom_ACR_Phantom_20240716_Filter_Test','Vp_low','-v7.3')
    cd ../
else
    cd (current_file)
    cd data/
    load('data_phantom_ACR_Phantom_20240716_Filter_Test','Vp_low')
    cd ../
end
disp(' ********** Finish Reading Data ***********')

%%
display = 0; 
if(display)
figure(1)

    
    for ip = 1:length(angles)
        im_1(:,:) = Vp_low(:,:,ip);
        imagesc(im_1); 
        s = sprintf('Low, %d/%d',ip,nb_proj)
        title(s)
        pause(0.01)

    end

end

disp(' ********** Finish Displaying Data ***********')

%% Bad pixel detection map
ring_correction  = 1;

if (ring_correction)

    clear im;
    clear Vnc_low;
    close all
    
    Slow = A_low;
    threshold=3;
    for iiz = 1:size(Vp_low,1)
        iiz
        ie = 1; % low 
        peak_locations
        Vp_low(iiz,locations,:)  = nan;
        A_low(iiz,locations)  = nan;
        B_low(iiz,locations)  = nan;
    end
end

tic
parfor i = 1:nb_proj
    Vp_low(:,:,i) = correctBadPixels(Vp_low(:,:,i));
end
toc
A_low  = correctBadPixels(A_low);
B_low  = correctBadPixels(B_low);
disp(' ********** Finish bad pixel map locator ***********')

disp('bad pixel correction time:')
toc

% figure
% imagesc(Vp_low(:,:,end))  ; axis off; axis tight; axis equal;  colormap gray
disp(' ********** Finish bad pixel correction ***********')

%% This is to do the same correction for the calibration data

Vn_low = (1./B_low).*log(Vp_low./A_low);
Vn_low = real(Vn_low);

%% Remove vertical strips 
ring_correction  = 1;

if (ring_correction)

    clear im;
    clear Vnc_low;
    close all
    Vnc_low = Vn_low;
    mask_3D = zeros(size(Vnc_low));
    
    Slow = mean(Vn_low,3);
    threshold=3;
    for iiz = 1:size(Vn_low,1)
        iiz
        ie = 1; % low 
        peak_locations
        Vnc_low(iiz,locations,:)  = nan;
    end
end
Vnc_low(55:60,:,:) = nan; 


%
tic
parfor i = 1:nb_proj
    Vnc_low(:,:,i) = correctBadPixels(Vnc_low(:,:,i));
end
toc

disp('bad pixel correction time:')
toc

%%
display = 0;
if(display)
    figure(2)
    for i = 1:nb_proj
        subplot(2,1,1)
        imagesc(Vn_low(:,:,i))  ; axis off; axis tight; %axis equal
        title(sprintf('projection# %d, low',i))
        subplot(2,1,2)
        imagesc(Vnc_low(:,:,i))  ; axis off; axis tight; %axis equal
        title(sprintf('projection# %d, low',i))
        pause(0.01)
    end
end 

%% Denoise projections 
sigma = 0.5; 
for is = 1:size(Vn_low,3) 
    image(:,:) = Vnc_low(:,:,is); 
    Vnc_low_proj = imgaussfilt(image,sigma); 
    Vnc_low_d(:,:,is) = Vnc_low_proj; 
end 

Vn_low = Vnc_low_d ; 

disp(' ********** Gaussian Denoised ***********')
figure; imagesc(Vn_low(:,:,100)); colormap gray; axis off;  
title('Denoised Gaussian Projection');

%% Ring correction again 
ring_correction  = 1;

if (ring_correction)

    clear im;
    close all
    mask_3D = zeros(size(Vnc_low));

    Slow = mean(Vnc_low,3);
    threshold=3;
    for iiz = 1:size(Vn_low,1)
        iiz
        ie = 1; % low 
        peak_locations
        Vnc_low(iiz,locations,:)  = nan;
    end
end
Vnc_low(55:60,:,:) = nan; 

%
tic
parfor i = 1:nb_proj
    Vnc_low(:,:,i) = correctBadPixels(Vnc_low(:,:,i));
end
toc

disp('bad pixel correction time:')
toc

%% Median filter 
% img_d = zeros(size(Vnc_low));
% for is = 1:size(Vnc_low, 3)
%     image = Vnc_low(:, :, is); 
%     slice_d = medfilt2(image, [5 5]); 
%     img_d(:, :, is) = slice_d;
% end
% 
% Vnc_low = img_d;


img_d = zeros(size(Vnc_low));
for is = 1:size(Vnc_low, 3)
    image = Vnc_low(:, :, is);           % Extract the 2D slice

    for row = 1:size(image, 1)
        img_d(row, :, is) = medfilt1(image(row, :), 7); % Apply median filter to each row
    end
end

Vnc_low = img_d;

disp('********** Median Filtered Along X-Direction ***********');

%% Sinogram Visulization 
% figure; 
% for slice_idx = 1:size(Vnc_low, 1)
%     sinogram = squeeze(Vnc_low(slice_idx, :, :));
%     sinogram = rot90(sinogram); 
%     imagesc(sinogram); colormap gray;axis off;
%     title(sprintf('Sinogram Slice %d', slice_idx));
%     pause(0.01); 
% end


figure();
slice_idx = size(Vnc_low, 1) / 2; 
sinogram = squeeze(Vnc_low(55, :,:));
sinogram = rot90(sinogram); 
imagesc(sinogram);
colormap gray; axis off;

title('Sinogram');

disp(' ********** Sinogram Visualization Complete ***********');



%% Wavelet Fourier for Ring Artifact 
wavelet = 1 ; 

if wavelet ==1 
    Vn_corrected = Vnc_low; 
    
    wavelet_type = 'db4'; % Daubechies wavelet with 4 coefficients
    level = 2; 
    
    % Wavelet filtering
    for i = 1:size(Vn_corrected, 3)
        [C, S] = wavedec2(Vn_corrected(:,:,i), level, wavelet_type);
        
        C_thresholded = wthresh(C, 's', 2 * max(abs(C)));
        
        Vn_corrected(:,:,i) = waverec2(C_thresholded, S, wavelet_type);
    end
    
    % Parameters for FFT filtering
    filter_radius = 15; 
    n = 2; 
    
    % FFT filtering
    for i = 1:size(Vn_corrected, 3)
        F = fft2(Vn_corrected(:,:,i));
        
        [rows, cols] = size(F);
        [X, Y] = meshgrid(1:cols, 1:rows);
        centerX = cols / 2;
        centerY = rows / 2;
        D = sqrt((X - centerX).^2 + (Y - centerY).^2);
        butterworth_filter = 1 ./ (1 + (D / filter_radius).^(2 * n));
        
        F_filtered = F .* butterworth_filter;
        
        Vn_corrected(:,:,i) = abs(ifft2(F_filtered));
    end
    
    window_size = 3;
    for i = 1:size(Vn_corrected, 3)
        filtered_image(:,:,i) = medfilt2(Vn_corrected(:,:,i), [window_size, window_size]);
    end
    
    figure;
    imagesc(Vnc_low(:,:,125)); colormap gray; axis tight; axis off; title('Before Filtering');
    
    figure;
    imagesc(filtered_image(:,:,125)); colormap gray; axis tight; axis off; title('After Combined Filtering');
end 


%% Titarenko Method 
% https://ieeexplore.ieee.org/document/7452593/authors

ring_correct = 1; 

if ring_correct ==1 
    image = Vnc_low; 
    F = fftshift(fft2(image));
    [rows, cols, channels] = size(F);
    
    % Define the regularization parameter beta
    beta = 0.001;  
    alpha = 4 * beta^2 / (1 - beta^2);
    G_j = zeros(rows, cols);
    for j = -cols/2:cols/2-1
        term = (2 + alpha - sqrt(alpha*(4 + alpha))) / 2;
        G_j(:, j + cols/2 + 1) = alpha / sqrt(alpha*(4 + alpha)) * term^abs(j);
    end
    G_j = repmat(G_j, [1, 1, channels]);
    filtered_F = F .* G_j;
    filtered_image = abs(ifft2(ifftshift(filtered_F)));
    
    figure();
    slice_idx = size(filtered_image, 1) / 2; 
    sinogram_ring = squeeze(filtered_image(55, :,:));
    sinogram_ring = rot90(sinogram_ring); 
    imagesc(sinogram_ring);
    colormap gray; axis off;
end 
%%
close all; 
clc; 
clear im_1 im_2 mask_3D
 
%-4- geometry
CreateGeometry
%geo.nDetector=[3001; 507];     
% alpha  = -0.885; %1.5;%0.25   ; -
% u0     = 6.57+ delta_c(2)*0.1   ;  %12.1;
alpha  = -0.962; %1.5;%0.25   ;
u0     = 6.7778+ delta_c(2)*0.1   ;  %12.1;
v0     = 0; %+ delta_c(1)*0.1  ; 
beta   = 0; 
phi    = 0;
geo.offDetector=[u0;v0];
geo.rotDetector=[alpha*(pi/180);beta*(pi/180);phi*(pi/180)];  % the rotation in radians

pixel_size = 0.1;
voxel_size = pixel_size*geo.DSO/geo.DSD;
size_image_z = size(Vnc_low,1)*pixel_size*geo.DSO/geo.DSD;
size_image_xy = size(Vnc_low,2)*pixel_size*geo.DSO/geo.DSD;

geo.sVoxel = [size_image_xy;size_image_xy;size_image_z]; %25  % number of voxels %(vx) <-- this for 100 um voxel size
geo.nVoxel = round(geo.sVoxel./ voxel_size);     


% geo.nVoxel = [1500;1500;197]; %25  % number of voxels %(vx) <-- this for 100 um voxel size
% geo.sVoxel = [150;150;19.7];      

geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)

type =  'hann';% 'hann'; %'shepp-logan'; %'ram-lak';
%% 
filter_test = 0;  
if filter_test ==1
    filters = {'ram-lak', 'shepp-logan', 'hann','hamming','cosine'};
    
    proj1 = single(Vnc_low);
    proj1(proj1 < 0) = 0;
    
    slice_num = 95;
    
    for idx = 1:length(filters)
        type = filters{idx};
    
        imgFDK = FDK(proj1, geo, angles, 'filter', type); 
        img_low = imgFDK .* 1000;
    
        disp([' ********** Finish IMAGE RECONSTRUCTION with ', type, ' filter ***********'])
    
        figure;
        imagesc(img_low(:, :, slice_num)); 
        colormap gray;
        axis off; axis tight; axis equal;
        title(['Reconstructed Image Slice ', num2str(slice_num), ' with ', type, ' filter']);
    end
end 

%% Filtered Back Projection  
proj1 = single(Vnc_low);
proj1(proj1<0) = 0;
size(proj1);
imgFDK = FDK(proj1,geo,angles,'filter',type); 
img_low = imgFDK.*1000;

disp(' ********** Finish IMAGE RECONSTRUCTION ***********')



%%
% Display all slices
% figure(5)
% for is = 1:size(img_low,3) 
%     imagesc(img_low(:,:,is))  ; axis off; axis tight; axis equal; colormap gray 
%     s = sprintf('Low %d slice',is);
%     title(s)
%     pause(0.01)
% end

im_low = mean(img_low(:,:,90:100),3);
figure; imagesc(im_low); colormap gray; colormap gray; axis off; axis tight; axis equal; 
title('Reconstructed Image');

figure; imagesc(img_low(:,:,95)); colormap gray; colormap gray; axis off; axis tight; axis equal; 
title('Reconstructed Image Slice');
% caxis([0,200]); 


%% Gaussian denoising slice by slice
sigma = 1;                       
img_d = zeros(size(img_low));   
for is = 1:size(img_low, 3)
    image = img_low(:, :, is);          
    slice_d = imgaussfilt(image, sigma); 
    img_d(:, :, is) = slice_d;           
end

img_lowg = img_d;

disp('********** Gaussian Filtered ***********');

im_lowg = mean(img_lowg(:, :, 90:100), 3);
figure; imagesc(im_lowg); colormap gray; axis off; axis tight; axis equal;
title('Denoised Image Gaussian');

figure; imagesc(img_lowg(:,:,95)); colormap gray; colormap gray; axis off; axis tight; axis equal; 
title('Reconstructed Image Slice Gaussian');


% cd Gel_Phantom/results/
% print(gcf, 'ACR_Phantom_Recon_20240716_Filter_Test', '-dpng', '-r300');
% cd ../
% cd ../

%% Binning 512x512 
img_low = img_lowg;
img_low_binned = rebinImage(img_low,512,512); 
im_low_binned = mean(img_low_binned(:, :, 90:100), 3);

figure;
imshow(img_low(:, :, 95), []);
title('Original Slice');

figure; 
imshow(img_low_binned(:, :, 95), []);
title('Binned Slice (512x512)');

% cd Gel_Phantom/results/
% print(gcf, 'ACR_Phantom_Recon_20240716_Filter_Test_Binned', '-dpng', '-r300');
% cd ../
% cd ../

%%
% cd Gel_Phantom/results/
% save('data_ACR_Phantom_Recon_20240716_Filter_Test',"img_low","im_low",'-v7.3')
% save('data_ACR_Phantom_Recon_Binned_20240716_Filter_Test','img_low_bin','im_low_binned','-v7.3')
% cd ../
% disp('Done Saving')