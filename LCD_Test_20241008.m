% Inputs 
% signal present 3D 
% signal absent 3D 
% ground truth 2D 

close all
clear all
clc
restoredefaultpath 

mfilename = '/home/sriharsha.marupudi/TIGRE-master/MATLAB';
addpath(genpath(mfilename));

current_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/';
addpath(genpath(current_file));

mfilename = '/home/sriharsha.marupudi/Desktop/PCD/1112023/Gel_Phantom/readingtools/';
addpath(genpath(mfilename));

addpath '/home/sriharsha.marupudi/Desktop/PCD/1112023';

results_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/Gel_Phantom/results/';

mfilename = '/home/sriharsha.marupudi/MATLAB Add-Ons/Collections/TV-L1 Image Denoising Algorithm/';
addpath(genpath(mfilename));

mfilename = '/home/sriharsha.marupudi/LCD_CT';
addpath(genpath(mfilename));

load_coords_filename =  '/home/sriharsha.marupudi/LCD_CT/data/coordinates/coordinates_20241008.mat';


%% Load your data
% load('data_Iodine_Water_Validation_3_20241015_Recon_Binned.mat', 'img_low_binned');
load('data_Iodine_Water_Validation_3_Noisy_20241015_Recon_Binned.mat','img_low_binned'); 
img = (img_low_binned(:,:,200:210));  

%% Convert to HU 
% P_water_actual = 93;  
% P_air_actual = -10;       
% HU_water = 0;
% HU_air = -1000;
% m = (HU_air - HU_water) / (P_air_actual - P_water_actual);
% b = HU_water - m * P_water_actual;
% convertToHU = @(P) m * P + b;
% img = convertToHU(img);

%% Add noise 
% img = mat2gray(img);
% noisy_img = imnoise(img, 'gaussian', 0.1);   
% 
% figure; 
% subplot(1, 2, 1);
% imshow(img(:,:,5), []); 
% title('Original Image');
% axis off; axis tight; axis equal;
% 
% subplot(1, 2, 2);
% imshow(noisy_img(:,:,5), []); 
% title('Noisy Image');
% axis off; axis tight; axis equal;
% 
% img= noisy_img; 

%% Mask Image to remove container 
[rows, cols, slices] = size(img);   
centerX = round(cols / 2);
centerY = round(rows / 2);
radius = min(rows, cols) / 2;

[X, Y] = meshgrid(1:cols, 1:rows);
circularMask = (X - centerX).^2 + (Y - centerY).^2 <= radius^2;

maskedImage = zeros(size(img));

for k = 1:slices
    currentSlice = img(:,:,k);
    maskedSlice = zeros(size(currentSlice)); 
    maskedSlice(circularMask) = currentSlice(circularMask);  
    maskedImage(:,:,k) = maskedSlice;  
end

figure;
imagesc(maskedImage(:,:,5));  
colormap gray; axis off; axis tight; axis equal; 
title('Masked Image');

%% Threshold masked image to get inserts 
% th = 100; %if not HU 
% th = 300;  % if HU Threshold value
th = 0.5;
diskSize = 5;  % Disk size for morphological operations
areaThreshold = 50;  % Minimum area for keeping objects in mask

maskedImage_ground_truth = mean(maskedImage(:,:,1:10), 3);
binaryMask = maskedImage_ground_truth > th;

% Remove small objects from the binary mask
cleanedMask = bwareaopen(binaryMask, areaThreshold);  

% Apply morphological closing to smooth the mask
se = strel('disk', diskSize);  
cleanedMask = imclose(cleanedMask, se);

cleanedGrayscaleImage = zeros(size(maskedImage_ground_truth));
cleanedGrayscaleImage(cleanedMask) = img(cleanedMask);

figure;
imshow(cleanedGrayscaleImage, []);
colormap(gray);
title('Ground Truth');

ground_truth = cleanedGrayscaleImage;

%% Signal Free Image 
% Apply insert locations to original image and set inserts to 0 

signalFreeMask = ~cleanedMask; 
signalFreeImage = zeros(size(img));  

for k = 1:slices
    currentSlice = maskedImage(:,:,k);
    maskedSlice = zeros(size(currentSlice)); 
    maskedSlice(signalFreeMask) = currentSlice(signalFreeMask);  
    signalFreeImage(:,:,k) = maskedSlice;  
end

figure;
imagesc(signalFreeImage(:,:,5));   
colormap(gray); axis equal; axis tight; axis off; 
title('Signal Free Image');


%% Specify observers to use
observers = {LG_CHO_2D()};
% observers = {LG_CHO_2D(), DOG_CHO_2D(), GABOR_CHO_2D(), ... };

%% Specify base directory and run the measurement
base_directory = '/your/base/directory/path';  % Adjust to your actual base directory

offset = 0; 
n_reader = 10; 
n_inserts = 6; 
insert_r = 30; 
res_table = measure_LCD1(maskedImage, signalFreeImage, observers, ground_truth, offset,n_reader, n_inserts,insert_r,load_coords_filename);
 

%% Plot and summarize results
% Plot results
custom_insert_names = {'5 Rod', '5 Solution','7.5 Rod', '7.5 Solution', '10 Rod','10 Solution'}; % mg/mL of iodine inserts 
set_ylim = [0 1.2];
plot_results1(res_table, [], custom_insert_names);

% Display results
res_table

% Summarize results
if ~is_octave
  groupsummary(res_table, ["observer", "insert_HU", "dose_level"],["mean", "std"])
end
%% or define a custom summary table by printing mean and standard deviation results
nreader = max(res_table.reader);
for i=1:n_inserts
    mean_AUC(i) = mean(res_table.auc([1:nreader]+(i-1)*nreader));
    std_AUC(i) = std(res_table.auc([1:nreader]+(i-1)*nreader));
    mean_snr(i) = mean(res_table.snr([1:nreader]+(i-1)*nreader));
    std_snr(i) = std(res_table.snr([1:nreader]+(i-1)*nreader));
end

insert_HU = res_table.insert_HU(1:nreader:end);
mean_AUC = mean_AUC(:);
std_AUC = std_AUC(:)
mean_snr = mean_snr(:);
std_snr = std_snr(:)
AUC_res = table(insert_HU, mean_AUC, std_AUC, mean_snr, std_snr);
AUC_res


