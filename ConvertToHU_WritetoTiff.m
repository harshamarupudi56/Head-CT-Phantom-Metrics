% Convert to HU and write to TIFF 

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


% Load your data
load('data_Iodine_Water_Validation_3_20241015_Recon_Binned.mat', 'img_low_binned');
img = (img_low_binned(:,:,200:210));  

% Convert to HU 
P_water_actual = 93;  
P_air_actual = -10;       

HU_water = 0;
HU_air = -1000;

m = (HU_air - HU_water) / (P_air_actual - P_water_actual);
b = HU_water - m * P_water_actual;

convertToHU = @(P) m * P + b;

img = convertToHU(img);

figure; imagesc(img(:,:,5)); colormap gray; axis off; axis tight; axis equal;

slice_to_save = img(:,:,5);

tiff_filename = '/home/sriharsha.marupudi/Desktop/PCD/1112023/data_Iodine_Water_Validation_3_20241015_Recon_Binned.tif';  % Change to your desired path and filename


% Convert the image to a suitable format for saving if necessary
% Ensure the data type is compatible for TIFF; typically double is fine
if isfloat(slice_to_save)
    % If the image contains HU values, it's typically a double.
    % You may want to scale if you need specific ranges.
    slice_to_save = uint16(slice_to_save); % Convert to uint16 if needed; adjust based on HU range.
end

% Save the image as a TIFF file without compression
imwrite(slice_to_save, tiff_filename, 'Compression', 'none');

% Display a message indicating that the image has been saved
disp(['Image saved as TIFF: ', tiff_filename]);