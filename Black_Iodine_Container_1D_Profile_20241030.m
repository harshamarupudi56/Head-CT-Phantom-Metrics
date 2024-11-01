close all
clear all
clc
restoredefaultpath
mfilename = '/home/sriharsha.marupudi/TIGRE-master/MATLAB';
addpath(genpath(mfilename));
mfilename = 'home/sriharsha.marupudi/Desktop/PCD/1112023/Gel_Phantom/readingtools/';
results_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/Gel_Phantom/results/';
current_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/' ;
addpath(genpath(current_file));
cd(current_file); 
%% Load data 
load('data_Black_Iodine_Container_Recon_20241029', 'im_low')
im_low_water = im_low; 
HU_image_water = im_low_water; 
%% Plot profile for Water Phantom

center_y_water = size(HU_image_water, 1) / 2; 
original_profile_width = 2200;  
profile_height = 30; 
vertical_offset = -180; 
center_y_water = center_y_water + vertical_offset; 
center_x_water = size(HU_image_water, 2) / 2;

left_extension = 50;  
right_extension = 75;  

new_center_x_water = center_x_water - left_extension + right_extension / 2;
new_profile_width = original_profile_width + left_extension + right_extension;
rect_water = [new_center_x_water - new_profile_width / 2, round(center_y_water - profile_height / 2), new_profile_width, profile_height];
im_rect_water = imcrop(HU_image_water, rect_water);
profile_avg_water = mean(im_rect_water, 1);


save('profile_data_Black_Iodine_Container_Recon_20241030.mat', 'rect_water', 'profile_avg_water')
%% Apodization 

winSize = 5; 
hammingWin = hamming(winSize);
hammingWin = hammingWin / sum(hammingWin);
profile_avg_water = conv(profile_avg_water, hammingWin, 'same');

%%

figure
imshow(HU_image_water, [])
rectangle('Position', rect_water, 'EdgeColor', 'r', 'LineWidth', 2)  
title('Water Phantom with Profile Region')

figure
plot(profile_avg_water, 'b')
% ylim([0 120]);
% xlim([0 2000]);
xlabel('Pixel Position')
ylabel('Intensity')
title('Water Phantom 1D Profile')



cd Gel_Phantom/results 
print(gcf, 'Black_Iodine_Container_20241030_Analysis.png', '-dpng', '-r300');
cd ../ 
cd ../ 


 
