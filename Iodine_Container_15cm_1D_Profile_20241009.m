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
load('data_Water_Phantom_Recon_20240716', 'im_low')
im_low_water = im_low; 
figure; imagesc(im_low); axis off; axis tight; axis equal; colormap gray; colorbar  
load('data_Iodine_Container_15cm_Phantom_Recon_20241003', 'im_low')
im_low_gel = im_low; 


P_water_actual = 90.1597;  
P_air_actual = 2.39;       

HU_water = 0;
HU_air = -1000;

m = (HU_air - HU_water) / (P_air_actual - P_water_actual);
b = HU_water - m * P_water_actual;

convertToHU = @(P) m * P + b;

HU_image_water = convertToHU(im_low_water);

figure;
imagesc(HU_image_water); axis off; axis tight; axis equal; colormap gray; colorbar
title('Water Hounsfield Units');
% 
HU_image_gel = convertToHU(im_low_gel); 
figure;
imagesc(HU_image_gel); axis off; axis tight; axis equal; colormap gray; colorbar
title('Gel Hounsfield Units');
 
%% Plot profile for Wax Phantom

% Define the center and the size of the extraction region
center_y_gel = size(HU_image_gel, 1) / 2; 
profile_width = 50; 
profile_height = 30; 

vertical_offset = 100;  % Positive value moves the line down, negative moves it up
new_center_y_gel = center_y_gel + vertical_offset;

rect_gel = [1, round(new_center_y_gel - profile_height/2), size(HU_image_gel, 2), profile_height];
im_rect_gel = imcrop(HU_image_gel, rect_gel);

profile_avg_gel = mean(im_rect_gel, 1);

save('profile_data_iodine_container_15cm_1D_Profile_20241009_analysis.mat', 'rect_gel', 'profile_avg_gel');

figure
imshow(HU_image_gel, [])
rectangle('Position', rect_gel, 'EdgeColor', 'r', 'LineWidth', 2)  
title('Gel Phantom with Profile Region')

figure
plot(profile_avg_gel, 'r')
% ylim([0 120]);
% xlim([0 2000]);
xlabel('Pixel Position')
ylabel('Intensity')
title('Gel Phantom 1D Profile')

%% Plot profile for Water Phantom

center_y_water = size(HU_image_water, 1) / 2; 
rect_water = [1, round(center_y_water - profile_height/2), size(HU_image_water, 2), profile_height];
im_rect_water = imcrop(HU_image_water, rect_water);
profile_avg_water = mean(im_rect_water, 1); 
save('profile_data_Water_Phantom_Recon_20241009.mat', 'rect_water', 'profile_avg_water');

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

%% Apodization 

winSize = 10; 
hammingWin = hamming(winSize);
hammingWin = hammingWin / sum(hammingWin);
profile_avg_gel = conv(profile_avg_gel, hammingWin, 'same');

profile_avg_water = conv(profile_avg_water, hammingWin, 'same');

%% Plot profiles together

figure;
hold on;
plot(profile_avg_gel, 'r', 'DisplayName', 'Gel Phantom');
plot(profile_avg_water, 'b', 'DisplayName', 'Water Phantom');

if length(profile_avg_gel) ~= length(profile_avg_water)
    max_length = max(length(profile_avg_gel), length(profile_avg_water));
    x_gel = linspace(1, length(profile_avg_gel), length(profile_avg_gel));
    x_water = linspace(1, length(profile_avg_water), length(profile_avg_water));
    x_new = linspace(1, max_length, max_length);
    
    profile_avg_gel_interp = interp1(x_gel, profile_avg_gel, x_new, 'linear', 'extrap');
    profile_avg_water_interp = interp1(x_water, profile_avg_water, x_new, 'linear', 'extrap');
    
    rmse = sqrt(mean((profile_avg_gel_interp - profile_avg_water_interp).^2));
else
    rmse = sqrt(mean((profile_avg_gel - profile_avg_water).^2));
end

fprintf('Root Mean Squared Error (RMSE): %.4f\n', rmse);

xlabel('Pixel Position');
ylabel('Intensity');
title('1D Profile of Gel Phantom vs Water Phantom');
legend('show');
hold off;

cd Gel_Phantom/results 
print(gcf, 'Gel_Water_comp_20241009_Analysis.png', '-dpng', '-r300');
cd ../ 
cd ../ 


 
