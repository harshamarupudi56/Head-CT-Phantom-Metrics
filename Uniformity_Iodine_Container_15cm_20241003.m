%% CT Number and Uniformity 
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
%% Convert to HU 
% threshold_low = 80;  
% threshold_high = 95; 

P_water_actual = 90.1597;  
P_air_actual = 2.39;       

HU_water = 0;
HU_air = -1000;

m = (HU_air - HU_water) / (P_air_actual - P_water_actual);
b = HU_water - m * P_water_actual;

convertToHU = @(P) m * P + b;

HU_image_water = convertToHU(im_low_water);

% HU_image(abs(im_low_water - P_water_actual) <= (threshold_high - threshold_low) / 2) = 0;
figure;
imagesc(HU_image_water); axis off; axis tight; axis equal; colormap gray; colorbar
title('Water Hounsfield Units');
% 
HU_image_gel = convertToHU(im_low_gel); 
figure;
imagesc(HU_image_gel); axis off; axis tight; axis equal; colormap gray; colorbar
title('Gel Hounsfield Units');
 
%%
[mean_values1,sd_values1,error1] = Calc_Uniformity('coordinations_Uniformity_Water_Phantom_20240716',HU_image_water);
[mean_values2,sd_values2,error2] = Calc_Uniformity('coordinations_Uniformity_Water_Phantom_20240716',HU_image_gel);

%%
barWidth = 0.35; 
data = [sd_values1; sd_values2]; 

std_data = [error1; error2]; 

uniformity_labels = {'Phantoms'};
numLabels = numel(uniformity_labels);
x_values = 1:numLabels;
water_color = [0, 0.4470, 0.7410];      % Blue color for Water
gel_color = [0.8500, 0.3250, 0.0980];   % Orange color for Gel
%%
figure;
hold on;
handles = [];

for i = 1:size(data, 2)
    y = data(:, i);
    err = std_data(:, i);
    xPositions = x_values(i) + barWidth * (0:size(data, 1)-1) - barWidth/2; 
    for j = 1:size(data, 1)
        if j == 2
            barColor = gel_color;
            displayName = 'Gel';
        else
            barColor = water_color;
            displayName = 'Water';
        end
        
        h = bar(xPositions(j), y(j), barWidth/2, 'FaceColor', barColor, 'DisplayName', displayName);
        if i == 1
            handles = [handles, h];
        end
        
        errorbar(xPositions(j), y(j), err(j), 'k', 'linestyle', 'none');
    end
end
%%
hold off;
xticks(x_values);
xticklabels(uniformity_labels);
ylabel('Standard Deviation', 'FontSize', 12);
legend(handles, 'Location', 'eastoutside');
title('Standard Deviation of Water and Gel Phantoms');

%%
[mean_values1,sd_values1,error1] = Calc_Uniformity('coordinations_Uniformity_Water_Phantom_20240716',HU_image_water);
[mean_values2,sd_values2,error2] = Calc_Uniformity('coordinations_Uniformity_Water_Phantom_20240716',HU_image_gel);


barWidth = 0.35; 
data = [mean_values1; mean_values2];  

num_values_CT = length(mean_values1); 
x_values_CT = 1:num_values_CT;             

water_color = [0, 0.4470, 0.7410];  % Blue color for Water
gel_color = [0.8500, 0.3250, 0.0980]; % Orange color for Gel

figure;
hold on;

for j = 1:num_values_CT
   
    barWater_CT = bar(x_values_CT(j) - barWidth/2, data(1, j), barWidth, 'FaceColor', water_color, 'DisplayName', 'Water');
    barGel_CT = bar(x_values_CT(j) + barWidth/2, data(2, j), barWidth, 'FaceColor', gel_color, 'DisplayName', 'Gel');
end

hold off;

xticks(x_values_CT);
xticklabels(arrayfun(@(x) sprintf('ROI %d', x), 1:num_values_CT, 'UniformOutput', false));
ylabel('CT Number', 'FontSize', 12);
legend({'Water', 'Gel'}, 'Location', 'eastoutside');
title('CT Number Accuracy of Water and Gel Phantoms');