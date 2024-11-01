close all
clear all
clc
restoredefaultpath
mfilename = '/home/sriharsha.marupudi/TIGRE-master/MATLAB';
addpath(genpath(mfilename));
current_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/' ;
addpath(genpath(current_file));
mfilename = 'home/sriharsha.marupudi/Desktop/PCD/1112023/Gel_Phantom/readingtools/';
addpath(genpath(mfilename));
addpath '/home/sriharsha.marupudi/Desktop/PCD/1112023'
results_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/Gel_Phantom/results/';

%%
use_saved_coords = 0; 
if use_saved_coords == 1
    try
        % Try loading saved coordinates
        load('coordinations_Gd_Gel_Test_2024930', "y4","x4","y3", "x3", "y2", "x2", "y1", "x1", "x_noise", "y_noise");
    catch
        % If loading fails, reset flag to manually pick new coordinates
        use_saved_coords = 0; 
    end
end

load('data_Gadolinium_Gel_Insert_Test_Phantom_Recon_Binned_20241001',"img_low_binned","im_low_binned")


%%
% figure(5)
% for is = 1:size(img_low_binned,3) 
%     imagesc(img_low_binned(:,:,is))  ; axis off; axis tight; axis equal; colormap gray 
%     s = sprintf('Low %d slice',is);
%     title(s)
%     pause(0.01)
% end

%% Crop Image
uSize = [0.2430 0.2430 0.2430]/3; % 0.1 mm voxel size of image
% figure; imagesc(img_low_binned(:,:,74)); colormap gray; axis off; axis tight; axis equal;
figure; imagesc(img_low_binned(:,:,150)); colormap gray; axis off; axis tight; axis equal;
clim([0 300]); 
 
img = img_low_binned;
radius = 5; 
slice_number = 150; 

%% If coordinates are not loaded, perform manual selection using ginput
if use_saved_coords == 0
    title("Gd Rod 117")
    M = round(ginput(1)); 
    x1 = M(:,1); y1 = M(:,2);

    title("Gd Solution 117")
    M = round(ginput(1)); 
    x2 = M(:,1); y2 = M(:,2);

    title("Gd 96")
    M = round(ginput(1)); 
    x3 = M(:,1); y3 = M(:,2);

    title("Background noise")
    M = round(ginput(1)); 
    x_noise = M(:,1); y_noise = M(:,2);

    % Save coordinates for future use
    save('coordinations_Gd_Gel_Test_20241001',"y3", "x3", "y2", "x2", "y1", "x1", "x_noise", "y_noise"); 
end

%% Crop Circular ROI
cropCircularROI = @(img, centerX, centerY, radius) ...
    (sqrt((repmat((1:size(img, 1))', 1, size(img, 2)) - centerY).^2 + ...
    (repmat((1:size(img, 2)), size(img, 1), 1) - centerX).^2) <= radius);

img = img_low_binned(:,:,slice_number);  

ROI1 = img .* cropCircularROI(img, x1, y1, radius);  
ROI2 = img .* cropCircularROI(img, x2, y2, radius);   
ROI3 = img .* cropCircularROI(img, x3, y3, radius);   

ROInoise = img .* cropCircularROI(img, x_noise, y_noise, radius);  % background noise

mean_ROI1 = mean(ROI1(ROI1 > 0));  
mean_ROI2 = mean(ROI2(ROI2 > 0));  
mean_ROI3 = mean(ROI3(ROI3 > 0));  
mean_ROInoise = mean(ROInoise(ROInoise > 0));  

%%
figure; imagesc(img); colormap gray; axis off; axis tight; axis equal;
hold on;

viscircles([x1, y1], radius, 'Color', 'b', 'LineWidth', 1); 
viscircles([x2, y2], radius, 'Color', 'g', 'LineWidth', 1); 
viscircles([x3, y3], radius, 'Color', 'r', 'LineWidth', 1);  
viscircles([x_noise, y_noise], radius, 'Color', 'y', 'LineWidth', 1); 

% legend({'I rod', 'I solution', 'Gd rod', 'Noise'});
title('ROI Locations on Image');
hold off;

cd (results_file)
print(gcf, 'Gel_Gd_Phantom_202401001_ROIs', '-dpng', '-r300');
cd ../
cd ../

%% Calculate Mean Pixel Intensity
mean_Gd_Rod_117= abs(mean_ROI1); 
mean_Gd_Solution_117 = abs(mean_ROI2);       
mean_Gd_Solution_96 = abs(mean_ROI3);    

mean_values = [mean_Gd_Rod_117,mean_Gd_Solution_117,mean_Gd_Solution_96] ;  

%% Plotting the Mean Pixel Intensity for Each Insert
bar_colors = [0 0.4470 0.7410; 0 0.4470 0.7410;0.8500 0.3250 0.0980]; 

figure;
b = bar(mean_values);

b.FaceColor = 'flat';
b.CData = bar_colors;

set(gca, 'XTickLabel', {'Gd Rod117', 'Gd Solution117', 'Gd Solution96'});
ylabel('Mean Pixel Intensity');
title('Mean Pixel Intensity Inserts')

cd (results_file)
print(gcf, 'Gel_Gd_Phantom_202401001_Mean_Pixel_Intensity', '-dpng', '-r300');
cd ../
cd ../


%% Calculate Contrast
contrast_Gd_Rod_117 = abs(mean_ROI1 - mean_ROInoise) / mean_ROInoise; 
contrast_Gd_Solution_117 = abs(mean_ROI2 - mean_ROInoise) / mean_ROInoise;       
contrast_Gd_96 = abs(mean_ROI3 - mean_ROInoise) / mean_ROInoise;   

contrast_values = [contrast_Gd_Solution_117, contrast_Gd_Rod_117, contrast_Gd_96];

%% Plotting the Contrast for Each Insert
bar_colors = [0 0.4470 0.7410; 0 0.4470 0.7410; 0.8500 0.3250 0.0980]; 

figure;
b = bar(contrast_values);

b.FaceColor = 'flat';
b.CData = bar_colors;

set(gca, 'XTickLabel', {'Gd Rod117', 'Gd Solution117', 'Gd Solution96'});
ylabel('Contrast');
title('Contrast Inserts')

cd (results_file)
print(gcf, 'Gel_Gd_Phantom_202401001_Contrast', '-dpng', '-r300');
cd ../
cd ../

%% SNR 
std_noise = std(ROInoise(ROInoise > 0));  
SNR_Gd_Rod_117 = mean_ROI1 / std_noise;       
SNR_Gd_Solution_117 = mean_ROI2 / std_noise;   
SNR_Gd_96 = mean_ROI3 / std_noise;  

SNR_values = [SNR_Gd_Solution_117, SNR_Gd_Rod_117, SNR_Gd_96];

%% Plotting the SNR for Each Insert
bar_colors = [0 0.4470 0.7410; 0 0.4470 0.7410; 0.8500 0.3250 0.0980];  % Same colors as before

figure;
b = bar(SNR_values);

b.FaceColor = 'flat';
b.CData = bar_colors;

set(gca, 'XTickLabel', {'Gd Rod117', 'Gd Solution117', 'Gd Solution96'});
ylabel('SNR');
title('SNR Inserts');

cd Gel_Phantom/results/
print(gcf, 'Gel_Gd_Phantom_202401001_SNR', '-dpng', '-r300');
cd ../
cd ../