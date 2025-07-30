close all
clear all
clc
restoredefaultpath
mfilename = '/home/sriharsha.marupudi/TIGRE-master/MATLAB';
addpath(genpath(mfilename));
current_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/';
addpath(genpath(current_file));
mfilename = 'home/sriharsha.marupudi/Desktop/PCD/1112023/Gel_Phantom/readingtools/';
addpath(genpath(mfilename));
addpath '/home/sriharsha.marupudi/Desktop/PCD/1112023'
results_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/Gel_Phantom/results/NPS';


%%
load('data_Iodine_Container_20241218_1_Recon_20241218',"img_low");
im_low_1 = img_low; 
im_low_1 = img_low(:,:,50:60);
load('data_Iodine_Container_20241218_2_Recon_20241218',"img_low");
im_low_2 = img_low; 
im_low_2 = img_low(:,:,50:60); 
load('data_Iodine_Container_20241218_3_Recon_20241218',"img_low");
im_low_3 = img_low; 
im_low_3 = img_low(:,:,50:60); 



% load('data_Iodine_Container_20250113_1_Recon_20250113',"im_low");
% im_low_1 = im_low; 
% load('data_Iodine_Container_20250113_2_Recon_20250113',"im_low");
% im_low_2 = im_low; 
% load('data_Iodine_Container_20250113_3_Recon_20250113',"im_low");
% im_low_3 = im_low; 
% % 
% 
% load('data_Iodine_Container_Water_Wax_Recon_1_20241218_2',"im_low");
% im_low_1 = im_low; 
% load('data_Iodine_Container_Water_Wax_Recon_2_20241218_2',"im_low");
% im_low_2 = im_low; 
% load('data_Iodine_Container_Water_Wax_Recon_3_20241218_2',"im_low");
% im_low_3 = im_low; 
% 


%%
b_data = cell(1,10);
a_data = cell(1,10);
index = 1; 
kernel = ones(3,3)/9;  

im_low_1 = im_low_1 - mean(im_low_1, "all");  
im_low_2 = im_low_2 - mean(im_low_2, "all");
im_low_3 = im_low_3 - mean(im_low_3, "all");



% cutoff_frequency = 3;
% filter_size =10;
% gaussian_filter = fspecial('gaussian', filter_size, cutoff_frequency);
% im_low_1 = im_low_1 - imfilter(im_low_1, gaussian_filter, 'same', 'replicate');
% im_low_2 = im_low_2 - imfilter(im_low_2, gaussian_filter, 'same', 'replicate');
% im_low_3 = im_low_3 - imfilter(im_low_3, gaussian_filter, 'same', 'replicate');


%% Calculate NPS 

[f_im_low_1,b_avg_im_low_1] = NPS_function2(im_low_1);
[f_im_low_2,b_avg_im_low_2] = NPS_function2(im_low_2);
[f_im_low_3,b_avg_im_low_3] = NPS_function2(im_low_3);


%% Apodization 
window = hanning(length(b_avg_im_low_1));
% window = hamming(length(b_avg_im_low_1));   


% Apply the window to each NPS data
b_avg_im_low_1_windowed = b_avg_im_low_1 .* window';
b_avg_im_low_2_windowed = b_avg_im_low_2 .* window';
b_avg_im_low_3_windowed = b_avg_im_low_3 .* window';


NPS_average_windowed = (b_avg_im_low_1_windowed + b_avg_im_low_2_windowed + b_avg_im_low_3_windowed) / 3;
f_im_low_avg = (f_im_low_1 + f_im_low_2 + f_im_low_3) / 3;
save(fullfile(results_file, 'Average_NPS_Iodine_Container_20241218.mat'), 'NPS_average_windowed', 'f_im_low_avg');

figure;
plot(f_im_low_avg, NPS_average_windowed, '-*', 'LineWidth', 1.5);
title('Average NPS');
xlabel('Frequency (cycle/mm)');
ylabel('Average NPS');

cd(results_file);
print(gcf, 'Average_NPS_Iodine_Container_20241218', '-dpng', '-r300');
cd ../;
cd ../;