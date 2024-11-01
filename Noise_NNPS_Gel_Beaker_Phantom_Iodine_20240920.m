%% Noise and NPS Analysis 

close all
clear all
clc
restoredefaultpath
mfilename = '/home/sriharsha.marupudi/TIGRE-master/MATLAB';
addpath(genpath(mfilename));
mfilename = 'home/sriharsha.marupudi/Desktop/PCD/1112023/Gel_Phantom/readingtools/';
addpath(genpath(mfilename));
use_par = 1; %use parralelisation
addpath '/home/sriharsha.marupudi/Desktop/PCD/1112023'
results_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/Gel_Phantom/results/';
current_file_MTF = '/home/sriharsha.marupudi/Desktop/PCD/1112023/3DMTF/MTFTools-main/examples';
cd(current_file_MTF);
run Setup; 
current_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/' ;
addpath(genpath(current_file));
cd(current_file); 

cd Gel_Phantom/results/
load('data_Gel_Beaker_Phantom_Iodine__Recon_20240920.mat',"img_low")
cd ../ 

%% 

img_low_noise = img_low - mean(img_low, "all");
figure; imagesc(img_low(:,:,110)); colormap gray; axis off; axis tight; axis equal; 
figure; imagesc(img_low_noise(:,:,110)); colormap gray; axis off; axis tight; axis equal; 

% Apply high pass filter (adjust based on your existing filter settings)
cutoff_frequency = 6;
filter_size = 25;
gaussian_filter = fspecial('gaussian', filter_size, cutoff_frequency);
img_low_filt = img_low_noise - imfilter(img_low_noise, gaussian_filter, 'same', 'replicate');
[f_low, b_avg_low, a_avg_low, NNPS_low] = NPS_function1(img_low_filt);

plot(f_low,NNPS_low)
xlabel('Spatial Frequency');
ylabel('NNPS');