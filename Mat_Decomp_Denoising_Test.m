close all
clear all
mfilename='/home/sriharsha.marupudi/TIGRE-master/MATLAB';
addpath(genpath(mfilename));
current_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023';
addpath(genpath(current_file));
cd (current_file)
current_file = '/home/sriharsha.marupudi/Desktop/PCD/1112023/readingtools/';
addpath(genpath(current_file));
mfilename =  '/home/sriharsha.marupudi/MATLAB Add-Ons/Collections/TV-L1 Image Denoising Algorithm/';
addpath(genpath(mfilename))
results_folder = '/home/sriharsha.marupudi/Desktop/PCD/1112023/Material_Decomposition_Results/';
cd results/
load('2D_reconstruction_material_decomp_SP_25ma_2',"im_low","im_high")
cd ../
load('coordinations_material_decomp_SP_25ma_STC',"y4","x4","y3","x3","y2","x2","y1","x1","yw","xw")


%radial_profile
%% Minimize MSE with Wiener filter 
window = [5,5];
im_low  = wiener2(im_low,window);
im_high = wiener2(im_high,window);

im_low  = imgaussfilt(im_low,window); % Gaussian seems to reduce more noise
im_high = imgaussfilt(im_high,window);

im_low = wavelet_denoising(im_low, 'db4', 3, 'soft', 'universal');

im_high = wavelet_denoising(im_low, 'db4', 3, 'soft', 'universal');

 
figure(1) 
subplot(1,2,1)
    imagesc(im_low); axis tight; axis equal; colormap gray;
subplot(1,2,2)
    imagesc(im_high); axis tight; axis equal; colormap gray;



%%
%voxel_size = 0.1;
Nx = size(im_high,1);
Ny = size(im_high,2);
Nz = size(im_high,3); %2 mm along z
radius_iodine  = 80; %mm diamter adjust this based on measurement
voxel_size = 1;

%% Make this loop 


% water
    Nxh = xw;
    Nyh = yw;
    get_sig_low_high
    sig_wt_low  = sig_low  ;
    sig_wt_high = sig_high ;

    % id 2.5
    Nxh = x1;
    Nyh = y1;
    get_sig_low_high
    sig_id1_low  = sig_low  ;
    sig_id1_high = sig_high ;

    % id 5
    Nxh = x2;
    Nyh = y2;
    get_sig_low_high
    sig_id2_low = sig_low  ;
    sig_id2_high = sig_high ;

    % id 7.5
    Nxh = x3;
    Nyh = y3;
    get_sig_low_high
    sig_id3_low = sig_low  ;
    sig_id3_high = sig_high ;

    % id 10
    Nxh = x4;
    Nyh = y4;
    get_sig_low_high
    sig_id4_low = sig_low  ;
    sig_id4_high= sig_high ;

MU_id_low  = sig_id4_low  - sig_wt_low ; %Subtract out water 
MU_id_high = sig_id4_high - sig_wt_high;

MU_wt_low  =  sig_wt_low ;
MU_wt_high =  sig_wt_high ;


M_R = [MU_wt_low MU_id_low; MU_wt_high MU_id_high];
M_R_inv = inv(M_R)                                                     ;
im_wt    = M_R_inv(1,1)*im_low  + M_R_inv(1,2)*im_high ;
im_id   = M_R_inv(2,1)*im_low  + M_R_inv(2,2)*im_high ;

im_id = im_id*10;% Max mg/ml

%%
close all

subplot(1,2,1)
    imagesc(im_wt); axis tight; axis equal; colormap gray;
    axis off
    title('Water')
    caxis([0 4])
    colorbar
subplot(1,2,2)
    imagesc(im_id); axis tight; axis equal; colormap gray;
    caxis([0 15])
    axis off
    title('Iodine')
    colorbar

 


%% Phantom analysis:


    i = 1
    % water
    Nxh = xw;
    Nyh = yw;
    get_sig
    C(i)  = sig;
    C_std(i) = sig_std;
    i = i + 1;

    % id 0.5
    Nxh = x1;
    Nyh = y1;
    get_sig
    C(i)  = sig;
    C_std(i)  = sig_std;
    i = i + 1;

    % id 2
    Nxh = x2;
    Nyh = y2;
    get_sig
    C(i)  = sig;
    C_std(i)  = sig_std;
    i = i + 1;

    % id 3
    Nxh = x3;
    Nyh = y3;
    get_sig
    C(i)  = sig;
    C_std(i)  = sig_std;
    i = i + 1;

    % id 4
    Nxh = x4;
    Nyh = y4;
    get_sig
    C(i)  = sig;
    C_std(i)  = sig_std;




%%

        C_true = [0 2.5 5 7.5 10];
        figure(10)
        x0=25;
        y0=25;
        width=0.75*550;
        height=0.75*400;
        set(gcf,'position',[x0,y0,width,height])
        hold on
        plot(C_true,C,'--*','linewidth',2);ylim([0,12]);xlim([0,10])
        xlabel('True iodine concentration[mg/mL]')
        ylabel('Estimated iodine concentration[mg/mL]')
         

        C_true_mean = mean(C_true);
        C_mean = mean(C);
        SS_tot = sum((C_true - C_true_mean).^2);
        SS_res = sum((C - C_true).^2);
        r_squared = 1 - (SS_res / SS_tot);

        C_true = cast(C_true, class(C));

        mse_value = immse(C_true, C);
        rmse_value = sqrt(mse_value);
 
        % Display or use the computed r_squared value
        disp(['R^2 value: ', num2str(r_squared)]);
        disp(['Root Mean Squared Error (RMSE) value: ', num2str(rmse_value)]);


        
         