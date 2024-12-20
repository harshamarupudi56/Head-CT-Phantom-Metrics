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
