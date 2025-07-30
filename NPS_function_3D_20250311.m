function [bin_centers, radial_avg] = compute_3D_NPS_with_padding(input_volume)

img_stack = input_volume;  
[Nx, Ny, Nz] = size(img_stack);

% Define desired VOI size
VOI_size = [256, 256, 256];  
DeltaX = 0.1; DeltaY = 0.1; DeltaZ = 0.1;
scaling_factor = (DeltaX / 256)^3;

% Define ROI coordinates for VOI placement (change coordinates as needed)
ROI_coordinates_2D = [
    1250, 1250; 1250, 1300; 1250, 1300; 1250, 1350;
    1300, 1250; 1300, 1300; 1300, 1300; 1300, 1350;
    1300, 1250; 1300, 1300; 1300, 1300; 1300, 1350;
    1350, 1250; 1350, 1300; 1350, 1300; 1350, 1350
];

% Initialize the average NPS matrix
average_NPS = zeros(VOI_size(1), VOI_size(2));
overlay_image = mean(img_stack, 3);  

% Process each VOI
for idx = 1:size(ROI_coordinates_2D, 1)
    x_start = ROI_coordinates_2D(idx, 1);
    y_start = ROI_coordinates_2D(idx, 2);
    z_start = 1;
    z_end = Nz;  % Image's Z-dimension ends at 50

    % Extract the VOI from the image volume
    x_end = x_start + VOI_size(1) - 1;
    y_end = y_start + VOI_size(2) - 1;

    voi = img_stack(x_start:min(x_end, Nx), y_start:min(y_end, Ny), z_start:z_end);

    % Pad VOI to desired size (256 x 256 x 256) using symmetric padding
    padded_voi = padarray(voi, [VOI_size(1) - size(voi, 1), VOI_size(2) - size(voi, 2), VOI_size(3) - size(voi, 3)], 'symmetric', 'post');

    % Detrend the VOI 
    detrended_voi = padded_voi;
    % for z = 1:size(padded_voi, 3)
    %     voi_slice = padded_voi(:, :, z);
    %     detrended_voi(:, :, z) = voi_slice - mean2(voi_slice);  
    % end

    % Initialize NPS for this VOI
    nps = zeros(VOI_size(1), VOI_size(2), VOI_size(3));
    for z = 1:VOI_size(3)
        voi_slice = detrended_voi(:, :, z);

        % Perform FFT on the slice
        fftData = fftshift(fftn(voi_slice));
        d_value = abs(fftData);
        nps(:, :, z) = (((d_value.^2) / 2)./16) * scaling_factor;
    end

    % Integrate NPS across the Z-dimension
    integrated_NPS = sum(nps, 3);

    % Accumulate the results into the average NPS
    average_NPS = average_NPS + integrated_NPS;
end

% Average the NPS across all VOIs
average_NPS = average_NPS / size(ROI_coordinates_2D, 1);

% Radial Frequency Calculation
[fx, fy] = meshgrid(...
    (-VOI_size(2)/2:VOI_size(2)/2-1) / (VOI_size(2) * DeltaX), ...
    (-VOI_size(1)/2:VOI_size(1)/2-1) / (VOI_size(1) * DeltaY));
r = sqrt(fx.^2 + fy.^2);  

% Flatten radial frequency and NPS data
r_flat = r(:);
nps_flat = average_NPS(:);

% Define radial bins
n_bins = 100; 
max_r = max(r_flat);  
bin_edges = linspace(0, max_r, n_bins + 1);

% Compute radial averages
bin_centers = zeros(1, n_bins);
radial_avg = zeros(1, n_bins);
for b = 1:n_bins
    bin_mask = (r_flat >= bin_edges(b)) & (r_flat < bin_edges(b + 1));
    bin_centers(b) = (bin_edges(b) + bin_edges(b + 1)) / 2;
    radial_avg(b) = mean(nps_flat(bin_mask), 'omitnan'); % Exclude NaN values
end
