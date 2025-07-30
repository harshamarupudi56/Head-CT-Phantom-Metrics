function [f, b_avg, a_avg, NNPS] = compute_NPS_and_plot(image)
    index = 1; 
    all_b_data = [];  % Store all b values for averaging
    all_a_data = [];  % Store all a values for averaging
    all_LAS = [];     % Store LAS values for normalization

    % Defining Pixel, voxel size, and the corresponding f according to Nyquist Criterion
    DeltaX = 0.1;  % Voxel size of x and y 100 microns 
    f = (1:256/2)/256/DeltaX;

    % Loop through slices independently
    for idx = 1:10
        DI = image(:,:,idx);  % Process each slice independently
        
        % Display image with ROI locations
        figure(1);
        imshow(DI, []);
        title(['Selected ROIs on Slice ', num2str(idx)]);
        hold on;
        
        % Extract 16 Overlapped ROIs
        for row = 1170:50:1220  % Corrected range
            for col = 1032:50:1082  % Corrected range
                rect = [col, row, 255, 255];
                
                % Draw rectangle to indicate ROI location
                rectangle('Position', rect, 'EdgeColor', 'r', 'LineWidth', 1.5);
                
                % Extract ROI
                DX{row, col} = imcrop(DI, rect);
                Detrending_ROI = DX{row, col} - mean2(DX{row, col});  % Detrend the ROI
                
                % Compute NPS
                NPS{row, col} = fft2(Detrending_ROI);
                abs_NPS{row, col} = (abs(NPS{row, col}).^2);
                shifted_NPS{row, col} = fftshift(abs_NPS{row, col});
                NPS_2D{row, col} = DeltaX^2 * shifted_NPS{row, col} / 256^2;
                
                % Compute NPS value for a region
                a{row, col} = mean(mean(NPS_2D{row, col}(128:148, 117:137)));  
            end
        end
        hold off; 
        
        % Compute the average NPS for the slice across all ROIs
        a_values = []; % Initialize array to store 'a' values
        for row = 1170:50:1220
            for col = 1032:50:1082
                if isfield(a, sprintf('%d,%d', row, col)) % Ensure field exists
                    a_values = [a_values, a{row, col}]; % Extract values into array
                end
            end
        end
        a_new = mean(a_values, 'all'); % Compute mean safely
    
        % Aggregate full 2D NPS for this slice
        NPS_sum = zeros(size(NPS_2D{1170,1032}));
        for row = 1170:50:1220
            for col = 1032:50:1082
                NPS_sum = NPS_sum + NPS_2D{row, col};
            end
        end
        a_new2 = NPS_sum / 16;  % Average across ROIs
        
        % Radial NPS measurement
        b = rscan(a_new2);
        b(1, 128) = b(1, 127);  % Fix potential edge issue
        
        % Store NPS values for averaging later
        all_b_data = [all_b_data; b];
        all_a_data = [all_a_data; a_new];
        all_LAS = [all_LAS; mean2(DI)];  % Store mean signal intensity for normalization

        % Plot individual slice NPS
        % figure(3);
        % plot(f, a_new);
        % hold on;
        % 
        % Store slice-wise NPS data
        b_data{index} = b;
        a_data{index} = a_new;
        index = index + 1;
    end
    hold off;

    % Compute final averaged NPS
    b_avg = mean(all_b_data, 1);
    a_avg = mean(all_a_data, 1);
    
    % Compute NNPS (normalized NPS)
    LAS_avg = mean(all_LAS);  % Properly averaged over slices
    NNPS = b_avg / LAS_avg^2;  
    
end

