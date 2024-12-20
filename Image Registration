function registered_img = register_image(results_file, coord_file, HU_image, HU_image_2)
    % Navigate to the results directory
    cd(results_file);

    % Normalize images for better visualization in cpselect
    display_gel_1 = mat2gray(HU_image);  % Convert to grayscale for visualization
    display_gel_2 = mat2gray(HU_image_2);

    % Check if control points file exists
    if isfile(coord_file)
        load(coord_file, 'fixedPoints', 'movingPoints');
        disp('Loaded saved control points for Image 2.');
    else
        disp('Selecting new control points for Image 2...');
        % Use the normalized images for control point selection
        [fixedPoints, movingPoints] = cpselect(display_gel_1, display_gel_2, 'Wait', true);
        save(coord_file, 'fixedPoints', 'movingPoints');
        disp('Control points for Image 2 saved.');
    end

    % Compute the transformation
    tform_2 = fitgeotrans(movingPoints, fixedPoints, 'similarity');
    output_view_2 = imref2d(size(HU_image));  % Reference for the 3D dataset

    % Apply the transformation to register the second image
    reg
