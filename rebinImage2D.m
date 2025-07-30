function rebinnedImage = rebinImage2D(image, newRows, newCols)
    [rows, cols] = size(image);
    
    rowScale = rows / newRows;
    colScale = cols / newCols;
    
    rebinnedImage = zeros(newRows, newCols);
    
    for i = 1:newRows
        for j = 1:newCols
            rowStart = round((i - 1) * rowScale) + 1;
            rowEnd   = min(round(i * rowScale), rows);
            colStart = round((j - 1) * colScale) + 1;
            colEnd   = min(round(j * colScale), cols);
            
            binRegion = image(rowStart:rowEnd, colStart:colEnd);
            rebinnedImage(i, j) = mean(binRegion(:));
        end
    end
end
