function rebinnedImage = rebinImage(image, newRows, newCols)
    [rows, cols, slices] = size(image);
    
    rowScale = rows / newRows;
    colScale = cols / newCols;
    
    rebinnedImage = zeros(newRows, newCols, slices);
    
    for k = 1:slices
        slice = image(:,:,k);
        
        for i = 1:newRows
            for j = 1:newCols
                rowStart = round((i-1) * rowScale) + 1;
                rowEnd = min(round(i * rowScale), rows);
                colStart = round((j-1) * colScale) + 1;
                colEnd = min(round(j * colScale), cols);
                binRegion = slice(rowStart:rowEnd, colStart:colEnd);
                rebinnedImage(i, j, k) = mean(binRegion(:));
            end
        end
    end
end
