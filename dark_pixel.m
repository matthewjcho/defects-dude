% - Rows
% - Columns
% - Clusters (partial rows and partial columns - rows and columns that do not streak across the entire image, and chunks - pixels that are clumped side-by-side)
% - Pixels

clearvars;
clc;

% Load Image
imageFolder = 'C:\Users\matthewjcho\OneDrive\MATLAB Scripts & Development\Test';
imageFile = dir(fullfile(imageFolder, '*.tif'));
I = imread(fullfile(imageFile.folder, imageFile.name));

% Crop image to ROI (consistent for all 8K x 8K images, but adjust if necessary)
croppedX_start = 1; croppedX_end = 4096;
croppedY_start = 1; croppedY_end = 4096;
I = I(croppedY_start:croppedY_end, croppedX_start:croppedX_end, :);

% Convert to double for accurate calculations
I = double(I);  

% Calculate mean pixel value and define bright defect threshold
meanValue = mean(I, 'all');
criterionDark = meanValue * 0.5;

% Identify bright defects using logical indexing
darkDefect = I < criterionDark;

% Extract defect locations
[defectRows, defectCols] = find(darkDefect);

% Offset cropped coordinates back to original image dimensions
originalCols = defectCols + (croppedX_start - 1);
originalRows = defectRows + (croppedY_start - 1);

% Convert to ImageJ coordinates (zero-based)
imageJCols = originalCols - 1;
imageJRows = originalRows - 1;
imageJOriginal = [imageJCols, imageJRows];

rowRanges = [];
colRanges = [];
partialRowRanges = [];
partialColRanges = [];
removePartialRows = [];
removePartialCols = [];
duplicatePartialYCount = [];


if ~isempty(imageJOriginal)
    sortedPointsX = sortrows(imageJOriginal, [1,2]);

    X2 = sortedPointsX(:,1);
    Y2 = sortedPointsX(:,2);

    potentialColCoordinatesIdx = find(diff(X2) == 0);
    if ~isempty(potentialColCoordinatesIdx)
        potentialColCoordinatesIdx = unique([potentialColCoordinatesIdx; potentialColCoordinatesIdx + 1]);
    end

    potentialColCoordinates = [X2(potentialColCoordinatesIdx), Y2(potentialColCoordinatesIdx)];

    if ~isempty(potentialColCoordinates)

        colX1 = potentialColCoordinates(:,1);
        colY1 = potentialColCoordinates(:,2);
        
        xSame = diff(colX1) == 0;
        col_changeIdx = find(xSame);
    
        if ~isempty(col_changeIdx)
            col_changeIdx = unique([1; col_changeIdx; col_changeIdx + 1]);
        end
        
        col_xTrue = colX1(col_changeIdx);
        col_yTrue = colY1(col_changeIdx);
        
        duplicateX = find(diff([col_xTrue; numel(col_xTrue)]) ~= 0);
        col_segmentLengths = diff([1; duplicateX]); % Count repeated X sequences
        
        col_validLengthIdx = find(col_segmentLengths > 100);
        
        validXValues = col_xTrue(duplicateX(col_validLengthIdx));
        
        matchingXIdx = ismember(colX1, validXValues);
        
        col_validX = colX1(matchingXIdx);
        col_validY = colY1(matchingXIdx);
        
        removePotentialCols = [col_validX, col_validY];
        removeAllCols = removePotentialCols;
    
        if ~isempty(removeAllCols)
    
            colX = removeAllCols(:,1);
            colY = removeAllCols(:,2);
    
            col_rangeIdx = find(diff(colX) >= 1);
    
            if ~isempty(col_rangeIdx)
    
                col_firstSegmentX = colX(1:col_rangeIdx(1));
                col_firstSegmentY = colY(1:col_rangeIdx(1));
                
                col_firstSegmentXStart = col_firstSegmentX(1);
                col_firstSegmentYStart = col_firstSegmentY(1);
                col_firstSegmentXEnd = col_firstSegmentX(end);
                col_firstSegmentYEnd = col_firstSegmentY(end);
        
                col_lastSegmentX = colX(col_rangeIdx(end)+1:end);
                col_lastSegmentY = colY(col_rangeIdx(end)+1:end);
        
                col_lastSegmentXStart = col_lastSegmentX(1);
                col_lastSegmentYStart = col_lastSegmentY(1);
                col_lastSegmentXEnd = col_lastSegmentX(end);
                col_lastSegmentYEnd = col_lastSegmentY(end);
        
                col_firstSegment = [col_firstSegmentXStart, col_firstSegmentYStart; col_firstSegmentXEnd, col_firstSegmentYEnd];
                col_lastSegment = [col_lastSegmentXStart, col_lastSegmentYStart; col_lastSegmentXEnd, col_lastSegmentYEnd];
        
                for d = 1:numel(col_rangeIdx) - 1
                    col_segmentX = colX(col_rangeIdx(d)+1:col_rangeIdx(d+1));
                    col_segmentY = colY(col_rangeIdx(d)+1:col_rangeIdx(d+1));
                    
                    col_segmentXStart = col_segmentX(1);
                    col_segmentYStart = col_segmentY(1);
                    col_segmentXEnd = col_segmentX(end);
                    col_segmentYEnd = col_segmentY(end);
        
                    colRanges = [colRanges; col_segmentXStart, col_segmentYStart; col_segmentXEnd, col_segmentYEnd];
                
                end
    
                colRanges = [col_firstSegment; colRanges; col_lastSegment];
    
            else
                col_segmentXStart = colX(1);
                col_segmentYStart = colY(1);
                col_segmentXEnd = colX(end);
                col_segmentYEnd = colY(end);
    
                colRanges = [colRanges; col_segmentXStart, col_segmentYStart; col_segmentXEnd, col_segmentYEnd];
    
            end          
        end
    
    
    
        if ~isempty(removeAllCols)
            %colsTable = array2table(removeAllCols, 'VariableNames', {'X', 'Y'});
            %writetable(colsTable, 'dark defects.xlsx', 'Sheet', 'Col Coordinates', 'WriteMode', 'Append')
    
            col_toRemove = ismember(imageJOriginal, removeAllCols, 'rows');
            imageJOriginal = imageJOriginal(~col_toRemove, :);
    
            if ~isempty(colRanges)
                numCols = numel(colRanges(:,1)) / 2;
    
                if numel(colRanges(:,1)) == 2
                    numCols = 1;
                end
    
                colRanges(1,3) = numCols;
                colRanges(2:end, 3) = NaN;
                colRanges(any(colRanges(:,3) == 0 & (1:size(colRanges,1))' ~= 1, 2), :) = [];
                
                colRangesTable = array2table(colRanges, 'VariableNames', {'Range X', 'Range Y', 'Count'});
                writetable(colRangesTable, 'dark defects.xlsx', 'Sheet', 'Column Ranges', 'WriteMode', 'Append');
            end
        end
    end

    sortedPointsY = sortrows(imageJOriginal, [2,1]);

    X1 = sortedPointsY(:,1);
    Y1 = sortedPointsY(:,2);

    potentialRowCoordinatesIdx = find(diff(Y1) == 0);
    if ~isempty(potentialRowCoordinatesIdx)
        potentialRowCoordinatesIdx = unique([potentialRowCoordinatesIdx; potentialRowCoordinatesIdx + 1]);
    end

    potentialRowCoordinates = [X1(potentialRowCoordinatesIdx), Y1(potentialRowCoordinatesIdx)];
    
    if ~isempty(potentialRowCoordinates)        
        rowX1 = potentialRowCoordinates(:,1);
        rowY1 = potentialRowCoordinates(:,2);
        
        ySame = diff(rowY1) == 0;
        changeIdx = find(ySame);
        
        xTrue = rowX1(changeIdx);
        yTrue = rowY1(changeIdx);
        
        duplicateY = find(diff([yTrue; numel(yTrue)]) ~= 0);
        segmentLengths = diff([1; duplicateY]); % Count repeated Y sequences
        
        validLengthIdx = find(segmentLengths > 100);
        
        validYValues = yTrue(duplicateY(validLengthIdx));
        
        matchingYIdx = ismember(rowY1, validYValues);
        
        validX = rowX1(matchingYIdx);
        validY = rowY1(matchingYIdx);
        
        removePotentialRows = [validX, validY];
        removeAllRows = removePotentialRows;
            
        if ~isempty(removeAllRows)
    
            rowX = removeAllRows(:,1);
            rowY = removeAllRows(:,2);
    
            rangeIdx = find(diff(rowY) >= 1);
    
            if ~isempty(rangeIdx)
            
                firstSegmentX = rowX(1:rangeIdx(1));
                firstSegmentY = rowY(1:rangeIdx(1));
                
                firstSegmentXStart = firstSegmentX(1);
                firstSegmentYStart = firstSegmentY(1);
                firstSegmentXEnd = firstSegmentX(end);
                firstSegmentYEnd = firstSegmentY(end);
        
                lastSegmentX = rowX(rangeIdx(end)+1:end);
                lastSegmentY = rowY(rangeIdx(end)+1:end);
        
                lastSegmentXStart = lastSegmentX(1);
                lastSegmentYStart = lastSegmentY(1);
                lastSegmentXEnd = lastSegmentX(end);
                lastSegmentYEnd = lastSegmentY(end);
        
                firstSegment = [firstSegmentXStart, firstSegmentYStart; firstSegmentXEnd, firstSegmentYEnd];
                lastSegment = [lastSegmentXStart, lastSegmentYStart; lastSegmentXEnd, lastSegmentYEnd];
        
                for j = 1:numel(rangeIdx) - 1
                    segmentX = rowX(rangeIdx(j)+1:rangeIdx(j+1));
                    segmentY = rowY(rangeIdx(j)+1:rangeIdx(j+1));
                    
                    segmentXStart = segmentX(1);
                    segmentXEnd = segmentX(end);
                    segmentYStart = segmentY(1);
                    segmentYEnd = segmentY(end);
        
                    rowRanges = [rowRanges; segmentXStart, segmentYStart; segmentXEnd, segmentYEnd];
                    
                end
        
                
                rowRanges = [firstSegment; rowRanges; lastSegment];
    
            else
                segmentXStart = rowX(1);
                segmentYStart = rowY(1);
                segmentXEnd = rowX(end);
                segmentYEnd = rowY(end);
    
                rowRanges = [rowRanges; segmentXStart, segmentYStart; segmentXEnd, segmentYEnd];
    
            end       
        end
    
    
        if ~isempty(removeAllRows)
            %rowsTable = array2table(removeAllRows, 'VariableNames', {'X', 'Y'});
            %writetable(rowsTable, 'dark defects.xlsx', 'Sheet', 'Row Coordinates');
    
            row_toRemove = ismember(imageJOriginal, removeAllRows, 'rows');
            imageJOriginal = imageJOriginal(~row_toRemove, :);
    
            if ~isempty(rowRanges)
                numRows = numel(rowRanges(:,1)) / 2;
    
                if numel(rowRanges(:,1)) == 2
                    numRows = 1;
                end
    
                rowRanges(1,3) = numRows;
                rowRanges(2:end, 3) = NaN;
                rowRanges(any(rowRanges(:,3) == 0 & (1:size(rowRanges,1))' ~= 1, 2), :) = [];
                
                rowRangesTable = array2table(rowRanges, 'VariableNames', {'Range X', 'Range Y', 'Count'});
                writetable(rowRangesTable, 'dark defects.xlsx', 'Sheet', 'Row Ranges', 'WriteMode', 'Append');
            end       
        end
    end

    
    % At this point, ImageJ has no rows and columns

    if ~isempty(imageJOriginal)
        potentialClusterCoordinates2 = [imageJOriginal(:,1), imageJOriginal(:,2)];

        sorted_potentialClusterX = sortrows(potentialClusterCoordinates2, [1,2]);
        
        x_sortedPotentialClusterX = sorted_potentialClusterX(:,1);
        y_sortedPotentialClusterX = sorted_potentialClusterX(:,2);

        findDuplicatePartialX = find(diff(x_sortedPotentialClusterX) == 0);

        if ~isempty(findDuplicatePartialX)
            findDuplicatePartialX = unique([1; findDuplicatePartialX; findDuplicatePartialX + 1]);
    
            duplicatePartialX = x_sortedPotentialClusterX(findDuplicatePartialX);
        
            findPartialColX = ismember(x_sortedPotentialClusterX, duplicatePartialX, 'rows');
        
            potentialPartialColX = x_sortedPotentialClusterX(findPartialColX);
            potentialPartialColY = y_sortedPotentialClusterX(findPartialColX);
        
            potentialPartialColY_gap = find(diff(potentialPartialColX) >= 1);
    
            if ~isempty(potentialPartialColY_gap)
                potentialPartialColY_gap = unique([1; potentialPartialColY_gap; potentialPartialColY_gap + 1; numel(potentialPartialColY)]);
                
                countY = diff(potentialPartialColY_gap) + 1;
                findCountY = find(countY > 4);
            
                if ~isempty(findCountY)
                    findCountY = unique([findCountY; findCountY + 1]);
                end

                findFilteredPartialColValue = potentialPartialColX(potentialPartialColY_gap(findCountY)); % Here is the issue
                
                findFilteredPartialCol = ismember(potentialPartialColX, findFilteredPartialColValue);

                filteredPartialColX = potentialPartialColX(findFilteredPartialCol);
                filteredPartialColY = potentialPartialColY(findFilteredPartialCol);

                incrementYIdx = find(abs(diff(filteredPartialColY)) <= 15);
                
                if ~isempty(incrementYIdx)
                    incrementYIdx = unique([incrementYIdx; incrementYIdx + 1]);
                end

                validPartialColX = filteredPartialColX(incrementYIdx);
                validPartialColY = filteredPartialColY(incrementYIdx);

                validXSame = find(diff(validPartialColX) >= 1);

                if ~isempty(validXSame)
                    validXSame = unique([1; validXSame; validXSame + 1; numel(incrementYIdx)]);

                    occurrencesX = diff(validXSame) + 1;
                    findOccurrencesX = find(occurrencesX > 4);
                    
                    if ~isempty(findOccurrencesX)
                        findOccurrencesX = unique([findOccurrencesX; findOccurrencesX + 1]);
                    end

                    partialColXValue = validPartialColX(validXSame(findOccurrencesX));
                    
                    findPartialColXValue = ismember(validPartialColX, partialColXValue, 'rows');

                    partialColX = validPartialColX(findPartialColXValue);
                    partialColY = validPartialColY(findPartialColXValue);

                    partialColCoordinates = [partialColX, partialColY];

                else
                    gapY = find(diff(validPartialColY) > 50);
                    
                    if numel(validPartialColY) > 4
                        if ~isempty(gapY)
                            partialColCoordinates = [];
                        else
                            partialColCoordinates = [validPartialColX, validPartialColY];
    
                        end

                    else
                        partialColCoordinates = [];
                    end
                end


                if ~isempty(partialColCoordinates)
                    partialCol_toRemove = ismember(imageJOriginal, partialColCoordinates, 'rows');
                    imageJOriginal = imageJOriginal(~partialCol_toRemove, :);

                    findPartialCols = find(diff(partialColCoordinates(:,1)) >= 1);
                    findPartialCols = [findPartialCols; numel(partialColCoordinates(:,1))];

                    numPartialCols = numel(findPartialCols);
                    partialColCoordinates(1,3) = numPartialCols;
                    partialColCoordinates(2:end, 3) = NaN;
                    partialColCoordinates(any(partialColCoordinates(:,3) == 0 & (1:size(partialColCoordinates,1))' ~= 1, 2), :) = [];
                        
                    colClustersTable = array2table(partialColCoordinates, 'VariableNames', {'X', 'Y', 'Count'});
                    writetable(colClustersTable, 'dark defects.xlsx', 'Sheet', 'Partial Column Clusters', 'WriteMode', 'Append');

                end
                
           end
        end
    end

    if ~isempty(imageJOriginal)
        potentialClusterCoordinates1 = [imageJOriginal(:,1), imageJOriginal(:,2)];
    
        sorted_potentialClusterY = sortrows(potentialClusterCoordinates1, [2,1]);
    
        x_sortedPotentialClusterY = sorted_potentialClusterY(:,1);
        y_sortedPotentialClusterY = sorted_potentialClusterY(:,2);
        
        findDuplicatePartialY = find(diff(y_sortedPotentialClusterY) == 0);

        if ~isempty(findDuplicatePartialY)
            findDuplicatePartialY = unique([1; findDuplicatePartialY; findDuplicatePartialY + 1]);
    
            duplicatePartialY = y_sortedPotentialClusterY(findDuplicatePartialY);
        
            findPartialRowY = ismember(y_sortedPotentialClusterY, duplicatePartialY, 'rows');
        
            potentialPartialRowX = x_sortedPotentialClusterY(findPartialRowY);
            potentialPartialRowY = y_sortedPotentialClusterY(findPartialRowY);
        
            potentialPartialRowX_gap = find(diff(potentialPartialRowY) >= 1);
    
            if ~isempty(potentialPartialRowX_gap)
                potentialPartialRowX_gap = unique([1; potentialPartialRowX_gap; potentialPartialRowX_gap + 1; numel(potentialPartialRowX)]);
            
                countX = diff(potentialPartialRowX_gap) + 1;
                findCountX = find(countX > 4);

                if ~isempty(findCountX)
                    findCountX = unique([findCountX; findCountX + 1]);
                end

                findFilteredPartialRowValue = potentialPartialRowY(potentialPartialRowX_gap(findCountX));
                
                findFilteredPartialRow = ismember(potentialPartialRowY, findFilteredPartialRowValue, 'rows');

                filteredPartialRowX = potentialPartialRowX(findFilteredPartialRow);
                filteredPartialRowY = potentialPartialRowY(findFilteredPartialRow);

                incrementXIdx = find(abs(diff(filteredPartialRowX)) <= 15);
                
                if ~isempty(incrementXIdx)
                    incrementXIdx = unique([incrementXIdx; incrementXIdx + 1]);
                end
                
                validPartialRowX = filteredPartialRowX(incrementXIdx);
                validPartialRowY = filteredPartialRowY(incrementXIdx);

                validYSame = find(diff(validPartialRowY) >= 1);
                if ~isempty(validYSame)
                    validYSame = unique([1; validYSame; validYSame + 1; numel(incrementXIdx)]);

                    occurrencesY = diff(validYSame) + 1;
                    findOccurrencesY = find(occurrencesY > 4);

                    if ~isempty(findOccurrencesY)
                        findOccurrencesY = unique([findOccurrencesY; findOccurrencesY + 1]);
                    end

                    partialRowYValue = validPartialRowY(validYSame(findOccurrencesY));
                    
                    findPartialRowYValue = ismember(validPartialRowY, partialRowYValue, 'rows');

                    partialRowX = validPartialRowX(findPartialRowYValue);
                    partialRowY = validPartialRowY(findPartialRowYValue);

                    partialRowCoordinates = [partialRowX, partialRowY];

                else
                    gapX = find(diff(validPartialRowX) > 50);
                    
                    if numel(validPartialRowX) > 4
                        if ~isempty(gapX)
                            partialRowCoordinates = [];
                        else
                            partialRowCoordinates = [validPartialRowX, validPartialRowY];
    
                        end

                    else
                        partialRowCoordinates = [];
                    end
                end
                

                if ~isempty(partialRowCoordinates)
                    partialRow_toRemove = ismember(imageJOriginal, partialRowCoordinates, 'rows');
                    imageJOriginal = imageJOriginal(~partialRow_toRemove, :);

                    findPartialRows = find(diff(partialRowCoordinates(:,2)) >= 1);
                    findPartialRows = [findPartialRows; numel(partialRowCoordinates(:,2))];

                    numPartialRows = numel(findPartialRows);
                    partialRowCoordinates(1,3) = numPartialRows;
                    partialRowCoordinates(2:end, 3) = NaN;
                    partialRowCoordinates(any(partialRowCoordinates(:,3) == 0 & (1:size(partialRowCoordinates,1))' ~= 1, 2), :) = [];

                    rowClustersTable = array2table(partialRowCoordinates, 'VariableNames', {'X', 'Y', 'Count'});
                    writetable(rowClustersTable, 'dark defects.xlsx', 'Sheet', 'Partial Row Clusters', 'WriteMode', 'Append');
                end
                
            end
        end
    end
  
   
    if ~isempty(imageJOriginal)
        pixelCoordinates = imageJOriginal;
        numPixels = numel(pixelCoordinates(:,1));
        pixelCoordinates(1,3) = numPixels;
        pixelCoordinates(2:end, 3) = NaN;
        pixelCoordinates(any(pixelCoordinates(:,3) == 0 & (1:size(pixelCoordinates,1))' ~= 1, 2), :) = [];
        pixelsTable = array2table(pixelCoordinates, 'VariableNames', {'X', 'Y', 'Count'});
        writetable(pixelsTable, 'dark defects.xlsx', 'Sheet', 'Pixels', 'WriteMode', 'Append');
        
        % Chunks
        sortedPixelCoordinatesX = sortrows(imageJOriginal, [1,2]);
        sortedPixelCoordinatesY = sortrows(imageJOriginal, [2,1]);
    
        x_sortedPixelCoordinatesX = sortedPixelCoordinatesX(:,1);
        y_sortedPixelCoordinatesX = sortedPixelCoordinatesX(:,2);
    
        x_sortedPixelCoordinatesY = sortedPixelCoordinatesY(:,1);
        y_sortedPixelCoordinatesY = sortedPixelCoordinatesY(:,2);
    
        findChunksXGap = find(diff(x_sortedPixelCoordinatesX) == 0);
        
        if ~isempty(findChunksXGap)
            findChunksXGap = unique([1; findChunksXGap; findChunksXGap + 1]);

            potentialChunkX = x_sortedPixelCoordinatesX(findChunksXGap);
            potentialChunkY = y_sortedPixelCoordinatesX(findChunksXGap);
        
            findChunksYGap = find(diff(potentialChunkY) == 1);
            findChunksYGap = unique([1; findChunksYGap; findChunksYGap + 1]);
        
            chunkX = potentialChunkX(findChunksYGap);
            chunkY = potentialChunkY(findChunksYGap);
        
            findChunkYValues = find(diff(chunkY) ~= 0);
            findChunkYValues = unique([1; findChunkYValues; findChunkYValues + 1]);
            
            validChunkX = chunkX(findChunkYValues);
            validChunkY = chunkY(findChunkYValues);
        
            countChunkY = arrayfun(@(x) sum(validChunkY == x), validChunkY);
        
            chunkYMask = countChunkY >= 2 & countChunkY < 10; 
        
            filteredChunkX = validChunkX(chunkYMask);
            filteredChunkY = validChunkY(chunkYMask);
        
            findChunkXValues = find(diff(filteredChunkX) <= 1);

            if ~isempty(findChunkXValues)
                findChunkXValues = unique([1; findChunkXValues; findChunkXValues + 1]);
            end
        
            trueChunkX = filteredChunkX(findChunkXValues);
            trueChunkY = filteredChunkY(findChunkXValues);
            trueChunkCoordinates = [trueChunkX, trueChunkY];

            if ~isempty(trueChunkCoordinates)
                findChunkXValues2 = find(diff(trueChunkX) == 1);
                pairedIdx = [findChunkXValues2; findChunkXValues2 + 1];
                pairedIdx = sort(pairedIdx);
            
                if mod(numel(pairedIdx), 2) ~= 0
                    error('pairedIdx must have an even number of elements.');
                end
            
                pairedIdx = reshape(pairedIdx, 2, []);
            
                pairedYValues = [trueChunkY(pairedIdx(1,:)), trueChunkY(pairedIdx(2,:))];
                validPairMask = (abs(pairedYValues(:,1) - pairedYValues(:,2)) <= 10) & (abs(pairedYValues(:,1) - pairedYValues(:,2)) > 0);
                pairedIdx = transpose(pairedIdx);
                validPairs = pairedIdx(validPairMask, :);
                
                rangeOffset = 5;  
                expandedMask = false(size(trueChunkY));
                
                for z = 1:size(validPairs, 1)
                    yVal = trueChunkY(validPairs(z,1));
                    lowVal  = yVal - rangeOffset;
                    highVal = yVal + rangeOffset;
                
                    expandedMask = expandedMask | (trueChunkY >= lowVal & trueChunkY <= highVal);

                end

                if ~isempty(expandedMask)
                    chunkX = trueChunkX(expandedMask);
                    chunkY = trueChunkY(expandedMask);

                    chunkCoordinates = [chunkX, chunkY];

                    chunks_toRemove = ismember(imageJOriginal, chunkCoordinates, 'rows');

                    if ~isempty(chunkCoordinates)
                        chunksTable = array2table(chunkCoordinates, 'VariableNames', {'X', 'Y'});
                        writetable(chunksTable, 'dark defects.xlsx', 'Sheet', 'Clusters', 'WriteMode', 'Append');
                    end
            
    
                    if ~isempty(chunks_toRemove)
                        pixelCoordinates = imageJOriginal(~chunks_toRemove, :);
                        numPixels = numel(pixelCoordinates(:,1));
                        pixelCoordinates(1,3) = numPixels;
                        pixelCoordinates(2:end, 3) = NaN;
                        pixelCoordinates(any(pixelCoordinates(:,3) == 0 & (1:size(pixelCoordinates,1))' ~= 1, 2), :) = [];

                        if ~isempty(pixelCoordinates)
                            pixelsTable = array2table(pixelCoordinates, 'VariableNames', {'X', 'Y', 'Count'});
                            writetable(pixelsTable, 'dark defects.xlsx', 'Sheet', 'Pixels', 'WriteMode', 'Overwrite');

                        end
                    end
                end
            end
        end
    end
end
