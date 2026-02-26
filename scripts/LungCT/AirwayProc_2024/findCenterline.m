function CL = findCenterline(A)



% Select largest connected component BEFORE skeletonization

% Find all connected components in 3D
CC = bwconncomp(A, 26);  % 26-connectivity for 3D

if CC.NumObjects == 0
    error('No airway components found in the airway mask.');
elseif CC.NumObjects > 1
    % Find the largest component
    numPixels = cellfun(@numel, CC.PixelIdxList);
    [maxPixels, largestIdx] = max(numPixels);
    
    % Log information about discarded components
    totalPixels = sum(numPixels);
    discardedPixels = totalPixels - maxPixels;
    fprintf('Connected component analysis:\n');
    fprintf('  Total components found: %d\n', CC.NumObjects);
    fprintf('  Largest component: %d voxels (%.1f%% of total)\n', maxPixels, 100*maxPixels/totalPixels);
    fprintf('  Discarded: %d voxels from %d smaller components\n', discardedPixels, CC.NumObjects-1);
    
    % Create new mask with only the largest component
    A_clean = false(size(A));
    A_clean(CC.PixelIdxList{largestIdx}) = true;
    A = A_clean;
    clear A_clean;
    
    % Save QC info about component selection
    comp_info.num_components = CC.NumObjects;
    comp_info.largest_size = maxPixels;
    comp_info.total_size = totalPixels;
    comp_info.percent_kept = 100*maxPixels/totalPixels;
    comp_info.discarded_sizes = numPixels(numPixels ~= maxPixels);
    save(fullfile(outD,[caseID,'_component_QC.mat']),'comp_info');
else
    fprintf('Airway mask has single connected component - no filtering needed.\n');
end

CL = bwskel(A, 'MinBranchLength', 10);  % Remove short spurious branches during skeletonization
