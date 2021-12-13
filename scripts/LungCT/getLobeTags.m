% Determine lobe tags for various lung segmentations
function lobe = getLobeTags(seg)
    utags = unique(seg(seg>0));
    if all(ismember(1:5,utags))
        % YACTA
        lobeTag = 1:5;
        lobeName = {'LUL','LLL','RUL','RLL','RML'};
    elseif all(ismember([11,12,13,21,22],utags))
        % ImBio
        lobeTag = [11,12,13,21,22];
        lobeName = {'RUL','RLL','RML','LUL','LLL'};
    elseif all(ismember(10:10:50,utags))
        % YACTA?
        lobeTag = 10:10:50;
        lobeName = {'LUL','LLL','RUL','RLL','RML'};
    else
        error('Could not match valid segmentation labeling schema: %s',num2str(utags'));
    end
    lobe = struct('name',lobeName,'val',num2cell(lobeTag));
end
