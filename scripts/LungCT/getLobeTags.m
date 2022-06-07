% Determine lobe tags for various lung segmentations
function lobe = getLobeTags(seg)
    utags = unique(seg(seg>0));
    
    if isempty(utags)
        lobeName = {}; lobeTag = {};
    elseif all(ismember(utags,[10,20]))
        lobeTag = [10 20];
        lobeName = {'RL', 'LL'};
    elseif all(ismember(utags,10:10:60))
        % YACTA?
        lobeTag = 10:10:60;
        lobeName = {'RUL','RML','RLL','LUL','LLi','LLL'};
    elseif all(ismember(utags,1:5))
        % YACTA
        lobeTag = 1:5;
        lobeName = {'LUL','LLL','RUL','RLL','RML'};
    elseif all(ismember(utags,[11,12,13,21,22,30]))
        % ImBio
        lobeTag = [11,12,13,21,22];
        lobeName = {'RUL','RLL','RML','LUL','LLL'};
    else
        error('Could not match valid segmentation labeling schema: %s',num2str(utags'));
    end
    lobe = struct('name',lobeName,'val',num2cell(lobeTag));
end
