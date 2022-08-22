% Determine lobe tags for various lung segmentations
function lobe = getLobeTags(seg)
    utags = unique(seg(seg>0));
    
    if isempty(utags)
        lobeName = {}; lobeTag = {};
    elseif all(ismember(utags,[10,20]))
        lobeTag = {10,20};
        lobeName = {'RL', 'LL'};
    elseif all(ismember(utags,10:10:60)) % <- ******
        % YACTA
        if all(ismember(10:10:60,seg))
            lobeTag = {10,20,30,40,50,60,[10,20],[40,50]};
            lobeName = {'RUL','RML','RLL','LUL','LLi','LLL','RULplus','LULplus'};
        else
            lobeTag = {[10,20,30],[40,50,60]};
            lobeName = {'RL','LL'};
        end
    elseif all(ismember(utags,[1:5,30,255]))
        % YACTA
        lobeTag = num2cell(1:5);
        lobeName = {'LUL','LLL','RUL','RLL','RML'};
    elseif all(ismember(utags,[11,12,13,21,22,30]))
        % ImBio
        lobeTag = {11,12,13,21,22};
        lobeName = {'RUL','RLL','RML','LUL','LLL'};
    else
        error('Could not match valid segmentation labeling schema: %s',num2str(utags'));
    end
    lobe = struct('name',lobeName,'val',lobeTag);
end
