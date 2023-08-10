% Determine lobe tags for various lung segmentations
function lobe = getLobeTags(seg)
    utags = unique(seg(seg>0));
    
    if isempty(utags)
        lobeName = {}; lobeTag = {};
    elseif islogical(seg) || all(utags == 1)
        lobeTag = {1};
        lobeName = {'WholeLung'};
    elseif all(ismember(utags,[10,20]))
        lobeTag = {10,20};
        lobeName = {'RL', 'LL'};
    elseif all(ismember(utags,10:10:60)) % <- ******
        % YACTA
        lobeTag = {10,20,30,40,50,60,[10,20],[40,50]};
        lobeName = {'RUL','RML','RLL','LUL','LLi','LLL','RULplus','LULplus'};
    elseif all(ismember(utags,[1:5,30,255]))
        % YACTA
        lobeTag = num2cell(1:5);
        lobeName = {'LUL','LLL','RUL','RLL','RML'};
    elseif all(ismember(utags,1:6))
        % YACTA -- OLD
        lobeTag = num2cell(1:6);
        lobeName = {'LUL','LLi','LLL','RUL','RML','RLL'};
    elseif all(ismember(utags,[11,12,13,21,22,30,255]))
        % ImBio
        lobeTag = {11,12,13,21,22};
        lobeName = {'RUL','RLL','RML','LUL','LLL'};
    else
        error('Could not match valid segmentation labeling schema: %s',num2str(utags'));
    end
    lobe = struct('name',lobeName,'val',lobeTag);
end
