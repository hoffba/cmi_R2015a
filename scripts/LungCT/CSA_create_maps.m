function [image] = CSA_create_maps(V)
%createCSA_maps.m Creates CSA maps from vessel segmentations
% Inputs:
%   V : (3D binary image) vessel segmentation
%
% Outputs: 
%   image : CSA skeleton 
    
    fprintf('   Finding connected segments\n');
    skel = bwskel(V > 0); %Create skeleton of vessels
    mask = bwconncomp(skel, 26); %Break into connected pieces
    
    pieces = mask.PixelIdxList; 
    np = numel(pieces);
    fprintf('   Found %u pieces.',np);

    D = bwdist(V == 0); %distance from each vessel pixel to the background ie radius
    csa = (D.*0.625).^2.*pi;
    
    %Skeleton to fill with CSA 
    image = zeros(size(V));
    
    stxt = sprintf('%%%uu (%%4.1f%%%% - %%s)',numel(num2str(np)));
    nbsp = 0;
    fprintf('   Processing ');
    t = tic;
    for i = 1:np %For each connected piece
        
        if numel(pieces{i})>1
        
            % Isolate just this piece on 3D matrix
            tmp = zeros(size(V));
            tmp( pieces{i}) = true;
            try
                [~,~,links] = Skel2Graph3D(tmp,0);
                for lidx = 1:length(links)
                    link = links(lidx).point;
                    image(link) = mean(csa(link));
                end
            catch E
                disp('check');
            end
            
        end
        
        % Print progress:
        if ~mod(i,round(np/1000))
            txt = sprintf(stxt,i,i/np*100,duration(0,0,toc(t)));
            fprintf([repmat('\b',1,nbsp),'%s'],txt);
            nbsp = numel(txt);
        end
    end
    fprintf('\n');
    image = single(image);
end