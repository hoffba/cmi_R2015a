function [prmind,C] = prm_seg(vals,vec,thresh,cutoff)
% Generates PRM indices for multi-dimensional PRM
% Inputs:
%       vals:       [basevals,yvals,...]; (nD columns)
%       vec:        vector of image #'s corresponding to columns of vals
%       thresh:     [n x 4] matrix : row = [d1 d2 m b]
%                               evaluating : d2 > d1*m + b
%       cutoff:     optional - [n x 3] matrix : row = [dim min max]
% Outputs:
%       prmind:     vector of PRM region indices
%       C:          logical matrix of unique (by row) region threshold stats

prmind = []; C = [];
if (nargin>2)
    [np,ndim] = size(vals);
    nth = size(thresh,1);
    if (ndim==length(vec)) && (nth>0)
        
        % Segment by thresholds
        ind = false(np,nth);
        for i = 1:nth
            if thresh(i,1)==0
                d1 = 1;
            else
                d1 = find(vec==thresh(i,1),1);
            end
            if thresh(i,2)==0
                d2 = 1;
            else
                d2 = find(vec==thresh(i,2),1);
            end
            ind(:,i) = ( vals(:,d2) > vals(:,d1)*thresh(i,3) + thresh(i,4) );
        end
        
        % Determine cropped values
        if (nargin<4) || ~isnumeric(cutoff) || (size(cutoff,2)~=3)
            cutoff = [];
        end
        cropind = false(np,1);
        % Replace 0 with first vec (dynamic image)
        [~,i,~] = unique(cutoff(:,1));
        cutoff = cutoff(i,:);
        for i = 1:size(cutoff,1)
            d1 = find(vec==cutoff(i,1),1);
            if ~isempty(d1)
                cropind = cropind | (vals(:,d1)<cutoff(i,2)) ...
                                  | (vals(:,d1)>cutoff(i,3));
            end
        end
        
        % Determine regions:
        [C,~,prmind] = unique(ind,'rows');
        prmind(cropind) = 0;
        
    end
end
    
 