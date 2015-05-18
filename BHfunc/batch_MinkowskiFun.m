%
function batch_MinkowskiFun(fnout,fov,img,r,thresh,ind,varargin)
% function batch_MinkowskiFun(fname,img,n,thresh,mask)
%   Runs a batch job to calculate local Minkowski Functionals over image
%   using defined image thresholds.
% Inputs:
%   fname
%   img
%   r
%   thresh
%   ind

[fpath,fname] = fileparts(fnout);
if nargin<7
    
    % This code starts the batch
    disp(['Starting batch Minkowski Functional analysis ... ',fname])
    batch(@batch_MinkowskiFun,0,{fnout,fov,img,r,thresh,ind,1});
    disp(' ... done')
    
else
    
    % This code is run inside the batch
    disp(fname);
    disp(fpath);
    disp([]);
    
    % Initialize data:
%     d = size(img);
    if isempty(ind)
%         ind = 1:prod(d);
        ind = 1:numel(img);
    end
    nth = length(thresh);
        
    % Loop over thresholds:
    for ith = 1:nth
        
        disp(['Threshold (',num2str(ith),'/',num2str(nth),'): ',num2str(thresh(ith))]);
        
        BW = img > thresh(ith);
        [MF,labels] = minkowskiFun(BW,r,ind);
%         nmf = length(labels);
        
        % Save results as we go:
        save(fullfile(fpath,[fname,'_th',num2str(thresh(ith)),'.mat']),'MF','labels');
%         MFout = zeros(d);
%         for j = 1:nmf
%             MFout(ind) = MF(:,j);
%             saveMHD(fullfile(fpath,[fname,'_th',num2str(thresh(ith)),'_',...
%                                     labels{j},'.mhd']),MFout,'',fov);
%         end
    end
    
end