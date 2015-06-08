%
% function batch_MinkowskiFun(fnout,fov,img,r,thresh,ind,varargin)
function batch_MinkowskiFun(fnout,img,mask,fov,ind,r,thresh,varargin)
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
if nargin<8
    
    % This code starts the batch
    disp(['Starting batch Minkowski Functional analysis ... ',fname])
    batch(@batch_MinkowskiFun,0,{fnout,img,mask,fov,ind,r,thresh,1});
    disp(' ... done')
    
else
    
    % This code is run inside the batch
    disp(fname);
    disp(fpath);
    disp([]);
    
    % Initialize data:
    d = size(img);
    if isempty(ind)
        ind = 1:prod(d);
    end
    nth = length(thresh);
        
    % Loop over thresholds:
    MFmeans = zeros(nth,4);
    for ith = 1:nth
        
        disp(['Threshold (',num2str(ith),'/',num2str(nth),'): ',num2str(thresh(ith))]);
        
        BW = img > thresh(ith);
        [MF,labels] = minkowskiFun(BW,r,ind);
        MFmeans(ith,:) = mean(MF,1);
        
        % Save results as we go:
        ofname = [fname,'_th',num2str(thresh(ith))];

        % Loop over MF results to interpolate and save
%         [Xq,Yq,Zq] = meshgrid(1:d(2),1:d(1),1:d(3));
%         [Y,X,Z] = ind2sub(d,ind);
%         F = scatteredInterpolant(X,Y,Z,MF(:,1),'linear','none');
%         for imf = 1:size(MF,2)
% 
%             disp(['Processing: ',labels{imf}]);
% 
%             % Interpolate MF map to original dimensions
%             t = tic;
%             F.Values = MF(:,imf);
%             MFimg = F(Xq,Yq,Zq);
%             MFimg(~mask | isnan(MFimg)) = 0;
% 
%             disp(['     ',num2str(toc(t))]);
% 
%             % Save as MHD
%             oname = fullfile(fpath,[ofname,'_',labels{imf},'.mhd']);
%             saveMHD(oname,MFimg,'',fov);
%         end
        save(fullfile(fpath,[ofname,'.mat']),'MF','labels','ind','d');
    end
    save(fullfile(fpath,[fname,'_means.mat']),'MFmeans');
end