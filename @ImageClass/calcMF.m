% ImageClass function
function [MF,p] = calcMF(self,vec,varargin)
% function [MF,labels] = calcMF(self,vec,varargin)
% Calculate Minkowski Functionals
% Inputs:   vec = 4D image index to perform analysis on
%           ithresh = vector of thresholds to use
%           n = window radius (voxels) (inf = global analysis)
%           optional Name/Value pairs:
%               'ApplyMask' = check for use of existing VOI
%               'tmode'     = string describing thresholding operator
%                               - starts with 'p' --> percentiles
%                               - {'<','<=','>','>=','==','~='}
%               'defVal'    = logical value of voxels outside the VOI

MF = [];
p = [];

if nargin<4 % GUI for user input:
    str = {'MF Data ("VOI", "PRM", or image index)',        'prm'   ;...
           'Values (for binarization)',                     ''      ;...
           'Mode (p if percentiles, then <,>,>=,<=,==)',    '=='    ;...
           'Window radius ("inf" for global analysis)',     'inf'   ;...
           'Grid spacing',                                  ''      ;...
           'Apply VOI to analysis? (y/n)',                  'y'     ;...
           'Non-VOI value (1/0)',                           '0'     };
    answer = inputdlg(str(:,1),'Minkowski Functionals',1,str(:,2));
    if isempty(answer)
        return;
    else
        vec = str2double(answer{1});
        if isnan(vec)
            vec = answer{1};
        end
        r = sscanf(answer{4},'%f');
        varargin = {'thresh',   str2num(answer{2}),...
                    'tmode',    answer{3},...
                    'n',        r,...
                    'gridsp',   sscanf(answer{5},'%f'),...
                    'defVal',   logical(str2double(answer{7}))};
        voichk = strcmpi(answer{6},'y');
    end
end

if self.check
    if ischar(vec)
        switch lower(vec)
            case 'voi'
                timg = self.mask.mat;
            case 'prm'
                timg = self.prm.mat;
            otherwise
                error('Invalid data selection.');
        end
    elseif ismember(vec,1:self.dims(4))
        timg = self.mat(:,:,:,vec);
    else
        error('Invalid image vector selected.');
    end
    if voichk
        mask = self.mask.mat;
    end
    gchk = any(isinf(r));
    varargin = [varargin,{'voxsz',self.voxsz,'mask',mask,'prog',true}];
    
    if gchk % Global analysis:

        % Calculate MF:
        [MF,p] = minkowskiFun(timg,varargin{:});
            
        if nargout==0
            % Save results to base workspace:
            p.MF = MF;
            assignin('base','MFresults',p);
        end
        
    else
        
        % Run in batch - will take too long to wait for
        
        % User select where to save results:
        fnout = fullfile(self.dir,[self.name,'_MF',datestr(now(),'yyyymmddHHMMSS'),'.mat']);
        [fnout,fpath] = uiputfile('*.mat','Save MF results as:',fnout);
        
        if ~isempty(fnout)
            batch_MinkowskiFun(fullfile(fpath,fnout),timg,varargin{:});
        end
    end
end

