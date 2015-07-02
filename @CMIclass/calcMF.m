% CMIclass function
function [MFout,labels] = calcMF(self,varargin)
% Calculate Minkowski Functionals

MFout = []; labels = {};
if (nargin==3) && ishandle(varargin{1})
    % If a GUI call, ask for user input:
    str = {'Image index (0=VOI)',                       num2str(self.vec);
           'Thresholds (vector of image values)',       '';...
           'Threshold operator (<,>,>=,<=,==)',         '>';...
           'Window Size ("inf" for overall analysis)',  'inf inf inf';...
           'Apply VOI to analysis? (y/n)',              'y';...
           'Non-VOI value (1/0)',                       '0'};
    answer = inputdlg(str(:,1),'Minkowski Functionals',1,str(:,2));
    if isempty(answer)
        return;
    else
        BWmode = true;
        if strcmp(answer{3},'<')
            BWmode = false;
        elseif ~strcmp(answer{3},'>')
            warning('Invalid input for threshold operator, using default (>).');
        end
        varargin = {'Vec',str2double(answer{1}),...
                    'Thresh',str2num(answer{2}),...
                    'BWmode',BWmode,...
                    'Window',str2num(answer{4}),...
                    'ApplyMask',strcmp(answer{5},'y'),...
                    'OutVal',logical(str2double(answer{6}))};
    end
end

p = inputParser;
addParameter(p,'Vec',self.vec,@isscalar);
addParameter(p,'Thresh',[],@isnumeric);
addParameter(p,'Window',[],@isvector);
addParameter(p,'ApplyMask',true,@islogical);
parse(p,varargin{:});
pp = p.Results;

[MFout,labels] = self.img.calcMF(self.vec,pp.Thresh,pp.Window,varargin{:});
assignin('base','MFvals',MFout);

