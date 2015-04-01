% CMIclass function
function [MFout,labels] = calcMF(self,varargin)
% Calculate Minkowski Functionals

MFout = []; labels = {};
if (nargin==3) && ishandle(varargin{1})
    % If a GUI call, ask for user input:
    str = {'Image index (0=VOI)',                       num2str(self.vec);
           'Thresholds (vector of image values)',       '';...
           'Window Size ("inf" for overall analysis)',  'inf inf inf';...
           'Apply VOI to analysis? (y/n)',              'y'};
    answer = inputdlg(str(:,1),'Minkowski Functionals',1,str(:,2));
    if isempty(answer)
        return;
    else
        varargin = {'Vec',str2double(answer{1}),...
                    'Thresh',str2num(answer{2}),...
                    'Window',str2num(answer{3}),...
                    'ApplyMask',strcmp(answer{4},'y')};
    end
end

p = inputParser;
addParamValue(p,'Vec',self.vec,@isscalar);
addParamValue(p,'Thresh',[],@isnumeric);
addParamValue(p,'Window',[],@isvector);
addParamValue(p,'ApplyMask',true,@islogical);
parse(p,varargin{:});
pp = p.Results;

[MFout,labels] = self.img.calcMF(self.vec,pp.Thresh,pp.Window,varargin{:});
if ndims(MFout)==5
    self.update4Dslider;
else
    a = [labels ; num2cell(MFout)];
    msgbox(sprintf('%s  =  % .4f\n',a{:}),'Minkowski Functionals','none');
end
