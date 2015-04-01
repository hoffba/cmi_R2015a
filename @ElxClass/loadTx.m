% ElxClass function
% Load Transformation from File
function stat = loadTx(self,fname)

stat = false;
if (nargin==1)
    [fname,fpath] = uigetfile('*.txt','Load TransformParameters.txt');
    if ischar(fname)
        fname = fullfile(fpath,fname);
    end
end
if ischar(fname) && exist(fname,'file')
    s = readElxTxt2Struct(fname);
    if ~isfield(s,'CenterOfRotationPoint')
        s.CenterOfRotationPoint = s.Origin;
    end
    if all(isfield(s,{'NumberOfParameters',...
                      'TransformParameters',...
                      'Size',...
                      'Spacing',...
                      'Origin',...
                      'CenterOfRotationPoint'}))
        self.Tx0 = s;
        stat = true;
    else
        error('Invalid parameter file.')
    end
end

