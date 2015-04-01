% ElxClass function
% Write Transform .txt.file
function fname = writeTx0(self,fname)
% fname = writeTx0;
%       GUI for .txt file placement, returns full file path
% fname = writeTx0(fname);
%       Write to file fname

if ~isempty(self.Tx0)
    if nargin==1
        [fname,path] = uiputfile('*.txt','Save TransformParameter.txt',...
            'InitialTransform.txt');
        if fname
            fname = fullfile(path,fname);
        else
            fname = '';
        end
    end
    fname = writeElxStruct2Txt(self.Tx0,fname);
else
    fname = '';
end