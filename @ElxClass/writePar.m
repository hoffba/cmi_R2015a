% ElxClass function
% Write Parameter.txt files
function fnout = writePar(self,fname)


fnout = {};
if ~isempty(self.Schedule)
    if nargin==1
        [fname,path] = uiputfile('*.txt','Save TransformParameter.txt',...
            'InitialTransform.txt');
        if fname
            fname = fullfile(path,fname);
        else
            fname = '';
        end
    end
    ns = length(self.Schedule);
    fnout = cell(1,ns);
    if ns>1
        nstr = cellstr(num2str((1:ns)'));
    else
        nstr = {''};
    end
    [path,fname,~] = fileparts(fname);
    for i = 1:ns
        fnout{i} = writeElxStruct2Txt(self.Schedule{i},...
                        fullfile(path,[fname,nstr{i},'.txt']));
    end
end
