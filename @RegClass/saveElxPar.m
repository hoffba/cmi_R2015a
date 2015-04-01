% RegClass function
function saveElxPar(self,ind,fname)
% saveElxPar
% saveElxPar(ind,'path/fname.txt')

ns = length(self.elxObj.Schedule);
if ns>0
    if ishandle(ind)
        ind = self.ind;
        if ns>1
            answer = questdlg('Save which parameters?','Save Elastix Parameters',...
                              'Selected','All','Cancel','Selected');
            if strcmp(answer,'All')
                ind = 1:length(self.elxObj.Schedule);
            elseif strcmp(answer,'Cancel')
                ind = [];
            end
        end
        if ~isempty(ind)
            [fname,path] = uiputfile('*.txt','Save Elastix Parameters',...
                                     fullfile(self.odir,'ElastixParameters.txt'));
            if fname==0
                ind = [];
            else
                fname = fullfile(path,fname);
            end
        end
    elseif ~isnumeric(ind) || ~(ischar(fname) || iscellstr(fname))
        ind = [];
    end
    if ~isempty(ind)
        self.elxObj.savePar(ind,fname);
    end
end