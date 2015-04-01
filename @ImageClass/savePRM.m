% ImageClass function
function stat = savePRM(self,hObject,~)
% Save PRM statistics to the base workspace
stat = false;
if self.prm.check
    answer = '';
    if nargin==1
        answer = questdlg('Save Stats or PRM Image?','Save PRM','Stats','PRM Image','Stats');
    elseif ishandle(hObject)
        str = get(hObject,'Tag');
        if strcmp(str,'analysis_saveprm')
            answer = 'PRM Image';
        elseif strcmp(str,'analysis_prmstats')
            answer = 'Stats';
        end
    end
    if ~isempty(answer)
        if strcmp(answer,'Stats')
            [labels,pcts] = self.prm.getStats;
            assignin('base','prmStats',[labels;num2cell(pcts)]');
            stat = true;
        elseif strcmp(answer,'PRM Image')
            stat = cmi_save(false,self.prm.mat,{'PRM'},self.voxsz.*self.dims(1:3));
        end
    end
end
