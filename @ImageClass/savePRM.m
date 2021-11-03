% ImageClass function
function stat = savePRM(self,hObject,~)
% Save PRM statistics to the base workspace
stat = false;
if self.prm.check
    answer = '';
    if nargin==1
        answer = questdlg('Save Stats or PRM Image?','Save PRM','Stats','Image','VOI','Stats');
    elseif ishandle(hObject)
        str = get(hObject,'Tag');
        if strcmp(str,'analysis_saveprm')
            answer = questdlg('Save Stats or PRM Image?','Save PRM','Image','VOI','Image');
        elseif strcmp(str,'analysis_prmstats')
            answer = 'Stats';
        end
    elseif ischar(hObject) && ismember(hObject,{'Stats','Image','VOI'})
        answer = hObject;
    end
    if ~isempty(answer)
        if strcmp(answer,'Stats')
            [labels,pcts] = self.prm.getStats;
            assignin('base','prmStats',[labels;num2cell(pcts)]');
            stat = true;
        elseif strcmp(answer,'Image')
            stat = cmi_save(false,self.prm.mat,{'PRM'},self.voxsz.*self.dims(1:3),self.orient);
        elseif strcmp(answer,'VOI')
            ofname = fullfile(self.dir,[self.prm.dlabels{2},'_PRM']);
            labl = self.prm.prmmap(:,2);
            for i = 1:self.prm.nprm
                stat = cmi_save(true,self.prm.mat==i,labl(i),self.voxsz.*self.dims(1:3),self.orient,...
                    [ofname,labl{i},'_VOI.hdr']);
            end
        end
    end
end
