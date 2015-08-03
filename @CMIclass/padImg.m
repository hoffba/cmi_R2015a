% CMIclass function
function padImg(self,~,~)
% Pad the current image

if self.img.check
    answer = inputdlg({'Dimension (1,2,3):','Length (voxels)','Pad Value (for each d4)'},...
        'Pad Image:',1,{'3','2','replicate'});
    if ~isempty(answer)
        
        dim = cell2mat(textscan(answer{1},'%u'));
        r = cell2mat(textscan(answer{2},'%u'));
        pval = strsplit(answer{3});
        tval = cellfun(@str2double,pval,'UniformOutput',false);
        nani = ~cellfun(@isnan,tval);
        pval(nani) = tval(nani);
        if length(r)==1
            r = r*ones(1,length(dim));
        end
        
        self.img.pad(dim,r,pval);
        
        rnew = zeros(1,3);
        rnew(dim) = r;
        self.slc(self.orient) = self.slc(self.orient) + rnew(self.orient);
        self.setView(self.orient);
    end
end