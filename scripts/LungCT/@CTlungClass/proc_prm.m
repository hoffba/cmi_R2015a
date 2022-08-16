function results = proc_prm(self)

[seg_ref,reg_seg] = self.getData('seg_ref','reg_seg');

% Check if file already exists:
fn_prm = self.getFileName('prm');
if exist(fn_prm,'file')
    writeLog(self.fn_log,'Loading PRM from file\n');
    self.dat.prm.mat = int8(niftiread(fn_prm));
else
    writeLot(self.fn_log,'Analyzing PRM\n')
    [ct_ref,reg_hom] = self.getData('ct_ref','reg_hom');
    mask = logical(seg_ref);
    if ~isempty(reg_seg)
        mask = mask & reg_seg;
    end
    [self.dat.prm.mat,~] = pipeline_PRM(ct_ref,self.dat.ct_ref.info,mask,reg_hom,...
        fullfile(self.procdir,[self.ID,'_PRM_Scatter']));
    % Save resulting PRM
    niftiwrite(self.dat.prm.mat,fn_prm(1:end-3),'Compressed',true);
end
if ~isempty(self.dat.prm.mat)
    % Tabulate 10-color PRM results
    writeLog(self.fn_log,'Tabulating 10-color PRM results...\n');
    results = lobeLoop(seg_ref,@(mask,prm,flag)tabulatePRM(mask,prm,flag),self.dat.prm.mat,1);

    % Map full PRM values (1:10) to (norm,fsad,emph,pd,ns)
    prmlabel = {'Norm', 'fSAD', 'Emph', 'PD',       'NS';...
                [1,2],  3,      [4,5],  [8,9,10],   6    };
    prm5 = self.dat.prm.mat;
    for i = 1:size(prmlabel,2)
        prm5(ismember(self.dat.prm.mat,prmlabel{2,i})) = i;
    end
    self.dat.prm.mat = prm5;
    clear prm5

    % QC PRM
    writeLog(self.fn_log,'Generating PRM Montage ...\n');
    QCmontage('prm',cat(4,self.dat.ct_ref.mat(:,:,ind),double(self.dat.prm.mat(:,:,ind))),...
        self.dat.ct_ref.info.PixelDimensions,...
        fullfile(self.procdir,[self.ID,'_PRM_Montage']));

    % Tabulate 5-color PRM results
    writeLog(self.fn_log,'Tabulating 5-color PRM results ...\n');
    T = lobeLoop(seg_ref,@(mask,prm,flag)tabulatePRM(mask,prm,flag),self.dat.prm.mat,0);
    results = addTableVarVal(results,T);
end

function T = tabulatePRM(mask,prm,flag)
if flag % 10-color
    vals = 1:10;
    tag = cellfun(@num2str,num2cell(vals),'UniformOutput',false);
else % 5-color
    vals = 1:5;
    tag = {'Norm', 'fSAD', 'Emph', 'PD', 'NS'};
end
nv = numel(vals);
vname = strcat('PRM_',tag);
T = table('Size',[1,nv],'VariableTypes',repmat({'double'},1,nv),...
    'VariableNames',vname);
np = nnz(mask);
for i = 1:numel(vals)
    T.(vname{i}) = nnz(prm(mask)==i)/np*100;
end


        