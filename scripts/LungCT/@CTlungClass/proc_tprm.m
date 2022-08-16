function results = proc_tprm(self)

results = [];

[prm,seg_ref] = self.getData('prm','seg_ref');

if ~(isempty(prm) || isempty(seg_ref))
    prmlabel = ["norm","fsad","emph","pd"];
    mflabel = ["v","s","b","x"];
    fn_tprm = fullfile(self.procdir,...
        string(self.ID) + ".tprm." + prmlabel + "." + mflabel' + string(self.fn_ext));
    if all(cellfun(@(x)exist(x,'file'),fn_tprm))
        writeLog(self.fn_log,'Loading tPRM from files ...\n');
        for iprm = 1:numel(prmlabel)
            for imf = 1:numel(mflabel)
                writeLog(self.fn_log,'   %s - %s\n',prmlabel(iprm),mflabel(imf))
                tprm = readNIFTI(fn_tprm(imf,iprm));
                str = prmlabel(iprm)+'_'+upper(mflabel(imf));
                T = lobeLoop(img.label,@(mask,tprm,str)tabulateTPRM(mask,tprm,str),tprm,str);
                results = addTableVarVal(results,T);
            end
        end
    else
        t = tic;
        writeLog(self.fn_log,'Generating tPRM ...\n');
    
        % Calculate MF values
        p = minkowskiFun(prm,...
            'thresh',   self.opts.tprm.thresh,...
            'tmode',    self.opts.tprm.tmode,...
            'n',        self.opts.tprm.n,...
            'gridsp',   self.opts.tprm.gridsp,...
            'voxsz',    self.dat.ct_ref.info.PixelSpacing,...
            'mask',     logical(seg_ref),...
            'prog',0);
    
        % Interpolate to maps
        for ithresh = 1:size(p.MF,1)
            for imf = 1:size(p.MF,2)
                writeLog(self.fn_log,'   %s - %s\n',prmlabel(ithresh),mflabel(imf));

                % Interpolate to image space
                writeLog(self.fn_log,'       Interpolating\n');
                tprm = grid2img(p.MF(ithresh,imf,:),p.ind,p.mask,3,1);

                % Save tPRM image
                writeLog(self.fn_log,'       Saving NIFTI\n');
                niftiwrite(tprm,char(fn_tprm(imf,ithresh)),self.dat.ct_ref.info)

                % Tabulate statistics
                writeLog(self.fn_log,'       Tabulating means\n');
                str = prmlabel(ithresh)+'_'+upper(mflabel(imf));
                T = lobeLoop(img.label,@(mask,tprm,str)tabulateTPRM(mask,tprm,str),tprm,str);
                results = addTableVarVal(results,T);
            end
        end
        writeLog(self.fn_log,'... tPRM complete (%s)\n',datestr(duration(0,0,toc(t)),'HH:MM:SS'));
    end
end 
    
    function T = tabulateTPRM(mask,tprm,str)
    vname = sprintf('tPRM_%s',str);
    T = table('Size',[1,1],'VariableTypes',{'double'},'VariableNames',{vname});
    T.(vname) = mean(tprm(mask));