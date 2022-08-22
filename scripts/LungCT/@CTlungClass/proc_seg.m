function stat = proc_seg(self,mod_str)

stat = false(1,numel(mod_str));
if nargin<2 || ~(iscellstr(mod_str) || ischar(mod_str) || isstring(mod_str)) || ~all(ismember(mod_str,{'seg_exp','seg_ins'}))
    mod_str = mod_str(ismember(mod_str,{'seg_exp','seg_ins'}));
end

% Loop over exp/ins
for i = 1:numel(mod_str)
    
    % Grab file name to look for / save as
    fn_save = self.getFileName(mod_str);
    
    ct_str = ['ct_',mod_str{i}(end-2:end)]; % hom/ref ct image corresponding to segmentation
    if exist(fn_save,'file')
        writeLog(self.fn_log,'Reading segmentation from file: %s',fn_save);
        info = niftiinfo(fn_save);
        self.dat.(mod_str{i}).mat = niftiread(fn_save);
        if isempty(self.dat.(sprintf('ct_%s',mod_str{i}(end-2:end))).info)
            self.dat.(ct_str).info = info;
        end
        stat(i) = true;
    else % if file not found, run segmentation
        ct = self.getData(ct_str);
        stat(i) = ~isempty(ct);
        if stat(i)
            writeLog(self.fn_log,'Segmentation: %s',mod_str{i});
            info = struct('fov',self.dat.(ct_str).info.PixelDimensions.*self.dat.(st_str).info.ImageSize,...
                          'orient',(self.dat.(ct_str).info.Transform.T * diag([-1 -1 1 1]))');
            id = [self.fn_base,'.',self.dat.(mod_str{i}).tag];
            ydir = self.dat.(sprintf('seg_%s',mod_str{i})).yacta;
            self.dat.(mod_str{i}).mat = CTlung_Segmentation(self.opts.seg_method,ct,info,id,ydir,self.fn_log);
            
            % Write resulting nifti file
            if ~isempty(self.dat.(mod_str{i}).mat)
                niftiwrite(self.dat.(mod_str{i}).mat,fn_save(1:end-3),self.dat.(ct_str).info,'Compressed',true);
            end
        end
    end
end

stat = all(stat);



