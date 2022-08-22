function stat = proc_scatnet(self)

stat = ~isempty(self.dat.scatnet.mat);
if ~stat
    fn_save = self.getFileName('scatnet');
    if exist(fn_save,'file')
        self.dat.scatnet.mat = niftiread(fn_save);
        stat = true;
    else
        [ct,seg] = self.getData('ct_ref','seg_ref');
        if ~isempty(ct) && ~isempty(seg)
            self.dat.scatnet.mat = ScatNet(ct,seg,0);
            
            % Save resulting AT map
            niftiwrite(self.dat.scatnet.mat,fn_save(1:end-3),self.dat.ct_ref.info,'Compressed',true);
            
            stat = true;
        end
    end
end