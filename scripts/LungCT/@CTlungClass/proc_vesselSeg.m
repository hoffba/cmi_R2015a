function results = proc_vesselSeg(self)

fn_re_ins = self.getFileName('re_ct_hom');
fn_re_seg = self.getFileName('re_seg_hom');

[ct_hom,seg_hom] = self.getData('ct_hom','seg_hom');

if ~(isempty(ct_hom) || isempty(seg_hom))
    writeLog(self.fn_log,'Vessel analysis ...\n');
    if exist(fn_re_ins,'file') && exist(fn_re_seg,'file')
        tinfo = niftiinfo(fn_re_ins);
        ct_hom = niftiread(tinfo);
        seg_hom = niftiread(fn_re_seg);
    else
        tinfo = self.dat.ct_hom.info;
    end
    results = vesselSeg_BH( ct_hom , seg_hom , tinfo , procdir );
end
