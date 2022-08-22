function results = run(self,modules)
% Process pipeline modules for requested measures

results = [];
dat_flds = fieldnames(self.dat);
if nargin==1
    % User GUI to select measurements
    [ind,ok] = listdlg('ListString',dat_flds);
    if ok
        modules = dat_flds(ind);
    else, return;
    end
else
    modules = dat_flds(ismember(dat_flds,modules));
end

% DICOM conversion
ind = ismember(modules,{'ct_ref','ct_hom'});
if any(ind)
    self.proc_read_CT(modules(ind));
end

% Lung segmentation
ind = ismember(modules,{'seg_ref','seg_hom'});
if any(ind)
    self.proc_seg(modules(ind));
end

% Loop over other modules
for i = 1:numel(modules)
    switch modules{i}
        case 'scatnet'
            self.proc_scatnet;
        case 'unreg'
            self.proc_unreg;
        case {'vessels','csa'}
            self.proc_vesselSeg;
        case 'airways'
            self.proc_airwaySeg;
        case {'jac','reg_ins'}
            self.proc_registration;
        case 'prm'
            self.proc_prm;
        case 'tprm'
            self.proc_tprm;
        otherwise
    end
end