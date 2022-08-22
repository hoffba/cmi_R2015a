function setcase(self,id,procdir,fn)
% Set case information before running any processes

if ischar(id) && ischar(fn_base) && ischar(procdir)
    % Set case parameters
    self.id = id;
    self.fn_base = fn_base;
    
    % Re-initialize dat structure:
    self.initializeDat;
    
    self.procdir = procdir;
    if ~isfolder(procdir)
        mkdir(procdir);
    end
    
    self.elxdir = fullfile(procdir,"elxreg_"+fn_base+"_Ins");
    if ~isfolder(self.elxdir)
        mkdir(self.elxdir);
    end
    
    self.dat.seg_exp.yacta = fullfile(procdir,"yacta_"+fn_base+"_Exp");
    self.dat.seg_ins.yacta = fullfile(procdir,"yacta_"+fn_base+"_Ins");
    
    self.fn_log = fullfile(procdir,"log_"+fn_base);
    
    if (nargin==4)
        if ischar(fn)
            fn = {fn};
        end
        ct_str = {'ct_ref','ct_hom'};
        for i = 1:min(numel(fn),2)
            self.(ct_str).srcpath = fn{i};
        end
    end
end