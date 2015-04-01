function sortMRSdata(fpath)
% Sorts MR Solutions image data into scan folders

if nargin==0
    fpath = uigetdir(pwd);
end
if ~isempty(fpath) && ischar(fpath) && exist(fpath,'dir')
    fnames = dir(fpath); fnames(1:2) = []; % remove . and ..
    nf = length(fnames);
    if nf>0
        C = {'Scan ID','Acq Date','#Slc','FOVx','FOVy','Thk',...
             'Nx','Ny','FFTx','FFTy','PPL','PPLdir'};
            hw = waitbar(0,'Sorting MR Solutions data ...');
            for i = 1:nf
                sname = fullfile(fpath,fnames(i).name);
                dname = fullfile(fpath,datestr(fnames(i).datenum,'yyyymmdd'),...
                                    fnames(i).name(1:4));
                if ~exist(dname,'dir')
                    mkdir(dname);
                    % Read scan info from file for catalog
                    [~,~,ext] = fileparts(fnames(i).name);
                    switch ext
                        case '.SUR'
                            [~,~,par] = Get_mrd_3D1(sname);
                            slsh = strfind(par.PPL,'\');
                            [~,pplname] = fileparts(par.PPL(slsh(end)+1:end));
                            C(end+1,:) = {fnames(i).name(1:4),...
                                          par.date,...
                                          par.IM_TOTAL_SLICES,...
                                          par.IM_FOV(1),...
                                          par.IM_FOV(2),...
                                          par.IM_FOV(3),...
                                          par.IM_RESOLUTION(1),...
                                          par.IM_RESOLUTION(2),...
                                          par.IM_FFT(1),...
                                          par.IM_FFT(2),...
                                          pplname,...
                                          par.PPL(1:slsh(end))};
                        case '.MRD'
                        case '.dcm'
                    end
                end
                movefile(sname,fullfile(dname,fnames(i).name));
                waitbar(i/nf,hw);
            end
            delete(hw);
        if size(C,1)>1
            cmi_csvwrite(fullfile(fpath,'MRSdata_catalog.csv'),C);
        end
    end
end