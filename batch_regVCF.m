function C = batch_regVCF(regObj)

% Structure containing data to coregister
tdata = struct('PatientName',...
                   {'13453','13894','13895','13896','13967'},...
               'Vertebra',{{'L2','T8','T9'},...
                           {'L3','T5','L4','T9'},...
                           {'L2','T12','L4','T10'},...
                           {'T7','T9'},...
                           {'L1','L3','T11'}},...
               'Stat',{[1,0,0],...
                       [1,1,1,1],...
                       [1,1,1,1],...
                       [1,1,],...
                       [1,1,0]});
tdir = '/Volumes/projects/Skeletal_Studies/clin_bone/Breast_Mets-Cathy_von_Poznak/FLDs/';
if ~exist(tdir,'dir')
    tdir = uigetdir(pwd,'Select folder containing patient folders.');
end

if tdir~=0
    % Initialize coregistration objects/settings:
    C = {};
    if (nargin==0) || ~isa(regObj,'RegClass')
        regObj = RegClass;
    end
    regObj.loadElxPar(fullfile(tdir,'VCF-ElastixParameters.txt'));
    set(regObj.h.checkbox_wait,'Value',0);
    regObj.cmiObj(1).img.prm.setOpts('thresh',[1,0,1,-100;...
                                               1,0,1, 100]);
    % Loop over patients and their vertebrae:
    for i = 1:length(tdata)
        for vi = 1:length(tdata(i).Vertebra)
            if tdata(i).Stat(vi)
                bdir = fullfile(tdir,tdata(i).PatientName,tdata(i).Vertebra{vi});
                
                % Find all time point images:
                fnames = dir(fullfile(bdir,'*.mhd')); fnames = {fnames(:).name};
                    % Ignore VOIs and already-coregistered images
                    ind = cellfun(@(x)isempty(strfind(x,'_RR'))&&isempty(strfind(x,'_VOI')),fnames);
                    fnames(~ind) = [];
                    nf = length(fnames);
                    
                % Load baseline image/VOI
                regObj.cmiObj(1).loadImg(false,fullfile(bdir,fnames{1}));
                regObj.cmiObj(1).loadMask(fullfile(bdir,[fnames{1}(1:end-4),'_VOI.mhd']));
                
                % Adjust schedule for low-res data:
                tpar = regObj.elxObj.Schedule{1}.FixedImagePyramidSchedule;
                R = regObj.cmiObj(1).img.voxsz(3)/regObj.cmiObj(1).img.voxsz(1);
                fZ = ceil(R./[1,2,4]);
                tpar([3,6,9]) = fZ;
                regObj.elxObj.setPar(1,'FixedImagePyramidSchedule',tpar);
                
                % Loop over coregistrations:
                for fi = 2:nf
                    
                    % Make sure elxreg_ folder exists
                    [~,bname] = fileparts(fnames{fi});
                    elxdir = fullfile(bdir,['elxreg_',bname]);
                    if ~exist(elxdir,'dir')
                        mkdir(elxdir);
                    end
                    
                    % Load moving image and VOI:
                    regObj.cmiObj(2).loadImg(false,fullfile(bdir,fnames{fi}));
                    regObj.cmiObj(2).loadMask(fullfile(bdir,[bname,'_VOI.mhd']));
                    
                    % Match VOIs for initial transform:
                    regObj.VOI2Tform;
                    regObj.setT('11',1);
                    regObj.setT('22',1);
                    regObj.setT('33',1);
                    
                    % Adjust schedule depending on resolutions:
                    tpar = regObj.elxObj.Schedule{1}.MovingImagePyramidSchedule;
                    R = regObj.cmiObj(2).img.voxsz(3)/regObj.cmiObj(2).img.voxsz(1);
                    fZ = ceil(R./[1,2,4]);
                    tpar([3,6,9]) = fZ;
                    regObj.elxObj.setPar(1,'MovingImagePyramidSchedule',tpar);
                    
                    % Start coregistration:
                    regObj.runElx;
                    
                    % Append result to Reference data set
                    regObj.cmiObj(1).loadImg(true,{fullfile(elxdir,[bname,'_RR.mhd'])});
                    
                end
                % Generate PRM statistics
                disp('check')
                
            end
        end
    end
end
