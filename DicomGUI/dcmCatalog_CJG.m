function C = dcmCatalog_CJG(filt)
% Catalogs desired DICOM header values

img=ImageClass;

hdrstr = {'PatientID','StudyID','StudyDate','SeriesNumber','AcquisitionNumber',...
    'StudyDescription','SeriesDescription',...
    'ScanOptions','FilterType','ConvolutionKernel','KVP',...
    'SliceThickness','Rows','Columns','PixelSpacing'};
nfld = length(hdrstr);
C = [{'Directory'},hdrstr(1:(end-1)),{'dy','dx','Slices'}];

% Find folders containing DICOMs
opath = uigetdir(pwd,'Select folder containing all DICOMS:');
try
    if opath~=0
        hw = waitbar(0,'Finding DICOM folders ...');
        D = dirtree(opath,filt); % D: cell array of strings (directories)
        ndir = length(D);
        
        % Loop over all DICOM directories
        for idir = 1:ndir
            waitbar(idir/ndir,hw,['Processing folder ',num2str(idir),...
                ' of ',num2str(ndir)]);
            
            fnames = dir(fullfile(D{idir},filt));
            if strcmp(fnames(1).name,'.')
                fnames(1:2) = [];
            end
            nf = length(fnames);
            if nf>1
                fields = {}; clear('dinfo');
                for ifn = nf:-1:1
                    
                    if ~isdir(fullfile(D{idir},fnames(ifn).name))&&...
                            isdicom(fullfile(D{idir},fnames(ifn).name))
                        
                        tinfo = dicominfo(fullfile(D{idir},fnames(ifn).name));
                        tfields = fieldnames(tinfo);
                        newfields = setdiff(tfields,fields); enf = ~isempty(newfields);
                        misfields = setdiff(fields,tfields); emf = ~isempty(misfields);
                        if (ifn==nf)
                            fields = tfields;
                        elseif (enf || emf)
                            % Need to adjust info structure fields to match
                            if emf
                                for k = 1:length(misfields)
                                    tinfo.(misfields{k}) = [];
                                end
                            end
                            if enf
                                for k = 1:length(newfields)
                                    dinfo(ifn).(newfields{k}) = [];
                                end
                                fields = [fields;newfields]; %#ok<*AGROW>
                                dinfo = orderfields(dinfo,tinfo);
                            else
                                tinfo = orderfields(tinfo,dinfo);
                            end
                        end
                        dinfo(ifn) = tinfo;
                        
                        % Find unique acquisition numbers
                        if isfield(dinfo,'AcquisitionNumber')
                            acqn = [dinfo(:).AcquisitionNumber];
                            acqu = unique(acqn);
                            if isempty(acqn)
                                acqn = 1; acqu = 1;
                            end
                        else
                            acqn = 1; acqu = 1;
                        end
                        % Store relevant parameters for each acquisition
                        for j = 1:length(acqu)
                            ii = find(acqn==acqu(j));
                            tC = cell(1,nfld);
                            
                            for ifld = 1:nfld
                                if isfield(dinfo,hdrstr{ifld})
                                    tC{ifld} = dinfo(ii(1)).(hdrstr{ifld}); % originally dinfo CJG
                                end
                            end
                            if isempty(tC{end})
                                dd = nan(1,2);
                            else
                                dd = tC{end};
                                
                                C = [C;D(idir),tC(1:(end-1)),{dd(1),dd(2),length(ii)}];
                                assignin('base','Catalog',C);
                                % CJG modifications
                                TempC=[D(idir),tC(1:(end-1)),{dd(1),dd(2),length(ii)}];
                                
                                % This is for data from other sites using
                                % the last folder name as the file name
                                % (e.g. is CF Stanford Data)
                                a=D{idir};
                                b=strfind(a,'/');
                                c=a((1+b(end)):end);
                                c(isspace(c))=[];
                                str_file=[c,'.mhd'];TempC{2}=c;
                                dir_Temp=fullfile(opath,'processed',TempC{2});
                                if ~isdir(dir_Temp)
                                    mkdir(dir_Temp);
                                end
%                                 img.loadImg(0,fullfile(D{idir},fnames(1).name))
%                              
                                'load files'
                                D{idir}
                                'what will be saved'
                                fullfile(opath,'processed',TempC{2},str_file)
%                                 img.saveImg(1,fullfile(opath,'processed',TempC{2},str_file),1)
                                
                                % below works for UM pulled data
                                %                 if (strcmp(TempC{9},'HELICAL MODE')|strcmp(TempC{9},'AXIAL MODE'))&j==1
                                %                     TempC{9}
                                %                     TempC{8}
%                                                     dir_Temp=fullfile(opath,'processed',TempC{2});
%                                                     if ~isdir(dir_Temp)
%                                                         mkdir(dir_Temp);
%                                                     end
                                %                     if (strcmp(TempC{8},'INSPIRATION')||strcmp(TempC{8},'HELICAL INSPIRATION'))
                                %                         str_file=[TempC{2},'_',TempC{3},'_',num2str(TempC{5}),'_',TempC{11},'_Acqu',num2str(length(acqu)),'_Ins.mhd']
%                                                         img.loadImg(0,fullfile(D{idir},fnames(1).name))
%                                                         img.saveImg(1,fullfile(opath,'processed',TempC{2},str_file),1)
                                %                     elseif (strcmp(TempC{8},'EXPIRATION')||strcmp(TempC{8},'HELICAL EXPIRATION'))
                                %                         str_file=[TempC{2},'_',TempC{3},'_',num2str(TempC{5}),'_',TempC{11},'_Acqu',num2str(length(acqu)),'_Exp.mhd']
                                %                         img.loadImg(0,fullfile(D{idir},fnames(1).name))
                                %                         img.saveImg(1,fullfile(opath,'processed',TempC{2},str_file),1)
                                %                     end
                                %                 end
                            end
                        end
                    end
                end
                clear info
            end
        end
        
        delete(hw);pause(10);
        [fname,path] = uiputfile('*.csv','Save Catalog Data:');
        if fname~=0
            cmi_csvwrite(fullfile(path,fname),C);
        end
    end
catch err
    'failed'
    
end
