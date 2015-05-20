function status = cmi_recon(fnames)
% function to reconstruct raw MRI fid files
% valid imaging seqences: fsems, sems_iadc, gems, sems

status = false;
if nargin==0
    reconpath = uigetdir(pwd,'Select directory for reconstruction:');
    if strcmp(reconpath(end-3:end),'.fid')
        fnames = {reconpath};
    else
        fnames = dir(fullfile(reconpath,'*.fid'));
        fnames = {fnames(:).name}';
        [sel,ok] = listdlg('ListString',fnames);
        if ok
            fnames = strcat(reconpath,filesep,fnames(sel));
        else
            fnames = {};
        end
    end
end
if ~isempty(fnames)
    nf = size(fnames,1);
    hw = waitbar(0,['Completed 0 of ' num2str(nf) ' reconstructions']);
    for i = 1:nf % loop over all .fid folders
        if exist(fnames{i},'dir')
            
            % Read parameters from procpar
            pp = readPROCPAR(fnames{i});
            
            if ~isempty(pp) && ~any(strcmp(pp.pslabel,{'scout',''}))
            
                % Read complex k-space data:
                kdata = readFID(fnames{i},pp);
                [d(1),d(2),d(3),necho,narr] = size(kdata);
            
                % Determine navigator data:
                navchk = false;
                if strncmp(pp.seqfil,'sems_iadc',9)
                    navchk = true;
                    kdata = flip(kdata,4);
                elseif isfield(pp,'navigator') && strcmp(pp.navigator,'y')
                    navchk = true;
                else
                    narr = narr*necho;
                    kdata = reshape(kdata,[d(1),d(2),d(3),narr]);
                end
                
                % Fourier reconstruction:
                d(1:2) = max(d(1:2));
                img = zeros([d,narr]);
                for iarr = 1:narr
                    for islc = 1:d(3)
                        if navchk
%                             img(:,:,islc,iarr) = navcorrect(kdata(:,:,islc,1,iarr),...
%                                                             kdata(:,:,islc,2,iarr));
                            img(:,:,islc,iarr) = abs(fftshift(fft2(kdata(:,:,islc,1,iarr),...
                                                                    d(1),d(2))));
                        else
                            img(:,:,islc,iarr) = abs(fftshift(fft2(kdata(:,:,islc,iarr),...
                                                                    d(1),d(2))));
                        end
                    end
                end
                img = flip(flip(img,1),2);
                
                % Determine array labels
                if isfield(pp,'gdiff') && (length(pp.gdiff)==narr) ...
                        && (length(pp.gdiff)>1)
                    if strncmp(pp.seqfil,'sems_iadc',9)
                        b = nan(1,narr);
                        b(ismember(pp.gdiff,[4, 4.5, 6, 7, 8, 8.1])) = 120;
                        b(ismember(pp.gdiff,[11, 20, 20.1, 21, 24])) = 1200;
                    else
                        b = num2cell(pp.bvalue);
                    end
                    label = cellstr(strcat('b=',num2str(b(:))));
                else
                    arrchk = ~(isempty(pp.array) || isempty(pp.array{1}));
                    if arrchk
                        for iarr = 1:length(pp.array)
                            if isfield(pp,pp.array{iarr})
                                arrchk = (length(pp.(pp.array{iarr}))==narr) && arrchk;
                            end
                        end
                    end
                    if arrchk
                        label = cell(narr,1);
                        for iarr = 1:length(pp.array)
                            if iarr>1
                                label = strcat(label,',');
                            end
                            label = strcat(label,pp.array{iarr},'=',...
                                cellstr(num2str(pp.(pp.array{iarr})')) );
                        end
                    else
                        label = repmat(pp.seqfil,1,narr);
                    end
                end
                
                % Determine FOV
                dz = (max(pp.pss)-min(pp.pss));
                fov = [pp.lro , pp.lpe , dz+dz/(pp.ns-1)]*10;
                
                % Save images as FLD
                [~,oname] = fileparts(fnames{i});
                saveFLD(fullfile(fnames{i},[oname,'.fld']),img,label(:)',fov);
                
                waitbar(i/nf,hw,['Completed ',num2str(i),' of ' num2str(nf) ' reconstructions']);
            end
        end
    end
    close(hw);
end


