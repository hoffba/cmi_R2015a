function seg = CTlung_Segmentation(method,ct,info,id,savepath,logfn)

seg = false;

if nargin<6
    logfn = '';
end

% Switch string input method to corresponding numeric index:
if ischar(method)
    segmethod = {'segLungHuman','getRespiratoryOrgans','DL_Craig','YACTA','TotalSegmentator','PTK'};
    method = find(strcmp(method,segmethod));
end

if method==4 && (nargin<5 || ~ischar(savepath) || ~isfolder(savepath))
    warning('YACTA must have a valid input path for saving results. Trying DL instead.');
    method = 3;
end

% Check for GPU (required if method==3)
if method==3 && ~gpuDeviceCount
    warning('DL_lung_segmentation was selected, but no GPU is available. Reverting to getRespiratoryOrgans method.');
    method = 2;
end

switch method
    case 1 % segLungHuman
        writeLog(logfn,'- Generating VOI from Step02_segLungHuman_cjg ...\n');
        seg = segLungHuman_cjg_bh(1,ct);
    case 2 % getRespiratoryOrgans
        writeLog(logfn,'- Generating VOI from getRespiratoryOrgans ...\n');
        seg = getRespiratoryOrgans(ct);
    case 3 % DL_Craig
        fprintf(logfn,'- Generating VOI from DL_lung_segmentation ...\n');
        seg = DL_lung_segmentation(ct);
    case 4 % YACTA
        writeLog(logfn,'- Generating VOI from YACTA ...\n');
        ydir = fullfile(savepath,['yacta_',id]);
        if ~isfolder(ydir)
            mkdir(ydir);
        end
        tname = fullfile(ydir,sprintf('%s.mhd',id));

        % Check for gapped data. Remove gaps for YACTA processing
        gapchk = false;
        orient = info.orient;
        ind = find(any(ct>-1000,[1,2]));
        if numel(ind)~=size(ct,3)
            gapchk = true;
            orient = orient*diag([1,1,(ind(2)-ind(1)),1]);
        end

        saveMHD(tname,ct(:,:,ind),id,info.fov,orient);
        yacta(tname,'wait');
%         yacta(tname,'wait','airways','renderer','hide','exportlabels','yactascp');

        resname = dir(fullfile(ydir,'*lung_lobes*explabels.mhd'));
        if ~isempty(resname)
            seg = cmi_load(1,[],fullfile(ydir,resname(end).name));
        end
        if nnz(seg)<10^3 || ~all(ismember([10,20,30,40,60],seg))
            writeLog(logfn,'- - Lobe segmentation failed. Checking for lungs...');
            resname = dir(fullfile(ydir,'*lung_right_left*explabels.mhd'));
            if ~isempty(resname)
                seg = cmi_load(1,[],fullfile(ydir,resname(end).name));
            end
            if nnz(seg)>10^3 && all(ismember([10,20],seg))
                writeLog(logfn,' Using R/L segmentation\n');
            else
                writeLog(logfn,' FAILED\n');
                seg = CTlung_Segmentation(2,ct,info,id,savepath,logfn);
            end
        end

        % Add gaps back into gapped segmentations
        if gapchk
            tseg = seg;
            seg = zeros(size(ct));
            seg(:,:,ind) = tseg;
        end

    case 5 % TotalSegmentator
        writeLog(logfn,'- Generating VOI from TotalSegmentator ...\n');
        seg = TotalSegmentator(ct,info,id,savepath);
        if ischar(seg) % Failed to run
            writeLog(logfn,'  TotalSegmentator: %s\n',seg);
            seg = CTlung_Segmentation(2,ct,info,id,savepath,logfn);
        end

    case 6 % PTK
        writeLog(logfn,'- Generating VOI from PTK ...\n');
        [label,res] = pipeline_PTKlobes(fn);

    otherwise
end

    