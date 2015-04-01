function stats = jobProcMouseLung(img,fov,labels,fname,svdir)
% Processes mouse lung CT images from the Bruker SkyScan 1176
% Performs the following steps:
%   - NOVA filter both images
%   - Save filtered images as DICOM
%   - Grow VOI to measure mean/mode HU and lung volume
%       * values are output in stats structure
%           use fetchOutputs function to retrieve stats

stats = [];
if ~exist(svdir,'dir')
    disp(svdir)
    disp('Save directory does not exist! ');
elseif ~all(ismember(labels,{'CT-Ins','CT-Exp'}))
    disp('Could not find correct labels!');
else
    
    stats = struct([]);
    
    % Determine file names:
    usi = strfind(fname,'_');
    fname = fname(1:(usi(end)-1));
    disp(['Dataset = ',fname]);
    iins = find(strcmp('CT-Ins',labels),1);
    iexp = find(strcmp('CT-Exp',labels),1);
    disp(['Ins = ',num2str(iins)])
    disp(['Exp = ',num2str(iexp)])
    disp(svdir)
    if ~isempty(iins)
        stats(iins).name = fullfile(svdir,[fname,'_Ins']);
    end
    if ~isempty(iexp)
        stats(iexp).name = fullfile(svdir,[fname,'_Exp']);
    end
    
    % Perform a NOVA filter
    n = 3;
    run = 2;
    sDev0 = 200;
    p = 0;
    d = 3;
    img = bFiltNOVA(img,n,run,sDev0,p,d);
    
    % Inspiration image:
    info = struct('PatientID',fname,...
                  'Modality','CT',...
                  'StudyDescription','Lung Disease',...
                  'StudyInstanceUID',dicomuid,...
                  'StudyID','Biphasic CT');
              
    for i = 1:size(img,4)
        stats(i).saved = saveMHD(stats(i).name,img(:,:,:,i),labels(i),fov,info);
    end
    
%     % Expiration image:
%     stats(iexp).saved = saveMHD(stats(iexp).name,...
%                                  img(:,:,:,iexp),labels(iexp),fov,info,0);
end
