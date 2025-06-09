function res = dipg_final(T)
% Input: 
%   T = table with columns:
%           CaseID  : <SubjectID>_<DateStamp>
%           ADC     : File name for ADC map
%           Anat    : File name for Anatomical image (for VOI and registration)
%           VOI     : File name for VOI (provided by Benison Lau)
%   * All files should be in the the path: <basedir>\<SubjectID>\<DateStamp>

    % Handle inputs
    if ~nargin
        [fn,path] = uigetfile('*.csv','Select Processing Catalog File');
        if fn
            T = fullfile(path,fn);
        else
            return;
        end
    end
    if ischar(T) && isfile(T)
        T = readtable(T);
    end
    if ~(istable(T) && all(ismember({'CaseID','ADC','Anat','VOI'},T.Properties.VariableNames)))
        error('Invalid input.');
    end

    % Initialize data path
    basedir = 'R:\CGalban_Lab\Cancer\DMG_DIPG\20230930_BLau\Processed\batch20241223';
    
    % Initialize CMI object for processing
    cmiobj = CMIclass;
    cmiobj.img.prm.setOpts('thresh',[2,3,1,-.55; 2,3,1,.55],...
                           'prmmap',{[false false], 'ADC_-';...
                                     [true  false], 'ADC_0';...
                                     [true  true], 'ADC_+'},...
                           'cutoff',[2,0.0001,3 ; 3,0.0001,3],...
                           'cmap',flip(eye(3)),...
                           'statchk',false);

    % Loop over subjects
    BLi = [];
    res = [];
    for i = 1:size(T,1)
        % Update baseline case if needed
        if isempty(BLi) || ~strcmp(extractBefore(T.CaseID{i},'_'),extractBefore(T.CaseID{BLi},'_'))
            BLi = i;
        end

        t = dipg_case_final(basedir,T(i,:),T(BLi,:),cmiobj); 
        if isempty(t)
            
        else
            res = addResultToTable(res,t);
        end
    end

    % Delete CMIclass object
    cmiobj.delete;

    % Write results to study folder
    if ~isempty(res)
        writetable(res,fullfile(basedir,'DIPG.AllResults.xlsx'));
    end
