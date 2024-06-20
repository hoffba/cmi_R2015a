function pipeline_report(procdir,res,opts,str)
% Function to generate analysis report for Mi-TAP Pipeline
% Inputs:
%   dp_body_name =  <char>  name of dptemplate body section in template file
%   procdir =       <char>  case processing directory
%   res =           <table> analysis results
%   opts =          <struct> Mi-TAP options
%   R =             (optional) mlreportgen.report.Report to append Document to

import mlreportgen.report.*;
import mlreportgen.dom.*;

chapters = {
    'ins',      'Inspiratory CT',                   'INSbody';...
    'exp',      'Expiratory CT',                    'EXPbody';...
    'airways',  'Airways (YACTA)',                  'AIRbody';...
    'vessels',  'Blood Vessels',                    'VESbody';...
    'prm',      'Parametric Response Map (PRM)',    'PRMbody';...
    'tprm',     'Topological PRM (tPRM)',           'tPRMbody'
    };

ind = 1:size(chapters,1);
if nargin==4
    ind = find(strcmpi(chapters(:,1),str),1);
end

if ~isempty(ind)
    % Determine PDF file name
    pdfname = fullfile(procdir,[opts.ID,'_Report.pdf']);

    % Initialize report
    R = [];
    if numel(ind)>1
        R = Report(pdfname,'pdf');
        tR = PageMargins();
        tR.Top = ".5in";
        tR.Bottom = ".5in";
        tR.Left = ".5in";
        tR.Right = ".5in";
        tR.Header = "0pt";
        R.Layout.PageMargins = tR;
        
        tR = TitlePage;
        tR.Title = 'Pipeline Results';
        tR.Author = opts.ID;
        append(R,tR);
        append(R,TableOfContents);
    end

    % Loop over analyses
    for i = ind
        D = genReport(chapters{i,:},procdir,res,opts);
    end

end


function D = genReport(tag,title_str,docprt,procdir,res,opts)

switch lower(str)
    case 'ins'
    case 'exp'
    case 'airways'
    case 'vessels'
    case 'prm'
        ch_name = 'Parametric Response Map (PRM)';
        dp_body_name = 'PRMbody';
    case 'tprm'
end

R_flag = nargin==5;
if R_flag
    % Input as Chapter to Report object
    D = Chapter(TemplateSrc=opts.report_path,Title=ch_name);
else
    % Generate Document object
    fname = fullfile(procdir,['PRM_Report_',opts.ID,'.pdf']);
    D = Document(fname,'pdf',opts.report_path);
    open(D);

    % Set up page margins
    D.CurrentPageLayout.PageMargins.Left = '0.5in';
    D.CurrentPageLayout.PageMargins.Right = '0.5in';
    D.CurrentPageLayout.PageMargins.Top = '0.5in';
    D.CurrentPageLayout.PageMargins.Bottom = '0.5in';
    D.CurrentPageLayout.PageMargins.Header = '0in';
    D.CurrentPageLayout.PageMargins.Footer = '0in';
end

% Set up header
dp = DocumentPart('pdf',opts.report_path,'PipelineHeader');
moveToNextHole(dp);
append(dp,ch_name);
moveToNextHole(dp);
append(dp,string(datetime("today")));
moveToNextHole(dp);
pos = find(opts.ID=='_',1,"last");
if isempty(pos)
    ID = opts.ID;
    tp = '';
else
    ID = opts.ID(1:pos-1);
    tp = opts.ID(pos+1:end);
end
append(dp,ID);
moveToNextHole(dp);
append(dp,tp);
append(D,dp);

% Set up body of report/chapter
dp = DocumentPart('pdf',opts.report_path,dp_body_name);
while ~strcmp(moveToNextHole(dp),'#end#')
    sethole(dp,res,opts);
end

% Add docpart to document
append(D,dp);

if R_flag
    % Append document to report
    append(R,D);
else % Save document to file
    close(D);

    % Save a copy to collation directory
    if ~isempty(opts.save_path)
        svdir = fullfile(opts.save_path,'Pipeline_Report');
        if ~isfolder(svdir)
            mkdir(svdir);
        end
        copyfile(fname,svdir);
    end
end


function sethole(dp,res,opts)
import mlreportgen.report.*;
import mlreportgen.dom.*;

switch dp.CurrentHoleId

% Exp
    case 'Exp_Montage'
    case 'Exp_Table'

% Ins
    case 'Ins_Montage'
    case 'Ins_Table'

% Airways
    case 'Air_Fig'
    case 'Air_Table'

% Vessels
    case 'Vessel_Fig'
    case 'Vessel_Table'

% PRM
    case 'PRM_Montage'
        Obj = Image(opts.fn.prmMontage);
        Obj.Style = {ScaleToFit(),Width('3.7in')};
    case 'PRM_Scatter'
        Obj = Image(opts.fn.prmScatter);
        Obj.Style = {ScaleToFit(),Width('3.7in')};
    case 'PRM_Table'
        resname = {'PRM_norm_pct', 'PRM_fsad_pct',   'PRM_emph_pct',    'PRM_pd_pct'};
        fmt =     {'%.1f',         '%.1f',           '%.1f',            '%.1f'};
        hdrname = {'NORMAL',       'FUNCTIONAL LDA', 'PERSISTENT LDA',  'INSPIRATION HDA'};
        bgcolor = {'#33cc33',      '#ffc40c',        '#ff3300',         '#cc00ff'};
        Nr = size(res,1);
        Nv = numel(resname);
        Obj = table('Size',[Nr,Nv],'VariableTypes',repmat({'cellstr'},1,Nv),'VariableNames',resname);
        for i = 1:Nv
            vname = resname{i};
            if ismember(vname,res.Properties.VariableNames)
                if isnumeric(res.(vname))
                    nanmask = isnan(res.(vname));
                    Obj.(vname) = arrayfun(@(x)sprintf(fmt{i},x),res.(vname),'UniformOutput',false);
                    Obj.(vname)(nanmask) = {'-'};
                end
            end
        end
        Obj.Properties.RowNames = res.ROI';
        Obj = MATLABTable(Obj);
        Obj.Style = {FontSize('10pt')};
        % Add Units of % to first header spot
        Obj.Children(2).Children.Entries(1).Children.Children.Content = '%';
        Obj.Children(2).Children.Entries(1).InnerMargin = '10px';
        Obj.Children(2).Children.Entries(1).Children.Bold = true;
        % Edit table values
        for i = (1:Nv)+1
            % Set Header String
            Obj.Children(2).Children.Entries(i).Children.Content = hdrname{i-1};
            Obj.Children(2).Children.Entries(i).InnerMargin = '10px';
            Obj.Children(2).Children.Entries(i).Children.Color = bgcolor{i-1};
            Obj.Children(2).Children.Entries(i).Children.Bold = true;
            % Remove extra quotes
            for j = 1:Nr
                Obj.Children(1).Children(j).Children(i).Children.Content = strip(Obj.Children(1).Children(j).Children(i).Children.Content,'''');
            end
        end
        for i = 2:2:Nr
            Obj.Children(1).Children(i).Style = {BackgroundColor('#cfcfcf')};
        end

% tPRM

    otherwise
        Obj = [];
end

% Append object to hole
if ~isempty(Obj)
    append(dp,Obj);
end




return


% Expiration
if isfile(opts.fn.expMontage)
    Ch = Chapter('Title','Expiration');
    Im = Image(opts.fn.expMontage);
    Im.Style = {Width('4.8in'),HAlign('center')};
    vars = {'Exp_Vol',          '%.1f',     0, 'Volume (L)';...
            'Exp_HU',           '%.0f',     0, 'Mean HU';...
            'Exp_856',          '%.1f',     0, 'HU<-856 (%)';...
            'scatnetAT_pct',    '%.1f',     0, 'SN_AT (%)'};
    T = pipeline_table(res,vars);
    add(Ch,Im);
    add(Ch,newline);
    add(Ch,T);
    append(R,Ch);
end

% Inspiration
if isfile(opts.fn.insMontage)
    Ch = Chapter('Title','Inspiration');
    Im = Image(opts.fn.insMontage);
    Im.Style = {Width('4.8in'),HAlign('center')};
    vars = {'Ins_Vol',        '%.1f',   0, 'Volume (L)';...
            'Ins_HU',         '%.0f',   0, 'Mean HU';...
            'Ins_950',        '%.1f',   0, 'Ins_950 (%)';...
            'Ins_810',        '%.1f',   0, 'Ins_810 (%)';...
            'Ins_810low',     '%.1f',   0, 'Ins_810low (%)';...
            'Ins_500',        '%.1f',   0, 'Ins_500 (%)';...
            'Ins_GGOI',       '%.1f',   0, 'GGOI (%)';...
            'Ins_FIBI',       '%.1f',   0, 'FIBI (%)';...
            'scatnetEmph_pct','%.1f',   0, 'SN_Emph (%)'};
    T = pipeline_table(res,vars);
    add(Ch,Im);
    add(Ch,newline);
    add(Ch,T);
    append(R,Ch);
end

% Airways table
vars = {'Wall_pct_1', '%.1f',   0, 'Wall_pct';...
        'Pi10',       '%.1f',   0, 'Pi10';...
        'Pi15',       '%.1f',   0, 'Pi15';...
        'BEI',        '%.1f',   0, 'BEI';...
        'BEI_gen',    '%.1f',   0, 'BEI_gen'};
if any(ismember(vars(:,1),res.Properties.VariableNames))
    Ch = Chapter(Title='YACTA Airways Results',Numbered=false);
    T = pipeline_table(res,vars);
    add(Ch,T);
    append(R,Ch);
end

% Vessel Analysis
vars = {'VOLUME_L',                 '%.1f',   0, 'Lung Volume (L)';...
        'VESSEL_VOLUME_L',          '%.1f',   3, 'Vessel Vol (mL)';...
        'PER_EMPH',                 '%.1f',   0, 'Emph (%)';...
        'NUM_VESSELS',              '%.0f',   0, '# Vessels';...
        'NUM_COMPONENTS',           '%.0f',   0, '# Components';...
        'NUM_ENDPOINTS',            '%.0f',   0, '# Endpoints';...
        'CSA_EXP_A',                '%.2f',   0, 'CSA Exp A';...
        'CSA_EXP_B',                '%.2f',   0, 'CSA Exp B';...
        'VESSEL_VOLUME_5DOWN_L',    '%.1f',   3, 'Vessel Vol 5down (mL)';...
        'VESSEL_VOLUME_5UP_L',      '%.1f',   3, 'Vessel Vol 5up (mL)'};
if any(ismember(vars(:,1),res.Properties.VariableNames))
    Ch = Chapter(Title='Blood Vessel Analysis',Numbered=false);
    T = pipeline_table(res,vars);
    add(Ch,T);
    append(R,Ch);
end

% Image Registration
if isfile(opts.fn.regMontage)
    Ch = Chapter(Title='Image Registration',Numbered=false);
    Im = Image(opts.fn.regMontage);
    Im.Style = {Width('4.8in'),HAlign('center')};
    add(Ch,Im);
    vars = {'Jac_mean',           '%.3f', 0, 'Mean Jacobian';...
            'scatnetEmphReg_pct', '%.1f', 0, 'SN-Emph-Reg';...
            'dBlood_mean',        '%.1f', 0, 'Mean dBlood (HU)'};
    if any(ismember(vars(:,1),res.Properties.VariableNames))
        T = pipeline_table(res,vars);
        add(Ch,newline);
        add(Ch,T);
    end
    append(R,Ch);
end

% PRM
if isfile(opts.fn.prmMontage)
    Ch = Chapter(Title='Parametric Response Map (PRM)',Numbered=false);
    T = {Image(opts.fn.prmMontage) , ' ' , Image(opts.fn.prmScatter)};
    T{1}.Style = {Width('3.7in'),HAlign('center')};
    T{3}.Style = {Width('3.7in')};
    T = Table(T);
    T.HAlign = 'center';
    T.entry(1,1).Style = {Width('3.7in'),Height('3.7in'),HAlign('center')};
    T.entry(1,1).Style = {Width('0.1in')};
    T.entry(1,1).Style = {Width('3.7in'),Height('3.7in')};
    add(Ch,T);
    vars = {'PRM_norm_pct', '%.1f', 0, 'Norm (%)';...
            'PRM_fsad_pct', '%.1f', 0, 'FSAD (%)';...
            'PRM_emph_pct', '%.1f', 0, 'Emph (%)';...
            'PRM_pd_pct',   '%.1f', 0, 'PD (%)';...
            'PRM_ns_pct',   '%.1f', 0, 'NS (%)'};
    if any(ismember(vars(:,1),res.Properties.VariableNames))
        T = pipeline_table(res,vars);
        add(Ch,newline);
        add(Ch,T);
    end
    append(R,Ch);
end

%tPRM
if any(startsWith(res.Properties.VariableNames,'tPRM_'))
    Ch = Chapter(Title='Topological PRM (tPRM)',Numbered=false);
    vars = {'tPRM_norm_V_mean', '%.1f',  2, 'Volume (V, %)';...
            'tPRM_norm_S_mean', '%.2f',  1, 'Surface (S, x10)';...
            'tPRM_norm_B_mean', '%.2f',  2, 'Breadth (B, x100)';...
            'tPRM_norm_X_mean', '%.2f',  3, 'Euler (X, x10^4)'};
    prmname = {'Norm','fSAD','Emph','PD'};
    T = cell(5,3);
    for iprm = 1:4
        irow = 3*mod(iprm-1,2) + 1;
        icol = 2*ceil(2*iprm/4) - 1;

        table_title = Text(['PRM-',prmname{iprm}]);
        table_title.Style = [table_title.Style,{Bold(),Underline('single'),FontSize('10pt')}];

        vars(:,1) = strcat('tPRM_',lower(prmname{iprm}),'_',{'V','S','B','X'}','_mean');
        tprm_table = pipeline_table(res,vars);

        T{irow,icol} = table_title;
        T{irow+1,icol} = tprm_table;
    end
    T = Table(T);
    T.Style = {HAlign('center')};
    T.Children(3).Height = '12pt';
    T.Children(1).En
    bgcolors = {'#98FB98','#EEFC5E','#EF9A9A','#BF92E4'};
    for iprm = 1:4
        irow = 3*mod(iprm-1,2) + 1;
        icol = 2*ceil(2*iprm/4) - 1;
        % irow = mod(uint64(i)*2-1,T.NRows);
        % icol = ceil(2*i/4);
        T.Children(irow).Entries(icol).Style = [T.Children(irow).Entries(icol).Style,...
                                                {BackgroundColor(bgcolors{iprm})}];
    end
    add(Ch,T);
    append(R,Ch);
end

close(R);



function T = pipeline_table(res,vars)
% Create Table with desired formatting
    import mlreportgen.report.*;
    import mlreportgen.dom.*;

    % Find desired variables in results table
    Nr = size(res,1);
    Nv = size(vars,1);
    T = table('Size',[Nr,Nv],'VariableTypes',repmat({'cellstr'},1,Nv),'VariableNames',vars(:,1)');
    for i = 1:Nv
        vname = vars{i,1};
        if ismember(vname,res.Properties.VariableNames)
            if isnumeric(res.(vname))
                nanmask = isnan(res.(vname));
                T.(vname) = arrayfun(@(x)sprintf(vars{i,2},x),res.(vname)*10^vars{i,3},'UniformOutput',false);
                T.(vname)(nanmask) = {'-'};
            end
        end
    end

    % Set ROI to row name and remove from table variables
    T.Properties.RowNames = res.ROI';

    % Create mlreportgen Table()
    T = Table(T);
    T.Style = {FontSize('10pt'),HAlign('center')};
    T.TableEntriesInnerMargin = '300pt';
    % Set Header Style
    T.Children(1).Style = {Bold(),TextOrientation('down'),VAlign('bottom')};
    for i = 1:Nv
        % Set Header String
        T.Children(1).Entries(i+1).Children.Content = [vars{i,4},'  |'];
        % Set Column Style
        T.ColSpecGroups(i+1).Style = {Width('40pt')};
        % Remove extra quotes
        for j = (1:Nr) + 1
            T.Children(j).Entries(i+1).Children.Content = strip(T.Children(j).Entries(i+1).Children.Content,'''');
        end
    end

function prm_page(fn,opts)
    import mlreportgen.report.*;
    import mlreportgen.dom.*;

    D = Document(fn,'pdf','PRM_Report.pdftx');
    open(D);

    D.CurrentPageLayout.PageMargins.Left = '0.5in';
    D.CurrentPageLayout.PageMargins.Right = '0.5in';
    D.CurrentPageLayout.PageMargins.Top = '0.5in';
    D.CurrentPageLayout.PageMargins.Bottom = '0.5in';

    dp = DocumentPart(D,'PageHeader');
    moveToNextHole(dp);
    append(dp,datetime("today"));
    moveToNextHole(dp);
    pos = find(opts.ID=='_',1,"last");
    append(dp,opts.ID(1:pos-1));
    moveToNextHole(dp);
    append(dp,opts.ID(pos+1:end))
    append(D,dp);
    
    dp = DocumentPart(D,'MainSection');

    moveToNextHole(dp);
    Im = Image(opts.fn.prmMontage);
    Im.Style = {Im.Style,Width('3.7in')};
    append(dp,Im);

    moveToNextHole(dp);
    Im = Image(opts.fn.prmScatter);
    Im.Style = {Im.Style,Width('3.7in')};
    append(dp,Im);

    moveToNextHole(dp);
    vars = {'PRM_norm_pct', '%.1f', 0, 'Norm (%)';...
            'PRM_fsad_pct', '%.1f', 0, 'FSAD (%)';...
            'PRM_emph_pct', '%.1f', 0, 'Emph (%)';...
            'PRM_pd_pct',   '%.1f', 0, 'PD (%)';...
            'PRM_ns_pct',   '%.1f', 0, 'NS (%)'};
    append(dp,pipeline_table(res,vars))

