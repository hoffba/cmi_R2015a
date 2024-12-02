function D = report_prm(procdir,res,opts,R)

try
    % Check that Matlab's report generator toolbox is available
    if ~license('test','matlab_report_gen')
        fprintf('MATLAB Report Generator not available.\n');
        D = [];
        return
    end

    import mlreportgen.report.*;
    import mlreportgen.dom.*;

    if nargin==1
        [~,ID] = fileparts(procdir);
        opts = struct('ID',ID,...
                      'save_path','',...
                      'fn',struct('prmMontage',[ID,'_PRM_Montage.tif'],...
                                  'prmScatter',[ID,'_PRM_Scatter.tif']));
        res = readtable(fullfile(procdir,[ID,'_PipelineResults.csv']));
        res.Properties.RowNames = res.ROI;
    end

    fn_template = fullfile(fileparts(which('cmi')),'scripts','LungCT','ReportTemplates',...
        '+Pipeline_Report','@Chapter','resources','templates','pdf','default.pdftx');

    R_flag = nargin==4;
    if R_flag
        % Input as Chapter to Report object
        D = Chapter(TemplateSrc=fn_template);
    else
        % Generate Document object
        fname = fullfile(procdir,['PRM_Report_',opts.ID,'.pdf']);
        D = Document(fname,'pdf',fn_template);
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
    dp = DocumentPart('pdf',fn_template,'PipelineHeader');
    moveToNextHole(dp);
    append(dp,'Parametric Response Map (PRM)');
    moveToNextHole(dp);
    append(dp,string(datetime("today")));
    moveToNextHole(dp);
    pos = find(opts.ID=='_',1,"last");
    if isempty(pos)
        ID = opts.ID;
        tp = 'Unknown';
    else
        ID = opts.ID(1:pos-1);
        tp = opts.ID(pos+1:end);
        % Re-format date for display
        if numel(tp)==8 % 'YYYYmmdd'
            tp = char(datetime(tp,'InputFormat','yyyyMMdd','Format','d, MMMM yyyy'));
        end
    end
    append(dp,ID);
    moveToNextHole(dp);
    append(dp,tp);
    append(D,dp);
    
    % Main body of the report
    dp = DocumentPart('pdf',fn_template,'PRMbody');

    % PRM montage
    moveToNextHole(dp);
    Im = Image(opts.fn.prmMontage);
    Im.Style = {ScaleToFit(),Width('3.7in')};
    append(dp,Im);

    % PRM Scatterplot
    moveToNextHole(dp);
    Im = Image(opts.fn.prmScatter);
    Im.Style = {ScaleToFit(),Width('3.7in')};
    append(dp,Im);

    % Table of PRM results
    moveToNextHole(dp);
    resname = {'PRM_norm_pct', 'PRM_fsad_pct',   'PRM_emph_pct',    'PRM_pd_pct'};
    fmt =     {'%.1f',         '%.1f',           '%.1f',            '%.1f'};
    hdrname = {'NORMAL',       'FUNCTIONAL LDA', 'PERSISTENT LDA',  'INSPIRATION HDA'};
    bgcolor = {'#33cc33',      '#ffc40c',        '#ff3300',         '#cc00ff'};
    Nr = size(res,1);
    Nv = numel(resname);
    T = table('Size',[Nr,Nv],'VariableTypes',repmat({'cellstr'},1,Nv),'VariableNames',resname);
    for i = 1:Nv
        vname = resname{i};
        if ismember(vname,res.Properties.VariableNames)
            if isnumeric(res.(vname))
                nanmask = isnan(res.(vname));
                T.(vname) = arrayfun(@(x)sprintf(fmt{i},x),res.(vname),'UniformOutput',false);
                T.(vname)(nanmask) = {'-'};
            end
        end
    end
    T.Properties.RowNames = res.ROI';
    T = MATLABTable(T);
    T.Style = {FontSize('10pt')};
    % Add Units of % to first header spot
    T.Children(2).Children.Entries(1).Children.Children.Content = '%';
    T.Children(2).Children.Entries(1).InnerMargin = '10px';
    T.Children(2).Children.Entries(1).Children.Bold = true;
    % Edit table values
    for i = (1:Nv)+1
        % Set Header String
        T.Children(2).Children.Entries(i).Children.Content = hdrname{i-1};
        T.Children(2).Children.Entries(i).InnerMargin = '10px';
        T.Children(2).Children.Entries(i).Children.Color = bgcolor{i-1};
        T.Children(2).Children.Entries(i).Children.Bold = true;
        T.Children(2).Children.Entries(i).Children.Style = [T.Children(2).Children.Entries(i).Children.Style,{HAlign('center')}];
        % Remove extra quotes
        for j = 1:Nr
            T.Children(1).Children(j).Children(i).Children.Content = strip(T.Children(1).Children(j).Children(i).Children.Content,'''');
            % Center table entries
            T.Children(1).Children(j).Children(i).Style = [T.Children(1).Children(1).Children(i).Style,{HAlign('center')}];
        end
    end
    for i = 2:2:Nr
        T.Children(1).Children(i).Style = [T.Children(1).Children(i).Style,{BackgroundColor('#cfcfcf')}];
    end
    append(dp,T);
    
    append(D,dp);

    if R_flag
        append(R,D);
    else
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
catch err
    fprintf('PRM report generation failed.\n');
end