function D = report_prm(procdir,res,opts,R)
    import mlreportgen.report.*;
    import mlreportgen.dom.*;

    fn_template = fullfile(opts.report_path);

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
        tp = '';
    else
        ID = opts.ID(1:pos-1);
        tp = opts.ID(pos+1:end);
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
        % Remove extra quotes
        for j = 1:Nr
            T.Children(1).Children(j).Children(i).Children.Content = strip(T.Children(1).Children(j).Children(i).Children.Content,'''');
        end
    end
    for i = 2:2:Nr
        T.Children(1).Children(i).Style = {BackgroundColor('#cfcfcf')};
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
