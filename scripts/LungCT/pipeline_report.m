function R = pipeline_report(R,ch_str,varargin)

import mlreportgen.report.*;
import mlreportgen.dom.*;

nv = numel(varargin);
N = nv/2;
if nargin==1 && (ischar(R) || isstring(R))
    % Initialize the Report with subject ID
    [R_path,ID] = fileparts(R);
    R = Report(fullfile(R_path,['PipelineOutput_',ID]),'pdf');
    tR = TitlePage;
    tR.Title = 'Pipeline Results';
    tR.Author = ID;
    append(R,tR);
    append(R,TableOfContents);
elseif isa(R,'mlreportgen.report.Report') && isvalid(R)
    chap = [];
    switch ch_str
% Segmentations
        case 'seg'
            if N>0 && N==round(N)
                chap = Chapter(Title='Segmentations',Numbered=false);
                chap.Layout.Landscape = true;
                lot = repmat({' '},1,2*N-1);
                for i = 1:N
                    fim = FormalImage( Image = varargin{2*i},...
                                       Caption = Text(varargin{2*i-1}) );
                    lot{2*i-1} = fim;
                end
                lot = Table(lot);
                append(chap,lot);
            end
% Airways table
        case 'airway'
            if nv==1 && istable(varargin{1})
                T = varargin{1};
                chap = Chapter(Title='YACTA Airways Results',Numbered=false);
                chap.Layout.Landscape = true;
                colstyle = {{'BEI','BEI_gen'},'%.0f'};
                append(chap,pipeline_table(T,colstyle));
            end
% ScatterNet for Air Trapping
        case 'scatnetAT'
            if nv==1 && istable(varargin{1})
                T = varargin{1};
                chap = Chapter(Title='Scatternet for Air Trapping on Expiration CT',Numbered=false);
                chap.Layout.Landscape = true;
                append(chap,pipeline_table(T));
            end
% ScatterNet for Emph on INS
        case 'scatnetEmph'
            if nv==1 && istable(varargin{1})
                T = varargin{1};
                chap = Chapter(Title='Scatternet for Emph on Inspiratory CT',Numbered=false);
                chap.Layout.Landscape = true;
                append(chap,pipeline_table(T));
            end
% Vessel Analysis
        case 'vessel'
            if nv==1
                T = varargin{1};
                chap = Chapter(Title='Blood Vessel Analysis',Numbered=false);
                chap.Layout.Landscape = true;
                colstyle = {T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'NUM_')),'%.0f'};
                append(chap,pipeline_table(T,colstyle));
            end
        case 'unreg'
            if nv>0 && N==round(N)
                chap = Chapter(Title='Unregistered Statistics',Numbered=false);
                chap.Layout.Landscape = true;
                lot = repmat({' '},2,2*N-1);
                for i = 1:N
                    txt = Text(varargin{2*i-1});
                    txt.Bold = true;
                    txt.Underline = 'single';
                    txt.FontSize = '12pt';
                    lot{1,2*i-1} = txt;
                    lot{2,2*i-1} = pipeline_table(varargin{2*i});
                end
                lot = Table(lot);
                lot.Style = {Width("100%")};
                lot.TableEntriesHAlign = 'center';
                lot.TableEntriesVAlign = 'center';
                append(chap,lot);
            end
        case 'jac'
            if nv==1
                T = varargin{1};
                chap = Chapter(Title='Spatial Jacobian',Numbered=false);
                chap.Layout.Landscape = true;
                append(chap,pipeline_table(T));
            end
        case 'scatnetEmphReg'
            if nv==1
                T = varargin{1};
                chap = Chapter(Title='ScatterNet for Emph from Registered Imspiratory CT',Numbered=false);
                chap.Layout.Landscape = true;
                append(chap,pipeline_table(T));
            end
        case 'dblood'
            if nv==1
                T = varargin{1};
                chap = Chapter(Title='Change in Blood Density',Numbered=false);
                chap.Layout.Landscape = true;
                append(chap,pipeline_table(T));
            end
        case 'prm10'
                T = varargin{1};
                chap = Chapter(Title='10-Color PRM Percents',Numbered=false);
                chap.Layout.Landscape = true;
                append(chap,pipeline_table(T));
        case 'prm'
                T = varargin{1};
                chap = Chapter(Title='Parametric Response Map (PRM) Percents',Numbered=false);
                chap.Layout.Landscape = true;
                append(chap,pipeline_table(T));
        case 'tprm'
                T = varargin{1};
                chap = Chapter(Title='Topological PRM (tPRM) Statistics',Numbered=false);
                chap.Layout.Landscape = true;
                append(chap,pipeline_table(T));
    end
    if ~isempty(chap)
        add(R,chap);
    end
end

function T = pipeline_table(T,colstyle_override)
% Create Table with desired formatting
    
    % Set ROI to row name and remove from table variables
    T.Properties.RowNames = T.ROI';
    T = removevars(T,'ROI');

    % Determine column number formatting
    colstyle = repmat({NumberFormat('%.2f')},1,size(T,2));
    vnames = T.Properties.VariableNames;
    if nargin==3
        for i = 1:size(colstyle_override,1)
            ind = ismember(vnames,colstyle_override{i,1});
            colstyle(ind) = {NumberFormat(colstyle_override{i,2})};
        end
    end

    % Create mlreportgen Table()
    T = MATLABTable(T);
    for i = 1:numel(colstyle)
        T.ColSpecGroups(i).Style = [T.ColSpecGroups(i).Style,colstyle(i)];
    end
    T.Style = {FontSize('10pt')};

