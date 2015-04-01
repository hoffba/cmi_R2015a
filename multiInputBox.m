function answer = multiInputBox(s)
% Creates UI box for input of linked variables
% Listbox on top linked to variable inputbox text
% INPUT:    s - structure containing default values,
%               first field contains cell array of (N) listbox options,
%               other fields contain either: single default value
%                                         or cell array of (N) values.
% OUTPUT:   answer - cell array of resulting user inputs

fields = fieldnames(s);
N = length(s.(fields{1}));
% First make sure structure s is in correct format
ok = iscell(s.(fields{1}));
for i = 1:length(fields)
    ok = ok && (isnumeric(s.(fields{i})) || ...
        (iscell(s.(fields{i})) && length(s.(fields{i}))==N));
end
if ok
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Create InputFig %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    FigWidth=175;
    FigHeight=100;
    FigPos(3:4)=[FigWidth FigHeight];  %#ok
    FigColor=get(0,'DefaultUicontrolBackgroundColor');

    InputFig=dialog(                     ...
      'Visible'          ,'off'      , ...
      'KeyPressFcn'      ,@doFigureKeyPress, ...
      'Name'             ,Title      , ...
      'Pointer'          ,'arrow'    , ...
      'Units'            ,'pixels'   , ...
      'UserData'         ,'Cancel'   , ...
      'Tag'              ,Title      , ...
      'HandleVisibility' ,'callback' , ...
      'Color'            ,FigColor   , ...
      'NextPlot'         ,'add'      , ...
      'WindowStyle'      ,WindowStyle, ...
      'Resize'           ,Resize       ...
      );


    %%%%%%%%%%%%%%%%%%%%%
    %%% Set Positions %%%
    %%%%%%%%%%%%%%%%%%%%%
    DefOffset    = 5;
    DefBtnWidth  = 53;
    DefBtnHeight = 23;

    TextInfo.Units              = 'pixels'   ;
    TextInfo.FontSize           = get(0,'FactoryUicontrolFontSize');
    TextInfo.FontWeight         = get(InputFig,'DefaultTextFontWeight');
    TextInfo.HorizontalAlignment= 'left'     ;
    TextInfo.HandleVisibility   = 'callback' ;

    StInfo=TextInfo;
    StInfo.Style              = 'text'  ;
    StInfo.BackgroundColor    = FigColor;


    EdInfo=StInfo;
    EdInfo.FontWeight      = get(InputFig,'DefaultUicontrolFontWeight');
    EdInfo.Style           = 'edit';
    EdInfo.BackgroundColor = 'white';

    BtnInfo=StInfo;
    BtnInfo.FontWeight          = get(InputFig,'DefaultUicontrolFontWeight');
    BtnInfo.Style               = 'pushbutton';
    BtnInfo.HorizontalAlignment = 'center';

    % Add VerticalAlignment here as it is not applicable to the above.
    TextInfo.VerticalAlignment  = 'bottom';
    TextInfo.Color              = get(0,'FactoryUicontrolForegroundColor');


    % adjust button height and width
    btnMargin=1.4;
    ExtControl=uicontrol(InputFig   ,BtnInfo     , ...
      'String'   ,getString(message('MATLAB:uistring:popupdialogs:Cancel'))        , ...
      'Visible'  ,'off'         ...
      );

    % BtnYOffset  = DefOffset;
    BtnExtent = get(ExtControl,'Extent');
    BtnWidth  = max(DefBtnWidth,BtnExtent(3)+8);
    BtnHeight = max(DefBtnHeight,BtnExtent(4)*btnMargin);
    delete(ExtControl);

    % Determine # of lines for all Prompts
    TxtWidth=FigWidth-2*DefOffset;
    ExtControl=uicontrol(InputFig   ,StInfo     , ...
      'String'   ,''         , ...
      'Position' ,[ DefOffset DefOffset 0.96*TxtWidth BtnHeight ] , ...
      'Visible'  ,'off'        ...
      );

    WrapQuest=cell(NumQuest,1);
    QuestPos=zeros(NumQuest,4);

    for ExtLp=1:NumQuest
      if size(NumLines,2)==2
        [WrapQuest{ExtLp},QuestPos(ExtLp,1:4)]= ...
          textwrap(ExtControl,Prompt(ExtLp),NumLines(ExtLp,2));
      else
        [WrapQuest{ExtLp},QuestPos(ExtLp,1:4)]= ...
          textwrap(ExtControl,Prompt(ExtLp),80);
      end
    end % for ExtLp

    delete(ExtControl);
    QuestWidth =QuestPos(:,3);
    QuestHeight=QuestPos(:,4);
    if ismac % Change Edit box height to avoid clipping on mac.
        editBoxHeightScalingFactor = 1.4;
    else 
        editBoxHeightScalingFactor = 1;
    end
    TxtHeight=QuestHeight(1)/size(WrapQuest{1,1},1) * editBoxHeightScalingFactor;
    EditHeight=TxtHeight*NumLines(:,1);
    EditHeight(NumLines(:,1)==1)=EditHeight(NumLines(:,1)==1)+4;

    FigHeight=(NumQuest+2)*DefOffset    + ...
      BtnHeight+sum(EditHeight) + ...
      sum(QuestHeight);

    TxtXOffset=DefOffset;

    QuestYOffset=zeros(NumQuest,1);
    EditYOffset=zeros(NumQuest,1);
    QuestYOffset(1)=FigHeight-DefOffset-QuestHeight(1);
    EditYOffset(1)=QuestYOffset(1)-EditHeight(1);

    for YOffLp=2:NumQuest,
      QuestYOffset(YOffLp)=EditYOffset(YOffLp-1)-QuestHeight(YOffLp)-DefOffset;
      EditYOffset(YOffLp)=QuestYOffset(YOffLp)-EditHeight(YOffLp);
    end % for YOffLp

    QuestHandle=[];
    EditHandle=[];

    AxesHandle=axes('Parent',InputFig,'Position',[0 0 1 1],'Visible','off');

    inputWidthSpecified = false;

    for lp=1:NumQuest,
      if ~ischar(DefAns{lp}),
        delete(InputFig);
        error(message('MATLAB:inputdlg:InvalidInput'));
      end


      EditHandle(lp)=uicontrol(InputFig    , ...
        EdInfo      , ...
        'Max'        ,NumLines(lp,1)       , ...
        'Position'   ,[ TxtXOffset EditYOffset(lp) TxtWidth EditHeight(lp)], ...
        'String'     ,DefAns{lp}           , ...
        'Tag'        ,'Edit'                 ...
        );

      QuestHandle(lp)=text('Parent'     ,AxesHandle, ...
        TextInfo     , ...
        'Position'   ,[ TxtXOffset QuestYOffset(lp)], ...
        'String'     ,WrapQuest{lp}                 , ...
        'Interpreter',Interpreter                   , ...
        'Tag'        ,'Quest'                         ...
        );

      MinWidth = max(QuestWidth(:));
      if (size(NumLines,2) == 2)
        % input field width has been specified.
        inputWidthSpecified = true;
        EditWidth = setcolumnwidth(EditHandle(lp), NumLines(lp,1), NumLines(lp,2));
        MinWidth = max(MinWidth, EditWidth);
      end
      FigWidth=max(FigWidth, MinWidth+2*DefOffset);

    end % for lp

    % fig width may have changed, update the edit fields if they dont have user specified widths.
    if ~inputWidthSpecified
      TxtWidth=FigWidth-2*DefOffset;
      for lp=1:NumQuest
        set(EditHandle(lp), 'Position', [TxtXOffset EditYOffset(lp) TxtWidth EditHeight(lp)]);
      end
    end

    FigPos=get(InputFig,'Position');

    FigWidth=max(FigWidth,2*(BtnWidth+DefOffset)+DefOffset);
    FigPos(1)=0;
    FigPos(2)=0;
    FigPos(3)=FigWidth;
    FigPos(4)=FigHeight;

    set(InputFig,'Position',getnicedialoglocation(FigPos,get(InputFig,'Units')));

    OKHandle=uicontrol(InputFig     ,              ...
      BtnInfo      , ...
      'Position'   ,[ FigWidth-2*BtnWidth-2*DefOffset DefOffset BtnWidth BtnHeight ] , ...
      'KeyPressFcn',@doControlKeyPress , ...
      'String'     ,getString(message('MATLAB:uistring:popupdialogs:OK'))        , ...
      'Callback'   ,@doCallback , ...
      'Tag'        ,'OK'        , ...
      'UserData'   ,'OK'          ...
      );

    setdefaultbutton(InputFig, OKHandle);

    CancelHandle=uicontrol(InputFig     ,              ...
      BtnInfo      , ...
      'Position'   ,[ FigWidth-BtnWidth-DefOffset DefOffset BtnWidth BtnHeight ]           , ...
      'KeyPressFcn',@doControlKeyPress            , ...
      'String'     ,getString(message('MATLAB:uistring:popupdialogs:Cancel'))    , ...
      'Callback'   ,@doCallback , ...
      'Tag'        ,'Cancel'    , ...
      'UserData'   ,'Cancel'       ...
      ); %#ok

    handles = guihandles(InputFig);
    handles.MinFigWidth = FigWidth;
    handles.FigHeight   = FigHeight;
    handles.TextMargin  = 2*DefOffset;
    guidata(InputFig,handles);
    set(InputFig,'ResizeFcn', {@doResize, inputWidthSpecified});

    % make sure we are on screen
    movegui(InputFig)

    % if there is a figure out there and it's modal, we need to be modal too
    if ~isempty(gcbf) && strcmp(get(gcbf,'WindowStyle'),'modal')
      set(InputFig,'WindowStyle','modal');
    end

    set(InputFig,'Visible','on');
    drawnow;

    if ~isempty(EditHandle)
      uicontrol(EditHandle(1));
    end

    if ishghandle(InputFig)
      % Go into uiwait if the figure handle is still valid.
      % This is mostly the case during regular use.
      uiwait(InputFig);
    end

    % Check handle validity again since we may be out of uiwait because the
    % figure was deleted.
    if ishghandle(InputFig)
      answer={};
      if strcmp(get(InputFig,'UserData'),'OK'),
        answer=cell(NumQuest,1);
        for lp=1:NumQuest,
          answer(lp)=get(EditHandle(lp),{'String'});
        end
      end
      delete(InputFig);
    else
      answer={};
    end
end