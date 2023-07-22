function plotxppaut4p4(varargin)
%% Plot bifurcation diagrams in Matlab that have been saved by XPPAUT
% How to use:
% In AUTO (XPPAUT) click-->File--->Write Pts
% and save your bifurcation diagram in a .dat file
% In Matlab run this program
% Select .dat file saved by AUTO
% A bifurcation diagram will be plotted.
% Tested with XPPAUT 8.0 and MATLAB 2017a.
% (c) Ting-Hao Hsu, 2018

%% Set Defaults
LineHeadings = {'Type','Class','LineStyle','Color','LineWidth',...
    'Marker','MarkerSize','MarkerEdgeColor','MarkerFaceColor','NrPts','Legend'};
colorRed = 'r'; colorBlue = 'b'; colorBlack = 'k'; colorNone = 'none';
colorGreen = 'g'; colorYellow = [.8 .8 0]; colorOlive = [0 .7 .7]; % colorCyan = 'c';
LineDefaultCell = {
    [1 0],'Stable Equilibria',...
        '-',colorRed,1,'none',6,colorRed,colorNone,5,'EP-';
    [2 0],'Unstable Equilibria',...
        '-',colorBlack,1,'none',6,colorBlack,colorNone,5,'EP+'; 
    [3 0],'Stable Periodic Orbits',...
        'none',colorGreen,1,'o',6,colorNone,colorGreen,5,'PO-';
    [4 0],'Unstable Periodic Orbits',...
        'none',colorBlue,1,'o',6,colorBlue,colorNone,5,'PO+';
    [2 5],'Branch Points',...
        '-',colorOlive,1,'none',6,colorOlive,colorNone,5,'BP';
    [2 3],'Hopft Bifurcation',...
        '-',colorBlue,1,'none',6,colorBlue,colorNone,5,'HB';
    [2 1],'Limit Points',...
        '-',colorRed,1,'none',6,colorRed,colorNone,5,'LP';
    [4 7],'Periodic Orbits',...
        '-',colorYellow,1,'none',6,colorYellow,colorNone,5,'PO';
    [2 6],'Periodic Doubling',...
        '-',colorOlive,1,'none',6,colorOlive,colorNone,5,'PD';
    [2 2],'Saddle-Node Limit Cycles',...
        '-',colorBlack,1,'none',6,colorBlack,colorNone,5,'SNLC'; 
    [NaN NaN],'Unknown Type',...
        '-',colorBlack,1,'none',6,colorBlack,colorNone,5,'??';
};

%% Initialize
LineDefaultTable = cell2table(LineDefaultCell);
LineDefaultTable.Properties.VariableNames = LineHeadings;
LineDefault = table2struct(LineDefaultTable);
TruncTolDefault = .1; % tolerence for truncation
LenTolDefault = 2; % tolerence for minimum length
LineStyleList = {'-','--',':','-.','none'};
MarkerStyleList = {'o','+','*','.','x','square','diamond',...
    '^','v','>','<','pentagram','hexagram','none'};
defaults = struct('LineDefault',LineDefault,...
    'TruncTol',TruncTolDefault,'LenTol',LenTolDefault,...
    'LineStyleList',{LineStyleList},'MarkerStyleList',{MarkerStyleList});
info = struct('Data',[],'Defaults',defaults,'Dir',pwd,...
    'LineData',[],'LineOptions',[],'Focus',[],...
    'TruncTol',[],'LenTol',[],'TruncTolOld',[],'LenTolOld',[],...
    'xlabel',[],'ylabel',[],'axlim',[],'hui',[]);
% Set Figure Units
set(0,'Units','characters')
tmp = get(0,'ScreenSize');
scrheight = tmp(4);
% Set figPlot
figPlot = figure('Units','characters','Visible','off');
pos0 = get(figPlot,'Position');
pos = [76 scrheight-pos0(4)-8.5 pos0(3:4)];
set(figPlot,'Position',pos,'Visible','on');
info.nfigPlot = get(figPlot,'Number');
% Set figPrompt
pos = [5 scrheight-45 70 40];
figPrompt = figure('Units','characters');
set(figPrompt,'MenuBar','none','ToolBar','none','NumberTitle','off');
str = 'Load to Start...';
set(figPrompt,'Position',pos,'Name',str);
info.nfigPrompt = get(figPrompt,'Number');
info.Defaults.posPrompt = pos;
% Call GUI
myPrompt(info,'START');
end

%% ====Draw Data====
function DrawData(info)
% Format Data
data = info.Data;
x = data(1,:);
y = data(2,:);
y2 = data(3,:);
styleA = data(4,:);
br = data(5,:);
styleB = data(6,:);
% Remove NaN
sd = sum(data);
ind = isnan(sd);
data(:,ind) = [];
info.Data = data;
% Detect the jumps of the curves
dx = abs(diff(x));
dy = abs(diff(y));
if isempty(info.TruncTol)
    info.TruncTol = info.Defaults.TruncTol;
end
TruncTol = info.TruncTol;
if isempty(info.LenTol)
    info.LenTol = info.Defaults.LenTol;
end
LenTol = info.LenTol;
tmp = find(dx>TruncTol | dy>TruncTol);
indj = [tmp numel(x)];
% Save Lines
dataOld = info.LineData;
info.LineData = [];
jline = 0;
ind0 = 1;
for ind2=indj
    dA = diff(styleA(ind0:ind2));
    dB = diff(styleB(ind0:ind2));
    dC = diff(br(ind0:ind2));
    tmp = ind0 - 1 + find(dA~=0 | dB~=0 | dC~=0);
    indc = [tmp ind2];
    for ind1=indc
        if ind1-ind0+1>=LenTol
            jline = jline + 1;
            xtmp = x(ind0:ind1);
            ytmp = y(ind0:ind1);
            y2tmp = y2(ind0:ind1);
            code = [styleA(ind1) br(ind1) styleB(ind1)];
            info.LineData(jline).x = xtmp;
            info.LineData(jline).y = ytmp;
            info.LineData(jline).y2 = y2tmp;
            info.LineData(jline).Code = code;
        end
        ind0 = ind1;
    end
    ind0 = ind2 + 1;
end
% Count Lines
numline = numel(info.LineData);
if numline>100
    str = sprintf('Total number of curves will be %d. Proceed anyway?',numline);
    stitle = 'Large Number of Curves';
    stroptions = {'Yes','No'};
    sdefault = 'No';
    choice  = questdlg(str,stitle,stroptions{:},sdefault);
    if strcmp(choice,'No')
        if ~isempty(info.TruncTolOld)
            info.TruncTol = info.TruncTolOld;
        end
        if ~isempty(info.LenTolOld)
            info.LenTol = info.LenTolOld;
        end
        info.Data = dataOld;
        myPrompt(info,'UPDATE');
        return;
    end
end
info.TruncTolOld = info.TruncTol;
info.LenTolOld = info.LenTol;
% Sort Lines
info.LineData = SortCodes(info.LineData);
% Set LineOptions and LineGroup
LineDefault = info.Defaults.LineDefault;
info.LineOptions = SetOption(info.LineData,LineDefault,info.LineOptions);
% Draw Lines
figPlot = figure(info.nfigPlot);
if isempty(info.axlim)
    clf(figPlot);
    plot(NaN,NaN);
    hold on;
    set(gca(figPlot),'FontSize',12);
else
    cla(figPlot);
    hold on;
end
for jline=1:numel(info.LineData)
    lineData = info.LineData(jline);
    code = info.LineData(jline).Code;
    type = code([1,3]);
    jop = find(arrayfun(@(c)isequal(c.Type,type),info.LineOptions),1);
    option = info.LineOptions(jop);
    x = lineData.x;
    y = lineData.y;
    y2 = lineData.y2;
    hp = plot(x,y,x,y2,'k-','linewidth',1);
    hp = SetMarkerIndices(x,y,y2,option,hp);
    info.LineData(jline).hp = hp;
    info.LineOptions(jop).jlines(end+1) = jline;
    headers = {'Color','LineStyle','LineWidth',...
        'Marker','MarkerSize','MarkerEdgeColor','MarkerFaceColor'}';
    for jhp=1:numel(hp)
        for jhd=1:numel(headers)
            ch = headers{jhd};
            hp(jhp).(ch) = option.(ch);
        end
    end
end
% Set Axes
axPlot = gca(figPlot);
xlim = get(axPlot,'XLim');
ylim = get(axPlot,'YLim');
if isempty(info.axlim)
    info.axlim = [xlim ylim];
end
set(axPlot,'XLim',xlim);
set(axPlot,'YLim',ylim);
if ~isempty(info.xlabel)
    xlabel(info.xlabel,'FontSize',16);
end
if ~isempty(info.ylabel)
    ylabel(info.ylabel,'FontSize',16);
end
% Make Legend
MakeLegend(info.LineOptions,info.LineData);
% Call GUI
myPrompt(info,'NEW');
end

%% ====SET MARKER INDICIES====
function hp=SetMarkerIndices(x,y,y2,option,hp)
NrPts = option.NrPts;
LineStyle = option.LineStyle;
MarkerIndices = 1:NrPts:numel(x);
if ~verLessThan('matlab','9.1') % MarkerIndices is only for MATLAB is 9.1(R2016b) or Newer
    for jhp=1:numel(hp)
        set(hp(jhp),'MarkerIndices',MarkerIndices);
        hp(jhp).MarkerIndices = MarkerIndices;
    end
elseif strcmp(LineStyle,'none')
    for jhp=1:numel(hp)
        delete(hp(jhp));
    end
    xtmp = x(MarkerIndices);
    ytmp = y(MarkerIndices);
    y2tmp = y2(MarkerIndices);
    hp = plot(xtmp,ytmp,xtmp,y2tmp,'ko','markersize',6);
end
end

%% ====SETOPTION====
function LineOptions = SetOption(LineData,LineDefault,LineOptionsOld)
tmp1 = arrayfun(@(c)num2str(c.Code([1,3])),LineData,'UniformOutput',false);
tmp2 = unique(tmp1);
types = arrayfun(@(j)str2double(strsplit(tmp2{j},' ')),...
    1:numel(tmp2),'UniformOutput',false);
numtypes = numel(types);
option = struct('Type',[],'jlines',[]);
LineOptions = repmat(option,1,numtypes);
for jop=1:numtypes
    type = types{jop};
    jOld = find(arrayfun(@(c)isequal(c.Type,type),LineOptionsOld),1);
    j0 = find(arrayfun(@(c)isequal(c.Type,type),LineDefault),1);
    if ~isempty(jOld)
        lnew = LineOptionsOld(jOld);
    elseif ~isempty(j0)
        lnew = LineDefault(j0);
    else
        lnew = LineDefault(end);
        lnew.Type = type;
    end
    lnew.jlines = [];
    for str = fieldnames(lnew)'
        LineOptions(jop).(str{1}) = lnew.(str{1});
    end
end
end

%% ====SORT CODES====
function lines = SortCodes(lines)
if ~isempty(lines)
    % Sort branch
    BR = arrayfun(@(c)c.Code(2),lines);
    [~,permC] = sort(BR);
    lines = lines(permC);
    % Sort the second code
    styleB = arrayfun(@(c)c.Code(3),lines);
    [~,permB] = sort(styleB);
    lines = lines(permB);
    % Sort the first code
    styleA = arrayfun(@(c)c.Code(1),lines);
    [~,permA] = sort(styleA);
    lines = lines(permA);
    % Count the number of visible lines
    numline = numel(lines);
    codepre = [];
    ncopy = 0;
    for jline=1:numline
        code = lines(jline).Code;
        if isequal(code,codepre)
            ncopy = ncopy+1;
        else
            ncopy = 1;
        end
        lines(jline).ncopy = ncopy;
        codepre = code;
    end
end
end

%% ====MAKE LEGEND====
function MakeLegend(options,lines)
numop = numel(options);
figs = NaN(1,numop);
names = cell(1,numop);
jfig = 0;
for jop=1:numel(options)
    jlines = options(jop).jlines;
    if ~isempty(jlines)
        jfig = jfig + 1;
        jline = options(jop).jlines(1);
        hp = lines(jline).hp;
        figs(jfig) = hp(1);
        names{jfig} = options(jop).Legend;
    end
end
if jfig>0
    figs(jfig+1:end) = [];
    names(jfig+1:end) = [];
    hlegend = legend(figs,names);
    set(hlegend,'FontSize',14);
end
end

%% ====MYPROMPT====
function myPrompt(info,status)
%% Preparation
MAXLINE = 23; % number of lines in each column
figPrompt = figure(info.nfigPrompt);
numOption = numel(info.LineOptions);
pos0 = get(figPrompt,'Position');
coord_top = pos0(4);
if isempty(info.Focus)
    info.Focus = 1;
end
%% Initialize
if strcmp(status,'START')
    clf(figPrompt);
    posbox = [pos0(3)/4 pos0(4)/2 30 4];
    str = 'Load File to Start...';
    uicontrol('Style','pushbutton','String',str,...
        'Units','characters','Position',posbox,...
        'Callback',@(s,e)setChoice(info,'Load'));
    %%-- Set Buttons
    posbox0 = [3 1.5 8 1.2];
    strs = {'Quit','Help'};
    for jstr=1:numel(strs)
        str = strs{jstr};
        posbox = posbox0 + (jstr>1)*[9 0 0 0];
        uicontrol('Style','pushbutton','String',str,...
            'Units','characters','Position',posbox,...
            'Callback',@(s,e)setChoice(info,str));
    end
    return;
end

if numel(info.LineOptions)==0
    return;
end

%% Resize and Show Filename
if strcmp(status,'NEW')
    clf(figPrompt);
    str = sprintf('Setup %s',info.Filename);
    set(figPrompt,'Name',str);
    posPrompt = get(figPrompt,'Position');
    posPromptDefalut = info.Defaults.posPrompt;
    posPrompt(3:4) = posPromptDefalut(3:4);
    set(figPrompt,'Position',posPrompt);
end

%% Set Buttons
jui = 0;
strs = {'Quit','Help','Reset','Black','Load/Append'};
opts = {'Quit','Help','Reset','Black','Load'};
if strcmp(status,'NEW')
    posbox0 = [3 1.5 8 1.2];
    indL = 5; % Load/Preview is the 5th botton
    for jh=1:numel(strs)
        jui = jui + 1;
        str = strs{jh};
        posbox = posbox0 + (jh-1)*[9 0 0 0] ...
            + (jh==indL)*[0 0 6 0] + (jh>indL)*[6 0 0 0];
        info.hui{jui} = uicontrol('Style','pushbutton','String',str,...
            'Units','characters','Position',posbox);
    end
elseif strcmp(status,'UPDATE')
    for jh=1:numel(opts)
        jui = jui + 1;
        opt = opts{jh};
        set(info.hui{jui},'Callback',@(s,e)setChoice(info,opt));
    end
end

%% ==Line Properties==
vspace = 0;
% Set Heading
if strcmp(status,'NEW')
    vspace = vspace + 2;
    strtxt = '== CURVE PROPERTIES ==';
    postxt = [3 coord_top-vspace 20 1];
    uicontrol('Style','text','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
end
% Select Line
jui = jui + 1;
focus = info.Focus;
if strcmp(status,'NEW')
    vspace = vspace + 1.5;
    strbox = cell(1,numOption);
    for j=1:numOption
        prop = info.LineOptions(j);
        if strcmp(prop.Class,'Unknown Type')
            strbox{j} = sprintf('(%s) %s [%d,*,%d]',prop.Legend,prop.Class,prop.Type);
        else
            strbox{j} = sprintf('(%s) %s',prop.Legend,prop.Class);
        end
    end
    posbox = [2 coord_top-vspace 30 1];
    info.hui{jui} = uicontrol('Style','popup','String',strbox,'Value',focus,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    set(info.hui{jui},'Value',focus,...
        'Callback',@(s,e)setChoice(info,'SelectLine',s.Value));
end
% Preparation
if isempty(info.Focus)
    info.Focus = 1;
end
option = info.LineOptions(focus);
% Set LineStyle
jui = jui + 1;
list = info.Defaults.LineStyleList;
val = find(strcmp(list,option.LineStyle),1);
if strcmp(status,'NEW')
    vspace = vspace + 2;
    postxt = [3 coord_top-vspace 12 1];
    strtxt = 'LineStyle';
    uicontrol('Style','text','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
    posbox = [20 coord_top-vspace 19 1.2] - [1 0 0 0];
    info.hui{jui} = uicontrol('Style','popup','String',list,'Value',val,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    set(info.hui{jui},'String',list,'Value',val,...
        'Callback',@(s,e)setChoice(info,'LineStyle',list,s.Value));
end
% Set Color
color = option.Color;
if ~ischar(color)
    strbox = sprintf('%.1f, %.1f, %.1f',color);
else
    strbox = color;
end
if strcmp(status,'NEW')
    vspace = vspace + 2;
    postxt = [3 coord_top-vspace 16 1];
    strtxt = 'Color';
    uicontrol('Style','text','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
    posbox = [20 coord_top-vspace 6 1.2];
    strtxt = 'RGB';
    jui = jui + 1;
    info.hui{jui} = uicontrol('Style','pushbutton','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
    posbox = [27 coord_top-vspace 9 1.2];
    jui = jui + 1;
    info.hui{jui} = uicontrol('Style','edit','String',strbox,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    jui = jui + 1;
    set(info.hui{jui},...
        'Callback',@(s,e)setChoice(info,'Color',color,'RGB'));
    jui = jui + 1;
    set(info.hui{jui},'String',strbox,...
        'Callback',@(s,e)setChoice(info,'Color',s.String));
end
% Set LineWidth
jui = jui + 1;
strbox = num2str(option.LineWidth);
if strcmp(status,'NEW')
    vspace = vspace + 2;
    postxt = [3 coord_top-vspace 16 1];
    uicontrol('Style','text','String','LineWidth',...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
    posbox = [20 coord_top-vspace 7 1.2];
    info.hui{jui} = uicontrol('Style','edit','String',strbox,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    set(info.hui{jui},'String',strbox,...
        'Callback',@(s,e)setChoice(info,'LineWidth',s.String));
end
% Set Marker
jui = jui + 1;
list = info.Defaults.MarkerStyleList;
val = find(strcmp(list,option.Marker),1);
if strcmp(status,'NEW')
    vspace = vspace + 2;
    postxt = [3 coord_top-vspace 12 1];
    uicontrol('Style','text','String','Marker',...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
    posbox = [20 coord_top-vspace 19 1.2] - [1 0 0 0];
    info.hui{jui} = uicontrol('Style','popup','String',list,'Value',val,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    set(info.hui{jui},'Value',val,...
        'Callback',@(s,e)setChoice(info,'Marker',list,s.Value));
end
% Set MarkerEdgeColor
color = option.MarkerEdgeColor;
if ~ischar(color)
    strbox = sprintf('%.1f, %.1f, %.1f',color);
else
    strbox = color;
end
if strcmp(status,'NEW')
    vspace = vspace + 2;
    postxt = [3 coord_top-vspace 17 1];
    uicontrol('Style','text','String','MarkerEdgeColor',...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
    posbox = [20 coord_top-vspace 6 1.2];
    jui = jui + 1;
    info.hui{jui} = uicontrol('Style','pushbutton','String','RGB',...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
    posbox = [27 coord_top-vspace 9 1.2];
    jui = jui + 1;
    info.hui{jui} = uicontrol('Style','edit','String',[],...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    jui = jui + 1;
    set(info.hui{jui},...
        'Callback',@(s,e)setChoice(info,'MarkerEdgeColor',color,'RGB'));
    jui = jui + 1;
    set(info.hui{jui},'String',strbox,...
        'Callback',@(s,e)setChoice(info,'MarkerEdgeColor',s.String));
end
% Set MarkerFaceColor
color = option.MarkerFaceColor;
if ~ischar(color)
    strbox = sprintf('%.1f, %.1f, %.1f',color);
else
    strbox = color;
end
if strcmp(status,'NEW')
    vspace = vspace + 2;
    postxt = [3 coord_top-vspace 17 1];
    uicontrol('Style','text','String','MarkerFaceColor',...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
    posbox = [20 coord_top-vspace 6 1.2];
    jui = jui + 1;
    info.hui{jui} = uicontrol('Style','pushbutton','String','RGB',...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
    posbox = [27 coord_top-vspace 9 1.2];
    jui = jui + 1;
    info.hui{jui} = uicontrol('Style','edit','String',strbox,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    jui = jui + 1;
    set(info.hui{jui},...
        'Callback',@(s,e)setChoice(info,'MarkerFaceColor',color,'RGB'));
    jui = jui + 1;
    set(info.hui{jui},'String',strbox,...
        'Callback',@(s,e)setChoice(info,'MarkerFaceColor',s.String));
end
% Set MarkerSize
jui = jui + 1;
strbox = num2str(option.MarkerSize);
if strcmp(status,'NEW')
    vspace = vspace + 2;
    postxt = [3 coord_top-vspace 12 1];
    uicontrol('Style','text','String','MarkerSize',...
        'HorizontalAlignment','left','Units','characters','Position',postxt)
    posbox = [20 coord_top-vspace 7 1.2];
    info.hui{jui} = uicontrol('Style','edit','String',strbox,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    set(info.hui{jui},'String',strbox,...
        'Callback',@(s,e)setChoice(info,'MarkerSize',s.String));
end
% Set NrPts
jui = jui + 1;
strbox = num2str(option.NrPts);
if strcmp(status,'NEW')
    vspace = vspace + 2;
    postxt = [3 coord_top-vspace 12 1];
    strtxt = 'MarkerSpace';
    uicontrol('Style','text','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',postxt)
    posbox = [20 coord_top-vspace 7 1.2];
    info.hui{jui} = uicontrol('Style','edit','String',strbox,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    set(info.hui{jui},'String',strbox,...
        'Callback',@(s,e)setChoice(info,'NrPts',s.String));
end
% Set Legend
jui = jui + 1;
strbox = option.Legend;
if strcmp(status,'NEW')
    vspace = vspace + 2;
    postxt = [3 coord_top-vspace 12 1];
    strtxt = 'Legend';
    uicontrol('Style','text','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',postxt)
    posbox = [20 coord_top-vspace 16 1.2];
    info.hui{jui} = uicontrol('Style','edit','String',strbox,...
        'HorizontalAlignment','left','Units','characters','Position',posbox,...
        'Callback',@(s,e)setChoice(info,'Legend',s.String));
elseif strcmp(status,'UPDATE')
    set(info.hui{jui},'String',strbox,...
        'Callback',@(s,e)setChoice(info,'Legend',s.String));
end
%% ==General Settings==
% Set Heading
if strcmp(status,'NEW')
    vspace = vspace + 2;
    strtxt = '== GENERAL SETTINGS ==';
    postxt = [3 coord_top-vspace 20 1];
    uicontrol('Style','text','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
end
% Set x-label
jui = jui + 1;
strbox = info.xlabel;
if strcmp(status,'NEW')
    vspace = vspace + 1.5;
    postxt = [3 coord_top-vspace 12 1];
    posbox = [20 coord_top-vspace 16 1.2];
    strtxt = 'x-label';
    uicontrol('Style','text','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
    info.hui{jui} = uicontrol('Style','edit','String',strbox,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    set(info.hui{jui},'String',strbox,...
        'Callback',@(s,e)setChoice(info,'xlabel',s.String));
end
% Set y-label
jui = jui + 1;
    strbox = info.ylabel;
if strcmp(status,'NEW')
    vspace = vspace + 1.5;
    postxt = [3 coord_top-vspace 12 1];
    posbox = [20 coord_top-vspace 16 1.2];
    strtxt = 'y-label';
    uicontrol('Style','text','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
    info.hui{jui} = uicontrol('Style','edit','String',strbox,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    set(info.hui{jui},'String',strbox,...
        'Callback',@(s,e)setChoice(info,'ylabel',s.String));
end
% Set axix-limits
hd = {'xLo','xHi','yLo','yHi'};
vspace = vspace + 1.5;
postxt = [3 coord_top-vspace 4 1];
posbox = [8 coord_top-vspace 8 1.2];
for jax=1:numel(hd)
    strtxt = hd{jax};
    strbox = sprintf('%.2f',info.axlim(jax));
    jui = jui + 1;
    if strcmp(status,'NEW')
        postxt = postxt + (jax==2 || jax==4)*[17 0 0 0] + (jax==3)*[-17 -1.5 0 0];
        uicontrol('Style','text','String',strtxt,...
            'HorizontalAlignment','left','Units','characters','Position',postxt);
        posbox = posbox + (jax==2 || jax==4)*[17 0 0 0] + (jax==3)*[-17 -1.5 0 0];
        info.hui{jui} = uicontrol('Style','edit','String',strbox,...
            'HorizontalAlignment','left','Units','characters','Position',posbox);
    elseif strcmp(status,'UPDATE')
        set(info.hui{jui},'String',strbox,...
            'Callback',@(s,e)setChoice(info,'axlim',s.String,jax));
    end
end
%% ==Numerics==
% Heading
if strcmp(status,'NEW')
    vspace = vspace + 3.5;
    strtxt = '== NUMERICS ==';
    postxt = [3 coord_top-vspace 20 1];
    uicontrol('Style','text','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
end
% Set TruncTol
jui = jui + 1;
strbox = sprintf('%.4f',info.TruncTol);
if strcmp(status,'NEW')
    vspace = vspace + 1.5;
    postxt = [3 coord_top-vspace 22 1];
    strtxt = 'TruncTol (stray lines)';
    uicontrol('Style','text','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
    posbox = [25 coord_top-vspace 11 1.2];
    info.hui{jui} = uicontrol('Style','edit','String',strbox,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    set(info.hui{jui},'String',strbox,...
        'Callback',@(s,e)setChoice(info,'TruncTol',s.String)); 
end
% Set TruncTol
jui = jui + 1;
strbox = sprintf('%d',info.LenTol);
if strcmp(status,'NEW')
    vspace = vspace + 1.5;
    postxt = [3 coord_top-vspace 26 1];
    strtxt = 'LenTol (short pieces)';
    uicontrol('Style','text','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
    posbox = [29 coord_top-vspace 7 1.2];
    info.hui{jui} = uicontrol('Style','edit','String',strbox,...
        'HorizontalAlignment','left','Units','characters','Position',posbox);
elseif strcmp(status,'UPDATE')
    set(info.hui{jui},'String',strbox,...
        'Callback',@(s,e)setChoice(info,'LenTol',s.String)); 
end
%% ==Visibility==
vspace = 0;
% Set Heading
if strcmp(status,'NEW')
    vspace = vspace + 2;
    strtxt = '== CURVE VISIBILITY ==';
    postxt = [45 coord_top-vspace 20 1];
    uicontrol('Style','text','String',strtxt,...
        'HorizontalAlignment','left','Units','characters','Position',postxt);
end
% Set Checkbox
vspace = vspace + 1.2;
posbox0 = [45 coord_top-vspace 24 1];
numlines = numel(info.LineData);
for jline=1:numlines
    jui = jui + 1;
    lines = info.LineData(jline);
    type = lines.Code([1 3]);
    jop = find(arrayfun(@(c)isequal(c.Type,type),info.LineOptions),1);
    op = info.LineOptions(jop);
    ncopy = lines.ncopy;
    if ncopy==1
        strbox = sprintf('[%d,%d,%d] (%s)',lines.Code,op.Legend);
    else
        strbox = sprintf('[%d,%d,%d] (%s)(%d)',lines.Code,op.Legend,ncopy-1);
    end
    hp = lines.hp;
    val = strcmp(hp(1).Visible,'on');
    if strcmp(status,'NEW')
        mjline = mod(jline-1,MAXLINE);
        if (mjline==0 && jline>1)
            posbox0 = posbox0 + [25 0 0 0];
            posPrompt = posPrompt + [0 0 25 0];
            set(figPrompt,'Position',posPrompt);
        end
        posbox = posbox0 + mjline*[0 -1.5 0 0];
        info.hui{jui} = uicontrol('Style','checkbox','Value',val,'String',strbox,...
            'Units','characters','Position',posbox);
    elseif strcmp(status,'UPDATE')
        set(info.hui{jui},'Value',val,...
            'Callback',@(s,e)setChoice(info,'Visible',jline,s.Value));
    end
end
if strcmp(status,'NEW')
    myPrompt(info,'UPDATE');
end
end

%% ====SETCHOICE====
function setChoice(info,opt,varargin)
figure(info.nfigPlot);
%%-- Update Properties
focus = info.Focus;
switch opt
    case 'axlim'
        str = varargin{1};
        val = str2double(str);
        jax = varargin{2};
        info.axlim(jax) = val;
        axis(info.axlim);
        myPrompt(info,'UPDATE');
        return;
    case {'xlabel','ylabel'}
        str = varargin{1};
        feval(opt,str,'FontSize',16);
        info.(opt) = str;
        myPrompt(info,'UPDATE');
        return;
    case 'Visible'
        jline = varargin{1};
        val = varargin{2};
        hp = info.LineData(jline).hp;
        type = info.LineData(jline).Code([1,3]);
        jop = find(arrayfun(@(c)isequal(c.Type,type),info.LineOptions),1);
        jlines = info.LineOptions(jop).jlines;
        if val
            for jfig=1:numel(hp)
                fig = hp(jfig);
                set(fig,'Visible','on');
            end
            jlines = [jlines jline];
        else
            for jfig=1:numel(hp)
                fig = hp(jfig);
                set(fig,'Visible','off');
            end
            jlines = setdiff(jlines,jline);
        end
        info.LineOptions(jop).jlines = jlines;
        MakeLegend(info.LineOptions,info.LineData);
        myPrompt(info,'UPDATE');
        return;
    case 'SelectLine'
        val = varargin{1};
        info.Focus = val;
        myPrompt(info,'UPDATE');
        return;
    case {'LineStyle','Marker'}
        list = varargin{1};
        val = varargin{2};
        str = list{val};
        info.LineOptions(focus).(opt) = str;
        type = info.LineOptions(focus).Type;
        jlines = find(arrayfun(@(c)isequal(type,c.Code([1,3])),info.LineData));
        for jline=jlines
            hp = info.LineData(jline).hp;
            for jfig=1:numel(hp)
                fig = hp(jfig);
                set(fig,opt,str);
            end
        end
        myPrompt(info,'UPDATE');
        return;
    case {'LineWidth','MarkerSize'}
        vin = varargin{1};
        val = str2double(vin);
        info.LineOptions(focus).(opt) = val;
        type = info.LineOptions(focus).Type;
        jlines = find(arrayfun(@(c)isequal(type,c.Code([1,3])),info.LineData));
        for jline=jlines
            hp = info.LineData(jline).hp;
            for jfig=1:numel(hp)
                fig = hp(jfig);
                set(fig,opt,val);
            end
        end
        myPrompt(info,'UPDATE');
        return;
    case 'NrPts'
        vin = varargin{1};
        val = str2double(vin);
        info.LineOptions(focus).(opt) = val;
        option = info.LineOptions(focus);
        type = info.LineOptions(focus).Type;
        jlines = find(arrayfun(@(c)isequal(type,c.Code([1,3])),info.LineData));
        for jline=jlines
            line = info.LineData(jline);
            hp = SetMarkerIndices(line.x,line.y,line.y2,option,line.hp);
            info.LineData(jline).hp = hp;
        end
        myPrompt(info,'UPDATE');
        return;
    case {'Color','MarkerEdgeColor','MarkerFaceColor'}
        vin = varargin{1};
        if numel(varargin)>1 % RGB
            color = vin;
            cnew = uisetcolor(color);
        elseif ~sum(isstrprop(vin,'digit'))
            cnew = vin;
        else
            for ch={',','[',']'}
                vin = strrep(vin,ch{1},' ');
            end
            strs = strsplit(vin,' ');
            ind = ~cellfun('isempty',strs);
            strs = strs(ind);
            cnew = str2double(strs);
        end
        info.LineOptions(focus).(opt) = cnew;
        type = info.LineOptions(focus).Type;
        jlines = find(arrayfun(@(c)isequal(type,c.Code([1,3])),info.LineData));
        for jline=jlines
            hp = info.LineData(jline).hp;
            for jfig=1:numel(hp)
                fig = hp(jfig);
                set(fig,opt,cnew);
            end
        end
        myPrompt(info,'UPDATE');
        return;
    case 'Legend'
        str = varargin{1};
        info.LineOptions(focus).(opt) = str;
        MakeLegend(info.LineOptions,info.LineData);
        myPrompt(info,'UPDATE');
        return;
    case {'TruncTol','LenTol'}
        % Save old value
        valOld = info.(opt);
        optOld = [opt,'Old'];
        info.(optOld) = valOld;
        % Check input format
        str = varargin{1};
        tflag = 0;
        if strcmp(opt,'TruncTol')
            val = str2double(str);
            if ~(val>0)
                strtxt = sprintf('%s must be a positive number.',opt);
                tflag = 1;
            end
        elseif strcmp(opt,'LenTol')
            str = strrep(str,' ','');
            val = round(str2double(str));
            if ~prod(isstrprop(str,'digit')) || ~(val>0)
                strtxt = sprintf('%s must be a positive integer.',opt);
                strtitle = '';
                msgbox(strtxt,strtitle,'modal');
                myPrompt(info,'UPDATE'); % ignore the change
                return;
            end
        end
        if tflag~=0
            strtitle = '';
            msgbox(strtxt,strtitle,'modal');
            myPrompt(info,'UPDATE'); % ignore the change
            return;
        end
        % Update value
        info.(opt) = val;
        DrawData(info);
        return;
    case 'Quit'
        msgstr = 'Quit the program?';
        choice = questdlg(msgstr,'','Yes','No','Yes');
        if strcmp(choice,'No')
            return;
        end
        close([info.nfigPrompt info.nfigPlot]);
        return;
    case 'Load'
        str = '.dat file saved by AUTO (XPPAUT)';
        [filename,path] = uigetfile(fullfile(info.Dir,'*.dat'),str);
        if isequal(filename,0)
            return;
        end
        info.Dir = path;
        fp = fopen(fullfile(path,filename),'r');
        data = fscanf(fp,'%f',[6,inf]);
        fclose(fp);
        if ~isempty(info.Data)
            fkeeps = {'Data','Defaults','nfigPlot','nfigPrompt','Dir'};
            stitle = 'Loading';
            str = 'Erase or append to the current diagram?';
            stropt = {'Erase','Append','Cancel'};
            sdefault = 'Erase';
            choice  = questdlg(str,stitle,stropt{:},sdefault);
            if strcmp(choice,'Append')
                fkeeps = [fkeeps,{'xlabel','ylabel','LineOptions'}];
                data = [info.Data data];
            elseif strcmp(choice,'Cancel')
                myPrompt(info,'UPDATE');
                return;
            end
            fnames = fieldnames(info)';
            fdel = setdiff(fnames,fkeeps);
            for ch=fdel
                info.(ch{1}) = [];
            end
        end
        info.Data = data;
        info.Filename = filename;
        DrawData(info);
        return;
    case 'Help'
        doc();
        return;
    case 'Reset'
        fnames = fieldnames(info)';
        fkeeps = {'Data','Defaults','nfigPlot','nfigPrompt','Filename','Dir'};
        fdel = setdiff(fnames,fkeeps);
        for ch=fdel
            info.(ch{1}) = [];
        end
        DrawData(info);
        return;
    case 'Black'
        colorBlack = 'k';
        colorNone = 'none';
        for hd={'Color','MarkerFaceColor','MarkerEdgeColor'}
            field = hd{1};
            for jop=1:numel(info.LineOptions)
                color = info.LineOptions(jop).(field);
                if ~isequal(color,colorNone)
                    info.LineOptions(jop).(field) = colorBlack;
                end
            end
            for jline=1:numel(info.LineData)
                hp = info.LineData(jline).hp;
                for jfig=1:numel(hp)
                    fig = hp(jfig);
                    color = fig.(field);
                    if ~isequal(color,colorNone)
                        fig.(field) = colorBlack;
                    end
                end
            end
        end
        % Dashed Lines for Unstable Equilibria
        type0 = [2,0];
        for jop=find(arrayfun(@(c)isequal(c.Type,type0),info.LineOptions))
            info.LineOptions(jop).LineStyle = '--';
        end
        for jline=find(arrayfun(@(c)isequal(c.Code([1,3]),type0),info.LineData))
            hp = info.LineData(jline).hp;
            for jfig=1:numel(hp)
                fig = hp(jfig);
                set(fig,'LineStyle', '--');
            end
        end
        myPrompt(info,'UPDATE');
        return;
end
end

%% ====DOC====
function doc()
    str_title = 'Instruction';
    str_doc = {
        '==HOW TO USE==';
        'In AUTO (XPPAUT) click File-->Write Pts';
        'and save your bifurcation diagram in a .dat file.';
        'In Matlab run this program';
        'and select a .dat file saved by AUTO.';
        'A bifurcation diagram will be plotted.';
        '';
        '==CURVE PROPERTIES==';
        'Select from the pop-up menus to choose a curve.';
        'Click RGB to choose a color.';
        'MarkerSpace:';
        '  Number of data points between every two markers.';
        '';
        '==GENERAL SETTINGS==';
        'Location of legend can be moved by dragging or';
        'by right-click using a mouse.';
        '';
        '==NUMERICS==';
        'TruncTol (Minimum distance of two consecuitive points):';
        '  Increase to avoid undesired gaps or missing curves.';
        '  Decrease to remove undesired stray lines.';
        'LenTol (Minimum number of data points in a curve):';
        '  Increase to remove curves that are too short.';
        '  Decrease to avoid missing curves.';
        '';
        '==CURVE VISIBILITY==';
        'Curves are listed by their codes.';
        'Check boxes:';
        '  Un-check to hide the curve.';
        '  Re-check to show the curve.';
        ''
    };
    msgbox(str_doc,str_title);
end