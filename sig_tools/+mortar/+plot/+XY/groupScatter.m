function [h, lh, gn] = groupScatter(x,y,g,varargin)
%groupScatter   Scatter plot with grouping variable
%   GROUPSCATTER(X,Y,G) creates a scatter plot of the vectors X and Y grouped
%   by G.  Points with the same value of G are shown with the same color
%   and marker.  G is a grouping variable defined as a categorical
%   variable, vector, cell array of strings, or string matrix, and it must
%   have the same number of rows as X and Y.  Alternatively G can be a cell
%   array of grouping variables (such as {G1 G2 G3}) to group the values in
%   X by each unique combination of grouping variable values.  Use the data
%   cursor to read precise values and observation numbers from the plot.
%
%   GROUPSCATTER(X,Y,G,CLR,SYM,SIZ) specifies the colors, markers, and
%   size to use.  CLR is either a string of color specifications or
%   a three-column matrix of color specifications.  SYM is a string
%   of marker specifications.  Type "help plot" for more information.
%   For example, if SYM='o+x', the first group will be plotted with a
%   circle, the second with plus, and the third with x.  SIZ is a
%   marker size to use for all plots.  By default, the marker is '.'.
%
%   GROUPSCATTER(X,Y,G,CLR,SYM,SIZ,DOLEG) lets you control whether legends
%   are created.  Set DOLEG to 'on' (default) or 'off'.
%
%   GROUPSCATTER(X,Y,G,CLR,SYM,SIZ,DOLEG,XNAM,YNAM) specifies XNAM and
%   YNAM as the names of the X and Y variables.  Each must be a
%   character string.  If you omit XNAM and YNAM, GSCATTER attempts to
%   determine the names of the variables passed in as the first and
%   second arguments.
%
%   H = GROUPSCATTER(...) returns an array of handles to the objects
%   created.
%
%   Example:  Scatter plot of car data coded by country.
%      load carsmall
%      gpscatter(Weight, MPG, Origin)
%
%   See also GRPSTATS, GRP2IDX.

%   Copyright 1993-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2010/03/22 04:41:24 $

pnames = {'--color', '--symbol', '--size',...
    '--show_legend', '--xlabel', '--ylabel',...
    '--legend_location', '--fill_marker', '--palette'};
dflts = {'', '.', [],...
    true, inputname(1), inputname(2),...
    'best', false, ''};
arg_help = {'Marker color as either a string of color specifications or a three-column matrix of color specifications',...
            'Marker symbols as a string of marker specifications',...
            'Marker size that can be scalar value that is applied globally or a vector of the same size as X that is applied per data point',...
            'Display legend if true',...
            'X-axis label',...
            'Y-axis label',...
            'Adds legend at the secified location',...
            'Fill markers if true',...
            'Specify custom color palette as a dictionary of group members to a 3-member color code'
            };
config = struct('name', pnames,...
                'default', dflts,...
                'help', arg_help);
opt = struct('prog', mfilename, 'desc', 'Scatter plot with grouping variable');
args = mortar.common.ArgParse.getArgs(config, opt, varargin{:});

% Don't plot anything if either x or y is empty
if isempty(x) || isempty(y),
   if nargout>0
       h = [];
   end
   return
end

if (ndims(x)==2) && any(size(x)==1), x = x(:); end
if (ndims(y)==2) && any(size(y)==1), y = y(:); end

if ndims(x)>2 || ndims(y)>2
   error('stats:gpscatter:MatrixRequired',...
         'X and Y must be 2-D.');
end
if size(x,1)~=size(y,1)
   error('stats:gpscatter:InputSizeMismatch',...
         'X and Y must have the same length.');
end

if (nargin > 2) && ~isempty(g)
   [g,gn,ignore1,ignore2,maxgrp] = mgrp2idx(g,size(x,1),','); %#ok<ASGLU>
   ng = max(g);
else
   g = [];
   gn = [];
   ng = 1;
   maxgrp = 1;
end

if (~isempty(g)) && (length(g) ~= size(x,1)),
   error('stats:gpscatter:InputSizeMismatch',...
         'There must be one value of G for each row of X.');
end

if (isempty(args.size))
   args.size = repmat(get(0, 'defaultlinemarkersize'), size(args.symbol));
   if any(args.symbol=='.'),
      units = get(gcf,'units');
      set(gcf,'units','pixels');
      pos = get(gcf,'Position');
      set(gcf,'units',units);
      args.size(args.symbol=='.') = max(1,min(15, round(15*min(pos(3:4))/size(x,1))));
   end
end

newplot;
if isempty(args.color)
%     args.color = hsv(maxgrp);
    args.color = get_palette(maxgrp);
end
hh = iscatter(x, y, g, args.color, args.symbol, args.size);

% palette specified, override defaults.
if isequal(class(args.palette), 'containers.Map')
    for ii=1:length(gn)
        p = args.palette(gn{ii});
        if isa(p, 'double') && isequal(numel(p), 3)
            set(hh(ii), 'color', p);
        else
            set(hh(ii), 'color', get_color(args.palette(gn{ii})))
        end
    end
end

if args.fill_marker
    for ii=1:length(hh)
        set(hh(ii), 'markerfacecolor', get(hh(ii), 'color'));
    end
end

% Label plots
if (~isempty(args.xlabel)), xlabel(texify(deblank(args.xlabel))); end
if (~isempty(args.ylabel)), ylabel(texify(deblank(args.ylabel))); end

% Add behavior object to lines, to customize datatip text
dataCursorBehaviorObj = hgbehaviorfactory('DataCursor');
set(dataCursorBehaviorObj,'UpdateFcn',{@gpscatterDatatipCallback,args.xlabel,args.ylabel});
for i=1:ng
    hgaddbehavior(hh(i),dataCursorBehaviorObj);
    setappdata(hh(i),'group',i);
    if ~isempty(gn)
        setappdata(hh(i),'groupname',gn{i});
    end
    if ~isempty(g)
        gind = find(g==i);
        setappdata(hh(i),'gind',gind);
    end
end

% Create legend if requested
if (args.show_legend && ~isempty(gn))
   t = find(ismember(1:size(gn,1),g));
   lh = legend(hh(t), texify(gn(t,:)), 'location', args.legend_location);
else
    lh = [];
end

% Nudge X axis limits if points are too close
xlim = get(gca, 'XLim');
d = diff(xlim);
xlim(1) = min(xlim(1), min(min(x))-0.05*d);
xlim(2) = max(xlim(2), max(max(x))+0.05*d);
set(gca, 'XLim', xlim);

if (nargout>0), h = hh; end

% Store information for gname
set(gca, 'UserData', {'gpscatter' x y g});

% -----------------------------
function datatipTxt = gpscatterDatatipCallback(obj,evt,xnam,ynam)

target = get(evt,'Target');
ind = get(evt,'DataIndex');
pos = get(evt,'Position');

group = getappdata(target,'group');
groupname = getappdata(target,'groupname');
gind = getappdata(target,'gind');

if isempty(xnam)
    xnam = 'x';
end
if isempty(ynam)
    ynam = 'y';
end
if isempty (gind)
    % One group
    % Leave group name alone, it may be empty
    % Line index number is the same as the original row
    obsind = ind;
else
    % Multiple groups
    % If group name not given, assign it its number
    if isempty(groupname)
        groupname = num2str(group);
    end
    % Map line index to the original row
    obsind = gind(ind);
end

datatipTxt = {...
    [xnam ': ' num2str(pos(1))]...
    [ynam ': ' num2str(pos(2))]...
    ''...
    ['Observation: ' num2str(obsind)]
    };

if ~isempty(groupname)
    datatipTxt{end+1} = ['Group: ' groupname];
end

