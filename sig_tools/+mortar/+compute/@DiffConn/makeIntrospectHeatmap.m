function [ofname, status, result] = makeIntrospectHeatmap(varargin)
% makeIntrospectHeatmap Generate and save an introspect heatmap
% ofname = makeIntrospect(cs, out, varargin)
% cs: Path to GCT(x) file of ordered and annotated introspect matrix (Required)
% out: Path to output image file (Required)
% --title: Plot title
% --colormap_max: float, Maximal colormap value
% --colormap_cutoff: float, Cutoff value for colormap
% --column_color: cell array, Column annotations fields to display as tracks
% --row_color: cell array, Row annotations fields to display as tracks
% --column_text: cell array, Column annotations fields to display as text
% --row_text: cell array, Row annotations fields to display as text

[args, help_flag] = getArgs(varargin{:});

if ~help_flag
    %heatmap settings
    color_scheme = sprintf(...
        '\\-%0.2f:#0000FF,\\-%0.2f:#FFFFFF,%0.2f:#FFFFFF,%0.2f:#FF0000',...
        args.colormap_max, args.colormap_cutoff,...
        args.colormap_cutoff, args.colormap_max);
    
    %% Create heatmap
    [ofname, status, result] = mkheatmap(args.cs, args.out,...
        'title', args.title, ...
        'color_scheme', color_scheme, ...
        'column_color', args.column_color, ...
        'row_color', args.row_color, ...
        'column_text', args.column_text, ...
        'row_text', args.row_text);
    
end

end

function [args, help_flag] = getArgs(varargin)
pnames = {'cs';...
    'out';...
    '--title';...
    '--colormap_max';...
    '--colormap_cutoff';...
    '--column_color';...
    '--row_color';...
    '--column_text';...
    '--row_text'};

defaults = {'';...
    '';...
    '';...
    1;...
    0;...
    '';...
    '';...
    '';...
    ''};

help_str = {'Path to GCT(x) file of ordered and annotated introspect matrix (Required)';...
    'Path to output image file (Required)';...
    'Plot title';...
    'float, Maximal colormap value';...
    'float, Cutoff value for colormap';...
    'cell array, Column annotations fields to display as tracks';...
    'cell array, Row annotations fields to display as tracks';...
    'cell array, Column annotations fields to display as text';...
    'cell array, Row annotations fields to display as text'};

config = struct('name', pnames,...
    'default', defaults,...
    'help', help_str);
opt = struct('prog', mfilename, 'desc', 'Make Introspect heatmap');

[args, help_flag] = mortar.common.ArgParse.getArgs(config, opt, varargin{:});

assert(isfileexist(args.cs), 'cs file not found : %s', args.cs);
assert(~isempty(args.out), 'Required argument out not specified');
out_path = fileparts(args.out);
mkdirnotexist(out_path);

end