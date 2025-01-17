function [ofname, status, result] = mkheatmap(dsfile, outfile, varargin)
% MKHEATMAP Wrapper for heatmap

% Heatmap script
hm_script = fullfile(mortarconfig('cmap_bin_path'), 'heatmap');
pnames = {'format', 'column_text', 'row_text', ...
    'iszscore', 'debug', 'cluster_row',...
    'cluster_col', 'cluster_distance', 'cluster_linkage',...
    'column_color', 'row_color', 'title',...
    'color_scheme'};
dflts = {'png', {'id'}, {'id'}, ...
    false, true, false, ...
    false, 'correlation', 'complete',...
    '', '', '',...
    ''};

distance_dict = containers.Map({'city block', 'euclidean',...
                                'kendall', 'correlation',...
                                'spearman'},...
                        {'City-block distance', 'Euclidean distance',...
                        'Kendall''s tau', 'One minus pearson correlation', 'One minus spearman rank correlation',...
                            });
linkage_dict = containers.Map({'average', 'complete', 'single'},...
    {'Average Linkage', 'Complete Linkage', 'Single Linkage'});

color_scheme_dict = containers.Map({'rankpoint_30',...
                                    'rankpoint_50',...
                                    'rankpoint_70',...
                                    'rankpoint_80',...
                                    'rankpoint_90',...
                                    'rankpoint_95',....
                                    'zscore_pm2',...
                                    'binary',...
                                    'lfcvc'},...
                                    {'\-100.0:#0000FF,\-30.0:#FFFFFF,30.0:#FFFFFF,100.0:#FF0000',...
                                    '\-100.0:#0000FF,\-50.0:#FFFFFF,50.0:#FFFFFF,100.0:#FF0000',...
                                    '\-100.0:#0000FF,\-70.0:#FFFFFF,70.0:#FFFFFF,100.0:#FF0000',...
                                    '\-100.0:#0000FF,\-80.0:#FFFFFF,80.0:#FFFFFF,100.0:#FF0000',...
                                    '\-100:#0000FF,\-90:#FFFFFF,90:#FFFFFF,100:#FF0000',...                                    
                                    '\-100:#0000FF,\-95:#FFFFFF,95:#FFFFFF,100:#FF0000',...                                    
                                    '\-10:#0000FF,\-2:#FFFFFF,2:#FFFFFF,10:#FF0000',...
                                    '0:#FFFFFF,0.5:#1F77B4',...
                                    '\-3:#008837,\-1:#a6dba0,\-0.7:#f7f7f7,0.7:#f7f7f7,1.0:#c2a5cf,3:#7b3294'
                                    });
                                
args = parse_args(pnames, dflts, varargin{:});
[p, f] = fileparts(outfile);
ofname = fullfile(p, [f, '.', args.format]);

% Heatmap options
optflag = sprintf('--format %s', args.format);
if args.iszscore
    optflag = sprintf('%s --zscores', optflag);
end

% Clustergram options
distance_id = distance_dict(args.cluster_distance);
linkage_id = linkage_dict(args.cluster_linkage);
if args.cluster_col
    optflag = sprintf('%s --cluster-columns --column-dist ''%s'' --linkage ''%s''', optflag, distance_id, linkage_id);
end
if args.cluster_row
    optflag = sprintf('%s --cluster-rows --row-dist ''%s'' --linkage ''%s''', optflag, distance_id, linkage_id);
end

% title
if ~isempty(args.title)
   optflag = sprintf('--title %s', args.title); 
end
% column and row headers
if ~isempty(args.column_color)
    ccolor = print_dlm_line(args.column_color, 'dlm', ' ');
    optflag = sprintf('--column-color %s %s', ccolor, optflag);
end
if ~isempty(args.row_color)
    rcolor = print_dlm_line(args.row_color, 'dlm', ' ');
    optflag = sprintf('--row-color %s %s', rcolor, optflag);
end

if ~isempty(args.color_scheme)
    if color_scheme_dict.isKey(args.color_scheme)
        color_scheme = color_scheme_dict(args.color_scheme);
    else
        color_scheme = args.color_scheme;
    end
    optflag = sprintf('--color-scheme ''%s'' %s', color_scheme, optflag);
end

ctext = print_dlm_line(args.column_text, 'dlm', ' ');
rtext = print_dlm_line(args.row_text, 'dlm', ' ');
flags = sprintf('--column-text %s --row-text %s %s', ctext, rtext, optflag);

% run it
runstring = sprintf('%s --gct %s -o %s %s', hm_script, dsfile, ofname, flags);
dbg (args.debug, 'Saving heatmap to: %s...', ofname)
dbg (args.debug, '%s', runstring)
[status, result] = system(runstring, '-echo');
dbg (args.debug, 'done.')

end