function saveResult(res, outpath, varargin)
% saveResult Save results produced by runCmapQuery

config = struct('name', {'--save_minimal', '--appenddim', '--use_gctx'},...
    'default', {false, true, true},...
    'help', {'Save minimal results if true', 'Append dimension to name', 'Output as GCTX'});
opt = struct('prog', mfilename, 'desc', 'Save results produced by runCmapQuery');
[args, help_flag] = mortar.common.ArgParse.getArgs(config, opt, varargin{:});

mkdirnotexist(outpath);
gct_writer = ifelse(args.use_gctx, @mkgctx, @mkgct);

% save arguments
% print_args('query_tool', fullfile(outpath, 'query_tool_params.txt'), res.args);

% save genesets
mkgmt(fullfile(outpath, 'up.gmt'), res.uptag);
mkgmt(fullfile(outpath, 'dn.gmt'), res.dntag);

% save result matrices
result_fields = {'cs', 'cs_up', 'cs_dn', 'leadf_up', 'leadf_dn'};
to_save = true(size(result_fields));
if args.save_minimal
    dbg(1, 'save_minimal specified, saving just cs')
    to_save = strcmp('cs', result_fields);   
end

for ii=1:length(result_fields)    
    if isfield(res, result_fields{ii}) && isds(res.(result_fields{ii})) && to_save(ii)
        outfile = fullfile(outpath, sprintf('%s.gctx', result_fields{ii}));
        gct_writer(outfile, res.(result_fields{ii}), 'appenddim', args.appenddim);
    end
end

end
