function runAnalysis_(obj, varargin)
args = obj.getArgs;
obj.res_ = main(args);
end

function res = main(args)
% Main function

%% Run Cmap Query
dbg(args.verbose, '# Running Set enrichment');
score_ds = parse_gctx(args.score);
rank_ds = parse_gctx(args.rank);

has_nan = any(isnan(score_ds.mat));
if has_nan
    dbg(args.verbose, '# NaN scores found! runnning the GSEA after filtering them...');
    res = run_gsea_nan(score_ds, rank_ds, args);
else
    dbg(args.verbose, '# No NaN scores found, runnning the matrix-optimized algorithm...');
    res = run_gsea_nonan(score_ds, rank_ds, args);
end

end

% GSEA pre-ranked for score matrix w/o NaNs
function res = run_gsea_nonan(score_ds, rank_ds, args)
query_result = mortar.compute.Connectivity.runCmapQuery(...
    'score', score_ds, ...
    'rank', rank_ds, ...
    'up', args.up,...
    'es_tail', args.es_tail,...
    'metric', args.metric,...
    'sig_meta', args.sig_meta,...
    'query_meta', args.query_meta,...
    'max_col', args.max_col);

% transpose the cs file to [sets x num_cols_in_score]
query_result.cs = transpose_gct(query_result.cs);
set_sizes = [query_result.uptag.len]';

dbg(args.verbose, '# Normalizing scores and computing null distributions');
[nes_result, fdr_result] = mortar.compute.Gutc.normalizeQueryWithPermutedNull(...
    query_result.cs, score_ds, rank_ds, set_sizes,...
    'num_perm', args.num_perm);

res = struct('args', args,...
    'query_result', query_result,...
    'nes_result', nes_result,...
    'fdr_result', fdr_result);
end

% GSEA pre-ranked for score matrix with NaNs
function res = run_gsea_nan(score_ds, rank_ds, args)

[nr, nc] = size(score_ds.mat);
sets_gmt = parse_geneset(args.up);
nset = length(sets_gmt);
%src_set_size = [sets_gmt.len]';
%eff_size_all = nan(nset, nc);

for ii=1:nc
    this_score = ds_slice(score_ds, 'cidx', ii);
    this_score = ds_delete_missing(this_score);
    this_rank = ds_slice(rank_ds, 'cid', this_score.cid,...
                         'rid', this_score.rid);
    this_query_res = mortar.compute.Connectivity.runCmapQuery(...
                'up', sets_gmt,...
                'score', this_score,...
                'rank', this_rank,...
                'es_tail', args.es_tail,...
                'metric', args.metric,...
                'sig_meta', args.sig_meta,...
                'query_meta', args.query_meta,...
                'max_col', args.max_col);
    this_query_res.cs = transpose_gct(this_query_res.cs);
    set_sizes = [this_query_res.uptag.len]';
    [nes_ds, fdr_ds] = mortar.compute.Gutc.normalizeQueryWithPermutedNull(...
        this_query_res.cs, this_score, this_rank, set_sizes,...
        'num_perm', args.num_perm);
    
    if isequal(ii, 1)
        query_result = this_query_res;        
        nes_result = nes_ds;
        fdr_result = fdr_ds;
    else
        query_result = merge_query_result(query_result, this_query_res);
        nes_result = merge_two(nes_result, nes_ds);
        fdr_result = merge_two(fdr_result, fdr_ds);
    end
end

% save the input query
query_result.uptag = sets_gmt;
query_result.dntag = [];

res = struct('args', args,...
    'query_result', query_result,...
    'nes_result', nes_result,...
    'fdr_result', fdr_result);
end

function res1 = merge_query_result(res1, res2)

ds_fields = {'cs', 'cs_up', 'cs_dn', 'leadf_up', 'leadf_dn'};
nds = length(ds_fields);
for ii=1:nds
    res1.(ds_fields{ii}) = merge_two(res1.(ds_fields{ii}), res2.(ds_fields{ii}));
end
end
