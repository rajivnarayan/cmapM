function checkArgs_(obj)
% sanity check the parameters
args = obj.getArgs;

%% ADD INPUT VALIDATION HERE

assert(is_ds_or_file(args.score), 'Score not specified');
if isempty(args.rank)
    dbg(1, '# Rank file not specified, calculating ranks from scores. Use --rank to save time and specify custom ranks')
    args.score = parse_gctx(args.score, 'detect_numeric', false, 'has_missing_data', true);
    args.rank = score2rank(args.score);
else
    assert(is_ds_or_file(args.rank), 'Specified Rank not found');
    dbg(args.verbose, '# Using pre-computed ranks')
end

%assert(is_struct_or_file(args.sig_meta), 'Signature metadata not specified');
assert(is_gset_or_file(args.up), 'Up geneset expected');

args.up = parse_geneset(args.up);
args.up = setfilter(args.up, args.score.rid, args.min_set_size, args.max_set_size);

assert(~isempty(args.up), 'Up geneset after filtering is empty');

% Update args
obj.setArgs(args);

end

function tf = is_struct_or_file(s)
tf = isstruct(s) || isfileexist(s, 'file');
end

function tf = is_ds_or_file(s)
tf = isds(s) || isfileexist(s, 'file');
end

function tf = is_gset_or_file(s)
tf = isgeneset(s) || ~isempty(uri_type(s)) || isfileexist(s, 'file');
end
