function saveDataset(res, out_path)

req_input = struct('ds', 'gct',...
        'ds_high', 'gct',...
        'ds_low', 'gct',...
        'support', 'gct',...
        'support_high', 'gct',...
        'support_low', 'gct',...
        'support_pct', 'gct',...
        'support_pct_high', 'gct',...
        'support_pct_low', 'gct',...
        'row_meta', 'tbl',...
        'col_meta', 'tbl');

ds_type = struct2cell(req_input);
req_field = fieldnames(req_input);
assert(all(isfield(res, req_field)), 'Required fields not found in res')

nds = length(ds_type);
for ii=1:nds
    switch(ds_type{ii})
        case 'gct'
            mkgctx(fullfile(out_path, sprintf('%s.gct', req_field{ii})),...
                   res.(req_field{ii}));
        case 'tbl'
            mktbl(fullfile(out_path, sprintf('%s.txt', req_field{ii})),...
                   res.(req_field{ii}));
        otherwise
            error ('Unsupported datatype :%s', ds_type{ii});
    end

end

end