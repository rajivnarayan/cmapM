% compute eCDF for specified tail signed scores
function d = getSignedCDF(ncs, tail)

is_nz = sign(tail)*ncs(:) > 0;
if (nnz(is_nz)>1)
    [f, v] = cdfcalc(ncs(is_nz));
    d = [v, f(1:end-1)];
else
    d = nan(1, 2);
end

end