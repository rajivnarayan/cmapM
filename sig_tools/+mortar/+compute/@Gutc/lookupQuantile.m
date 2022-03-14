% Lookup quantiles from eCDF
function q = lookupQuantile(cdf_val, x)
interp_method = 'pchip';
v = cdf_val(:, 1);
f = cdf_val(:, 2);
min_v = min(v);
max_v = max(v);
q = interp1(v, f,...
    clip(x, min_v, max_v),...
    interp_method);
end