function [qval_ds, pnull_ds, pobs_ds] = computeFDRGseaDs(ncs_ds, is_null)
% computeFDRGseaDs compute FDR for each column in NCS dataset
% qval_ds = computeFDRGseaDs(ncs_ds, is_null)

[nr, nc] = size(ncs_ds.mat);
qval_mat = ones(nr, nc);
pnull_mat = ones(nr, nc);
pobs_mat = ones(nr, nc);
% FDR adjustments, see computeFDRGsea for details
apply_qval_adjust = true;
apply_smooth = true;
log_xform = true;

% Get CDFs for the Null scores
ds_null = ds_slice(ncs_ds, 'ridx', is_null);
cdf_null_pos = mortar.compute.Gutc.getSignedCDF(ds_null.mat, 1);
cdf_null_neg = mortar.compute.Gutc.getSignedCDF(ds_null.mat, -1);

% Compute FDR q-values for each ranked list
for ii=1:nr
    [qval_mat(ii, :), pnull_mat(ii,:), pobs_mat(ii,:)] =...
        mortar.compute.Gutc.computeFDRGsea(ncs_ds.mat(ii, :),...
                           cdf_null_pos, cdf_null_neg,...
                           apply_qval_adjust, apply_smooth, log_xform);
end
qval_ds = mkgctstruct(qval_mat, 'rid', ncs_ds.rid, 'cid', ncs_ds.cid);
pnull_ds = mkgctstruct(pnull_mat, 'rid', ncs_ds.rid, 'cid', ncs_ds.cid);
pobs_ds = mkgctstruct(pobs_mat, 'rid', ncs_ds.rid, 'cid', ncs_ds.cid);

end
