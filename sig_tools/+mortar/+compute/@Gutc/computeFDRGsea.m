function [qval, numer, denom] = computeFDRGsea(ncs_obs,...
    cdf_null_pos, cdf_null_neg, apply_qval_adjust, apply_smooth, log_xform)
% computeFDRGsea Compute FDR q-values as in GSEA (See Reference below).
% [qval, num, denom] = computeFDRGSEA(ncs_obs, cdf_null_pos, cdf_null_neg, apply_null_adjust, apply_smooth, log_xform)
% ncs_obs : 1d vector of observed normalized connectivity scores across all
% genesets tested
% cdf_null_pos : empirical CDF of positive NCS null scores for all genesets, 
% 2-column matrix [X, F] returned by the getSignedCDF function
% cdf_null_neg : empirical CDF of negative null scores for for all genesets,
% 2-column matrix returned by the getSignedCDF function
% apply_qval_adjust : linearly interpolate q-values for scores exceeding
% the null distribution
% apply_smooth : ensure q-values are monotonic with respect to the NCS scores
% log_xform: return surprisals i.e -log10(q_value) if true
%
% Reference: See the Multiple Hypothesis Testing section in
% Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based
% approach for interpreting genome-wide expression profiles. Proc. Natl.
% Acad. Sci. U. S. A. 102, 15545?15550 (2005)

assert(mortar.util.Array.is1d(ncs_obs), 'NCS_OBS must be a 1d array');
assert(isscalar(apply_qval_adjust) && islogical(apply_qval_adjust), 'APPLY_QVAL_ADJUST must be a boolean scalar');
assert(isscalar(apply_smooth) && islogical(apply_smooth), 'APPLY_SMOOTH must be a boolean, scalar');
assert(isscalar(log_xform) && islogical(log_xform), 'LOG_XFORM must be a boolean, scalar');

sz = size(ncs_obs);
qval = nan(sz);
numer = nan(sz);
denom = nan(sz);

is_pos_obs = ncs_obs > 0;
is_neg_obs = ncs_obs < 0;

% positive scores
% observed and null quantiles for each score
if nnz(is_pos_obs) && nnz(is_pos_obs)>2
    obs_pos = ncs_obs(is_pos_obs);
    cdf_obs_pos = mortar.compute.Gutc.getSignedCDF(obs_pos, 1);
    q_obs_pos = mortar.compute.Gutc.lookupQuantile(cdf_obs_pos, obs_pos);
    q_null_pos = mortar.compute.Gutc.lookupQuantile(cdf_null_pos, obs_pos);
    % proportion of positive null scores >= NCS*
    numer(is_pos_obs) = 1 - q_null_pos;
    % prop. of positive obs scores whose NCS(S) >= NCS*
    denom(is_pos_obs) = (1 - q_obs_pos) + eps;
    qraw = numer(is_pos_obs) ./ denom(is_pos_obs);
    qval(is_pos_obs) = clip(qraw, eps, 1);
    % Linear q-value adjustment for scores exceeding the null distribution.
    if apply_qval_adjust
        qval(is_pos_obs) = mortar.compute.Gutc.adjustQvalue(...
            qraw, obs_pos, numer(is_pos_obs),...
            denom(is_pos_obs), 1);
    end    
    % Smooth q-values to be monotonic wrt to ncs
    if apply_smooth
        qval(is_pos_obs) = mortar.compute.Gutc.smoothFDR(qval(is_pos_obs), obs_pos, 1);
    end    
end

% negative scores, same steps as for positive but with reversed sign
if nnz(is_neg_obs) && nnz(is_neg_obs)>2
    obs_neg = ncs_obs(is_neg_obs);
    cdf_obs_neg = mortar.compute.Gutc.getSignedCDF(obs_neg, -1);
    q_obs_neg = mortar.compute.Gutc.lookupQuantile(cdf_obs_neg, obs_neg);
    q_null_neg = mortar.compute.Gutc.lookupQuantile(cdf_null_neg, obs_neg);
    numer(is_neg_obs) = q_null_neg;
    denom(is_neg_obs) = (q_obs_neg) + eps;
    qraw = numer(is_neg_obs) ./ denom(is_neg_obs);
    qval(is_neg_obs) = clip(qraw, eps, 1);
    % Linear q-value adjustment for scores exceeding the null distribution.
    if apply_qval_adjust
        qval(is_neg_obs) = mortar.compute.Gutc.adjustQvalue(...
            qraw, obs_neg, numer(is_neg_obs),...
            denom(is_neg_obs), -1);
    end    
    % Smooth q-values to be monotonic wrt to ncs
    if apply_smooth
        qval(is_neg_obs) = mortar.compute.Gutc.smoothFDR(qval(is_neg_obs), obs_neg, -1);
    end    
end

if (log_xform)
    numer = -log10(numer + eps);
    denom = -log10(denom + eps);
    qval = -log10(qval + eps);
end

end




