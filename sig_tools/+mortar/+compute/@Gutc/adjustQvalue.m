function qval = adjustQvalue(qval, ncs, numer, denom, use_right_tail)

[y0, iminq] = min(qval);

if (use_right_tail>0)
    % indexes of positive scores with min qvalues    
    edge_idx = find(ncs > ncs(iminq));
    if any(edge_idx)
        ncs_edge = ncs(edge_idx);
        [x0, x0_idx] = min(ncs_edge);
        [x1, x1_idx] = max(ncs_edge);
        numer_at_x1 = numer(edge_idx(x1_idx));
        denom_at_x1 = denom(edge_idx(x1_idx));
        to_fix = ncs >= x0;
    end
else
    % indexes of negative scores with min qvalues
    edge_idx = find(ncs < ncs(imin(qval)));
    if any(edge_idx)
        ncs_edge = ncs(edge_idx);
        [x0, x0_idx] = max(ncs_edge);
        [x1, x1_idx] = min(ncs_edge);        
        numer_at_x1 = numer(edge_idx(x1_idx));
        denom_at_x1 = denom(edge_idx(x1_idx));
        to_fix = ncs <= x0;
    end
end

if any(edge_idx) && (denom_at_x1 < numer_at_x1)
    y1 = min(y0, max(numer_at_x1, denom_at_x1));
    m = (y1 - y0) / (eps+(x1 - x0));
    c = y1 - m*x1;
    qval_fixed = clip(m*ncs(to_fix)+c, eps, 1);
    qval(to_fix) = qval_fixed;
end

end