% DEPRECATED, Linearly scale q-values where scores exceed the null distribution.
function qval = adjustQvalue0(qval, xq, use_right_tail)
% In this region all the raw q-values are zero. 

if (use_right_tail>0)
    % indexes of positive scores with min qvalues
    %edge_idx = xq > min(xq(qval <= min(qval)));
    edge_idx = xq > xq(imin(qval));
    if any(edge_idx)
        [x0, min_idx] = min(xq(edge_idx));
        q_tmp = qval(edge_idx);
        y0 = q_tmp(min_idx);
        x1 = max(xq);
        y1 = eps;
        m = (y1 - y0) / (eps+(x1 - x0));
        c = -m*x1;
        to_fix = xq >= x0;
        qval(to_fix) = (m*xq(to_fix)+c);
    end
else
    % indexes of negative scores with min qvalues
    %edge_idx = xq < max(xq(qval <= min(qval)));
    edge_idx = xq < xq(imin(qval));
    if any(edge_idx)
        [x0, max_idx] = max(xq(edge_idx));
        q_tmp = qval(edge_idx);
        y0 = q_tmp(max_idx);
        x1 = min(xq);
        y1 = eps;
        m = (y1 - y0) / (eps+(x1 - x0));
        c = -m*x1;
        to_fix = xq <= x0;
        qval(to_fix) = (m*xq(to_fix)+c);
    end
end

end