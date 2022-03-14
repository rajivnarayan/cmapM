% Smooth q-values to be monotonic wrt to ncs
function qval = smoothFDR(qval, ncs, use_right_tail)

if (use_right_tail>0)
    sort_order = 'ascend';
else
    sort_order = 'descend';
end
[~, srt_idx] = sort(ncs, sort_order);
nval = numel(qval);
qval_sort = qval(srt_idx);
last_q = qval_sort(1);
for ii=2:nval
    if (qval_sort(ii) < last_q)
        last_q = qval_sort(ii);
    end
    if (last_q < qval_sort(ii))
        qval(srt_idx(ii)) = last_q;
    end
end

end