function cnstr = set_active_set_min(cnstr,x_min,x_min_activ,n)

if ~isempty(x_min_activ)

    % check that activation range is within variable range
    index = x_min_activ < x_min;
    x_min_activ(index) = x_min(index);

    cnstr.min_activ_set = 1;
    cnstr.min_activ_lim = x_min_activ;
    cnstr.min_activ_indicator = zeros(n,1);
else

    cnstr = empty_active_set_min(cnstr);
end

end