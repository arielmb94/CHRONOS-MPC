function cnstr = set_active_set_max(cnstr,x_max,x_max_activ,n)

if ~isempty(x_max_activ)

    % check that activation range is within variable range
    index = x_max_activ > x_max;
    x_max_activ(index) = x_max(index);

    cnstr.max_activ_set = 1;
    cnstr.max_activ_lim = x_max_activ;
    cnstr.max_activ_indicator = zeros(n,1);
else

    cnstr = empty_active_set_max(cnstr);
end

end