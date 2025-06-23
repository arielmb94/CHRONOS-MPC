function v_full = fill_vec(v,n,N,fill_style)

l_v = length(v);
step_v = l_v/n;      % number of steps filled already
step_full = N/n;     % number of steps to fill

miss_step = step_full-step_v; % steps needed to be filled

v_full = zeros(N,1);
v_full(1:l_v) = v;

if fill_style == 0 % fill with zeros
    v_full(l_v+1:end) = zeros(miss_step*n,1);
else                % fill repeat last

    last_v = v(l_v-n+1:l_v); % get last step vector

    for i = 1:miss_step % fill with repeated last vector
        v_full(l_v+(i-1)*n+1:l_v+i*n) = last_v;
    end
end

end

