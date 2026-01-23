function out = Obs(I,w_l,w_r,V_l,V_r,A,Bu,C,L,x_k)

    y = [w_l,w_r,I]';
    u = [V_l,V_r]';
    x_k1 = A*x_k+Bu*u+L*(y-C*x_k);
    
    out = [x_k;x_k1];

end

