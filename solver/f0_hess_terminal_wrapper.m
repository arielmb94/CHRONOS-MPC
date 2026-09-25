function mpc = f0_hess_terminal_wrapper(mpc)
mpc.Q_f0_ter = f0_hess_terminal_local(mpc.Q_f0_ter,mpc.P2,mpc.s_col);
end

function Q_f0_ter = f0_hess_terminal_local(Q_f0_ter,P2,s_col)
Q_f0_ter(s_col,s_col) = Q_f0_ter(s_col,s_col) + P2;
end
