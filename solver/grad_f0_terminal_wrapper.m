function mpc = grad_f0_terminal_wrapper(mpc)
mpc.grad_se_f0_ter=grad_f0_terminal_local(mpc.grad_se_f0_ter,mpc.s_col,mpc.P2,mpc.xN_ref,mpc.s_ter);
end
function grad_se_f0_ter = grad_f0_terminal_local(grad_se_f0_ter,s_idx,P2,xN_ref,s_ter)
grad_se_f0_ter(s_idx) = grad_se_f0_ter(s_idx) - P2*(xN_ref-s_ter);
end
