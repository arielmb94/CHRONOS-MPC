function mpc = update_mpc_f0_hess(mpc)

mpc.hessCost = zeros(mpc.Nu+mpc.Nx);

if ~isempty(mpc.Qe)
    mpc.hessCost = mpc.hessCost + mpc.hessErrTerm;
end

if ~isempty(mpc.Ru)
    mpc.hessCost = mpc.hessCost + mpc.hessCtrlTerm;
end

if ~isempty(mpc.Rdu)
    mpc.hessCost = mpc.hessCost + mpc.hessDiffCtrlTerm;
end

if ~isempty(mpc.P)
    mpc.hessCost = mpc.hessCost + mpc.hessTerminalCost;
end

if ~isempty(mpc.Qz)
    mpc.hessCost = mpc.hessCost + mpc.hessPerfTerm;
end

end