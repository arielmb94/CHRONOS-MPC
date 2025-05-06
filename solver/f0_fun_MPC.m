function J = f0_fun_MPC(Qe,err,N,ny,R,deltaU,nu,P,s_ter)

    J = 0;

    if ~isempty(Qe)

        Je = cost_Error(Qe,err,N,ny);

        J = J + 0.5*Je;

    end    

    if ~isempty(R)

        Ju = cost_U(R,deltaU,N,nu);

        J = J + 0.5*Ju;

    end  

    if ~isempty(P)

        Jter = cost_Terminal(P,s_ter);

        J = J + Jter;

    end  

        
end


function Je = cost_Error(Qe,err,N,ny)

Je = 0;
for k = 1:N

    errk = err((k-1)*ny+1:k*ny);

    Je = Je + errk'*Qe*errk;

end

end

function Ju = cost_U(R,deltaU,N,nu)

Ju = 0;
for k = 1:N

    delta_uk = deltaU((k-1)*nu+1:k*nu);

    Ju = Ju + delta_uk'*R*delta_uk;

end

end

function Jter = cost_Terminal(P,s_ter)

Jter = s_ter'*P*s_ter;

end

