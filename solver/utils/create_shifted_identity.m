%%%%%%%%%%%
% Generates a shifted indentity matrix from a boolean vector u based on
% where the booleans are placed on u
%
% Example:
% u = [1 0 1] --> I = [1 0 0; 0 0 1]
% u = [1 0 1 0 0 1] I = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 1]
%%%%%%%%%%
function I_shift = create_shifted_identity(u)
    nu = sum(u);
    Nu = length(u);

    I_shift = zeros(nu,Nu);
    
    n = 1;
    for i = 1:Nu
        if u(i)
            I_shift(n,i) = 1;
            n = n+1;
        end    
    end
end