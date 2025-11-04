T = readtable('Traj.csv');

X = T.Var1;
Y = T.Var2;

%% Compute heading angle

%for smoothing
i_ah = 15; %samples ahead = meters/dX
nl = length(X);

for i = 1:length(X)  
    if i-i_ah<1
        deltaXi = X(i+i_ah)-X(i-i_ah+nl);
        deltaYi = Y(i+i_ah)-Y(i-i_ah+nl);
        
    elseif i+i_ah>nl
        
        deltaXi = X(i+i_ah-nl)-X(i-i_ah);
        deltaYi = Y(i+i_ah-nl)-Y(i-i_ah);
        
    else
        
        deltaXi = X(i+i_ah)-X(i-i_ah);
        deltaYi = Y(i+i_ah)-Y(i-i_ah);

    end    

    Theta(i) = atan2(deltaYi,deltaXi);
end

%%
Track.X = X;
Track.Y = Y;
Track.Theta = Theta;