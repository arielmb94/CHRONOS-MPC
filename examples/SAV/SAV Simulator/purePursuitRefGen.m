function Ref = purePursuitRefGen(Xgl,Ygl,YawHead,Vx,Vy,tp,Track)
   
    Xpath = Track.X;
    Ypath = Track.Y;
    Yawpath = Track.Theta;

    if abs(Vx)<0.6
        Vx = 0.6;
    end
    Beta = atan(Vy/Vx); %Car body slip angle
    
    %O = (Vx-0.5)/(2-0.5);

    L = tp*Vx; %Lookahead distance
    %L = 0.5;
    %L = tp;
    global pathindex
    
    
    % Find new Path Index
    pathindex = findMinDist(Xgl,Ygl,Xpath,Ypath,pathindex);
    
     % Find Pursuit Point
     XL = Xgl+L*cos(YawHead+Beta); %X position at lookahead distance if actual velocity is maintened
     YL = Ygl+L*sin(YawHead+Beta); %Y position at lookahead distance if actual velocity is maintened
    %pursuitindex = findMinDist(XL,YL,Xpath,Ypath,pathindex);
    pursuitindex = findPointL(L,Xgl,Ygl,Xpath,Ypath,pathindex);
    Xpursuit = Xpath(pursuitindex);
    Ypursuit = Ypath(pursuitindex);
    Yawpursuit = Yawpath(pursuitindex);

    % Find yawrate       
    alpha = atan2(Ypursuit-Ygl,Xpursuit-Xgl)-YawHead; 
    k = 2*sin(alpha)/norm([Xgl,Ygl]-[Xpursuit,Ypursuit]);  %k = 2*sin(alpha)/L
    Yawrate_Ref1 = k*Vx;

    tau = 0.5;
    Yawrate_Ref2 = getAngleDifference(YawHead,Yawpursuit)/tau;

    mix = 0.9;
    Yawrate_Ref = mix*Yawrate_Ref1 + (1-mix)*Yawrate_Ref2;
    
    % Find ye
    % Find Pursuit Point
    ye = sin(alpha)*L;
    ye = sign(sin(alpha))*norm([XL,YL]-[Xpursuit,Ypursuit]);
    
    %Output References
    %k = Kpath(pursuitindex);
    Ref = [ye,Yawrate_Ref,k];
    %Ref = Yawrate_Ref;
end


function y = findMinDist(X,Y,Xpath,Ypath,index)
    DO = true;
    while DO
        di = norm([X,Y]-[Xpath(index),Ypath(index)]);
        index = index+1;
        if index > length(Xpath)
            index = 1;
        end
        dnew = norm([X,Y]-[Xpath(index),Ypath(index)]);
        DO = di>dnew;
    end
    y = index;
end

function y = findPointL(L,X,Y,Xpath,Ypath,index)
    DO = true;
    while DO
        di = norm([X,Y]-[Xpath(index),Ypath(index)]);
        index = index+1;
        if index > length(Xpath)
            index = 1;
        end
        DO = L>di;
    end
    y = index;
end

function r = getAngleDifference(b1, b2)
    r = mod((b2 - b1),2*pi);

    if r >= pi    %#180.0
    
        r = r - 2*pi;   %#360.0
    
    end
end
