function y = GridK_PL_SFV2(x,K)
    
    k_g = 0.2;
    l_t= 0.282;
    l_m = 0.246;
    l_b = 0.29;
    l_cb = 0.276;
    r_ts = 0.1;
    r_ms = 0.155;
    m_cb = 0.068;
    m_b = 0.022;
    m_m = 0.014;
    m_mr = 0.236;
    m_ms = 0.219;
    L_ah = 0.86e-3;
    L_av = 0.86e-3;
    k_a = 0.0202;
    B_tr = 0.0086;
    B_mr = 0.0026;
    k_tvp = 0.0168;
    k_tvn = 0.0155;
    k_fhp = 0.0566;
    k_fhn = 0.0660;
    k_fvp = 0.3819;
    k_fvn = 0.2197;
    k_chp = 0.0158;
    k_chn = 0.0111;
    k_cvp = 0.0623;
    m_t = 0.015;
    m_tr = 0.221;
    m_ts = 0.119;
    J_tr = 0.0059;
    J_mr = 0.0254;
    R_a = 8;
    k_t = 0;
    k_m = 0.0017;
    k_1 = 6.5;
    k_2 = 8.5;
    k_oh = 0.0185; 
    k_ov = 0.1026;
    k_cvn = 0.0563;
    g = 9.81;
    Thtv0 = -0.619426178368110;
    
    K_A = (0.5*m_t+m_tr+m_ts)*l_t;
    K_B = (0.5*m_m+m_mr+m_ms)*l_m;
    K_C = 0.5*m_b*l_b+m_cb*l_cb;
    K_D = (m_m/3+m_mr+m_ms)*(l_m^2)+(m_t/3+m_tr+m_ts)*(l_t^2);
    K_E = m_b*(l_b^2)/3+m_cb*(l_cb^2);
    K_F = m_ms*r_ms^2+0.5*m_ts*r_ts^2;
    K_H = K_A*l_t+K_B*l_m+0.5*m_b*(l_b^2)+m_cb*l_cb;

    
    wh = x(1);
    omh = x(2);
    thth = x(3);
    wv = x(4);
    thtv = x(6);
    
    TththRef = x(13);
    TthtvRef = x(14);
    
    Keval = K{1}+K{2}*wh+K{3}*wv+K{4}*thtv+K{5}*cos(thtv)+K{6}*sin(thtv);
   
    
    %Keval = K{1}+wh*K{2}+K{3}*omh+K{4}*thth+K{5}*wv+K{6}*thtv+K{7}*cos(thtv)+K{8}*sin(thtv);
    
    x(6) = x(6)-Thtv0;
        
    OmhRef = (TththRef-thth)/3;
    OmvRef = (TthtvRef-thtv)/5;
    
    f5r = (g*(K_C*sin(TthtvRef)+(K_B-K_A)*cos(TthtvRef))+k_ov*OmvRef+(OmvRef^2)*K_H*sin(TthtvRef)*cos(TthtvRef))/(l_m+k_g*OmhRef*cos(TthtvRef));
    if f5r >= 0
        WvRef = sqrt(f5r/k_fvp);
    else
        WvRef = -sqrt(-f5r/k_fvp);
    end
    f2r = (f3(TththRef)-f6(TthtvRef)+k_oh*OmhRef-...
            ((k_m*WvRef*sin(TthtvRef)*OmvRef*(K_D*(cos(TthtvRef)^2)-K_E*(sin(TthtvRef)^2)-K_F-2*K_E*(cos(TthtvRef)^2)))/(K_D*(cos(TthtvRef)^2)+K_E*(sin(TthtvRef)^2)+K_F)))/(l_t*cos(TthtvRef));
    if f2r >= 0
        WhRef = sqrt(f2r/k_fhp);
    else
        WhRef = -sqrt(-f2r/k_fhp);
    end
    
    wr = [WhRef;WvRef];
    
    urh = ((B_tr*R_a+(k_a^2))*WhRef+f1(WhRef)*R_a)/(k_a*k_1);
    urv = ((B_mr*R_a+(k_a^2))*WvRef+f4(WvRef)*R_a)/(k_a*k_2);
    
    
    x = x(1:12);
    %u = [urh;urv]+Keval*x;
    u = Keval*x;
    y = [u;wr];
end

function y = f1(Wh)

    k_thp = 0.0027;
    k_thn = 0.0028;

    if Wh>=0
        y = k_thp*Wh^2;
    else
        y = -k_thn*Wh^2;
    end
end

function y = f4(Wv)

    k_tvp = 0.0168;
    k_tvn = 0.0155;

    if Wv>=0
        y = k_tvp*Wv^2;
    else
        y = -k_tvn*Wv^2;
    end
end

function y = f3(Thth)

    k_chp = 0.0158;
    k_chn = 0.0111;

    if Thth>=0
        y = k_chp*Thth;
    else
        y = k_chn*Thth;
    end
end

function y = f6(Thtv)

    k_cvp = 0.0623;
    k_cvn = 0.0563;
    
    Thtv0 = -0.619426178368110;

    if Thtv>=Thtv0
        y = k_cvp*(Thtv-Thtv0)^2;
    else
        y = k_cvn*(Thtv-Thtv0)^2;
    end
end

