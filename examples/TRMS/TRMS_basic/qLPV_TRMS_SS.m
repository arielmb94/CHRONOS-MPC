function sys = qLPV_TRMS_SS(Wh,Omh,Thth,Wv,Thtv)

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
    
    Jv = m_mr*l_m^2+m_m*(l_m^2)/3+m_cb*l_cb^2+m_b*(l_b^2)/3+m_tr*l_t^2+m_t*(l_t^2)/3+...
        0.5*m_ms*r_ms^2+m_ms*l_m^2+0.5*m_ts*r_ts^2+m_ts*l_t^2;

    
    % Matrix Coeff
    
    %dt_Wh
    a11 = -(B_tr/J_tr+(k_a^2)/(J_tr*R_a))-f1hat(Wh)/J_tr;
    b11 = k_a*k_1/(J_tr*R_a);
    
    %dt_Omh
    a21 = (l_t*f2hat(Wh)*cos(Thtv))/(K_D*(cos(Thtv)^2)+K_E*(sin(Thtv)^2)+K_F);
    a22 = -k_oh/(K_D*(cos(Thtv)^2)+K_E*(sin(Thtv)^2)+K_F);
    a23 = -f3hat(Thth)/(K_D*(cos(Thtv)^2)+K_E*(sin(Thtv)^2)+K_F);
    a24 = k_m*cos(Thtv)*(-(B_mr+(k_a^2)/R_a)-f4hat(Wv))/(J_mr*(K_D*(cos(Thtv)^2)+K_E*(sin(Thtv)^2)+K_F));
    a25 = k_m*Wv*sin(Thtv)*(K_D*(cos(Thtv)^2)-K_E*(sin(Thtv)^2)-K_F-2*K_E*(cos(Thtv)^2))/((K_D*(cos(Thtv)^2)+K_E*(sin(Thtv)^2)+K_F)^2);
    a26 = f6hat(Thtv)/(K_D*(cos(Thtv)^2)+K_E*(sin(Thtv)^2)+K_F);
    
    b22 = k_m*cos(Thtv)*(k_a*k_2/R_a)/(J_mr*(K_D*(cos(Thtv)^2)+K_E*(sin(Thtv)^2)+K_F));
    
    %dt_Thth
    a32 = 1;
    
    %dt_Wv
    a44 = -(B_mr/J_mr+(k_a^2)/(J_mr*R_a))-f4hat(Wv)/J_mr;
    
    b42 = k_a*k_2/(J_mr*R_a);
    
    %dt_Omv
    a52 = k_g*f5(Wv)*cos(Thtv)/Jv-Omh*K_H*sin(Thtv)*cos(Thtv)/Jv;
    a54 = l_m*f5hat(Wv)/Jv;
    a55 = -k_ov/Jv;
    a56 = g*((K_A-K_B)*cos(Thtv)-K_C*sin(Thtv))/(Jv*(Thtv-Thtv0));
    
    %dt_Thtv
    a65 = 1;
    
    A = [a11 0 0 0 0 0;
         a21 a22 a23 a24 a25 a26;
         0 a32 0 0 0 0;
         0 0 0 a44 0 0;
         0 a52 0 a54 a55 a56;
         0 0 0 0 a65 0];
    B = [b11 0;
         0 b22;
         0 0;
         0 b42;
         0 0;
         0 0];
     
    C = [1 0 0 0 0 0;
         0 0 1 0 0 0;
         0 0 0 1 0 0;
         0 0 0 0 0 1];
     
    sys = ss(A,B,C,0);
        
end

function y = f1hat(Wh)

    k_thp = 0.0027;
    k_thn = 0.0028;

    if Wh>=0
        y = k_thp*Wh;
    else
        y = -k_thn*Wh;
    end
end

function y = f2hat(Wh)

    k_fhp = 0.0566;
    k_fhn = 0.0660;

    if Wh>=0
        y = k_fhp*Wh;
    else
        y = -k_fhn*Wh;
    end
end

function y = f3hat(Thth)

    k_chp = 0.0158;
    k_chn = 0.0111;

    if Thth>=0
        y = k_chp;
    else
        y = k_chn;
    end
end

function y = f4hat(Wv)

    k_tvp = 0.0168;
    k_tvn = 0.0155;

    if Wv>=0
        y = k_tvp*Wv;
    else
        y = -k_tvn*Wv;
    end
end

function y = f5(Wv)

    k_fvp = 0.3819;
    k_fvn = 0.2197;

    if Wv>=0
        y = k_fvp*Wv^2;
    else
        y = -k_fvn*Wv^2;
    end
end

function y = f5hat(Wv)

    k_fvp = 0.3819;
    k_fvn = 0.2197;

    if Wv>=0
        y = k_fvp*Wv;
    else
        y = -k_fvn*Wv;
    end
end

function y = f6hat(Thtv)

    k_cvp = 0.0623;
    k_cvn = 0.0563;
    
    Thtv0 = -0.619426178368110;

    if Thtv>=Thtv0
        y = k_cvp*(Thtv-Thtv0);
    else
        y = k_cvn*(Thtv-Thtv0);
    end
end

