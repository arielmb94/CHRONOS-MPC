function [Ac, Bc, Mc] = qLPV_Furuta_SS(x)
    % State vector: x = [theta1; theta2; dtheta1; dtheta2]
    % theta2 = 0 is the upright equilibrium, theta2 = pi is downward resting
    
    th1  = x(1);
    th2  = x(2);
    dth1 = x(3);
    dth2 = x(4);

    % Physical parameters
    g   = 9.81;
    m2  = 0.05;     % Pendulum mass (kg)
    L1  = 0.215;    % Rotary arm length (m)
    l2  = 0.15;     % Pendulum center of mass distance (m)
    I1  = 0.002;    % Arm moment of inertia (kg.m^2)
    I2  = 0.0004;   % Pendulum moment of inertia (kg.m^2)
    mu1 = 0.005;    % Arm viscous friction (N.m.s/rad)
    mu2 = 0.0005;   % Pendulum viscous friction (N.m.s/rad)

    % Alison parameters
    % m2 = 0.270/3;
    % L1 = 0.170/2;
    % l2 = 0.130;
    % I1 = 0.005;
    % I2 = 5.07e-4;
    % mu1 = 0.001;
    % mu2 = 0.001;

    % Lumped mechanical parameters
    alpha = I1 + m2 * L1^2;
    beta  = m2 * l2^2;
    gamma = m2 * L1 * l2;
    delta = I2 + m2 * l2^2;
    eta   = m2 * g * l2;

    % Unnormalised sinc function: sin(th2)/th2
    if abs(th2) < 1e-6
        sinc_th2 = 1.0;
    else
        sinc_th2 = sin(th2) / th2;
    end

    % Inertia Matrix M(q)
    Mq = [alpha + beta * sin(th2)^2, gamma * cos(th2);
          gamma * cos(th2),          delta];

    % Coriolis + Viscous Friction: N(q, dq) + F
    NF = [mu1 + beta * dth2 * sin(2*th2),  -gamma * dth2 * sin(th2);
         -0.5 * beta * dth1 * sin(2*th2),   mu2];

    % Gravity matrix -S(th2) such that -d(q) = [0; eta*sinc(th2)*th2]
    minus_S = [0, 0;
               0, eta * sinc_th2];

    % Assemble full 4x4 Continuous Descriptor System: Mc * dx = Ac * x + Bc * u
    Mc = [eye(2),   zeros(2,2);
          zeros(2,2), Mq];

    Ac = [zeros(2,2), eye(2);
          minus_S,   -NF];

    Bc = [0; 0; 1; 0];
end