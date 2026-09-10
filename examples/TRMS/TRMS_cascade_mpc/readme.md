# TRMS Cascade MPC Example

## Folder contents

- *TRMS_cascade_mpc_init.m*: defines and builds the three CHRONOS controllers.
- *TRMS_cascade_mpc_sim_lpv.m*: runs the cascade in closed loop with the nonlinear plant.
- *qLPV_TRMS_cascade_mpc_SS.m*: evaluates the outer and rotor LPV models.

## Cascade architecture

This example separates the TRMS body dynamics from the rotor dynamics. An outer
MIMO MPC controls the horizontal and vertical motion by computing rotor-speed
references. Two inner SISO MPCs track those references and compute the motor
voltages applied to the nonlinear plant.

```text
                                   /--> omega_v^ref(0:N-1) --> Main-rotor MPC --> u_v
Angle references --> Outer MIMO MPC
                                   \--> omega_h^ref(0:N-1) --> Tail-rotor MPC --> u_h
```

The outer MPC provides its complete predicted rotor-speed sequences to the
inner MPCs. The inner layers can therefore use all the information predicted
by the outer layer instead of discarding everything after its first action.
Each controller can also use an LPV model and a horizon suited to its dynamics.

## Prediction models

The outer MPC uses the reduced body state

```math
x_o=\begin{bmatrix}
\Omega_h&\theta_h&\Omega_v&\widetilde\theta_v
\end{bmatrix}^T,
\qquad
\widetilde\theta_v=\theta_v-\theta_{v0},
```

and treats the rotor speeds as its control inputs:

```math
u_o=\begin{bmatrix}\omega_h^{ref}&\omega_v^{ref}\end{bmatrix}^T.
```

Its LPV model is

```math
\dot x_o=
\begin{bmatrix}
a_{22}(\rho)&a_{23}(\rho)&a_{25}(\rho)&a_{26}(\rho)\\
a_{32}&0&0&0\\
a_{52}(\rho)&0&a_{55}(\rho)&a_{56}(\rho)\\
0&0&a_{65}&0
\end{bmatrix}x_o+
\begin{bmatrix}
a_{21}(\rho)&a_{24}(\rho)\\
0&0\\
0&a_{54}(\rho)\\
0&0
\end{bmatrix}u_o+
\begin{bmatrix}b_{22}(\rho)\\0\\0\\0\end{bmatrix}d.
```

The main-rotor voltage is used as the measured disturbance $d=u_v$ to retain
its coupling with the horizontal dynamics.

The scheduling vector is
$\rho=[\omega_h,\Omega_h,\theta_h,\omega_v,\theta_v]^T$.
The coefficients marked with $(\rho)$ are evaluated from the measured TRMS
state before each set of outer and inner solves.

The inner MPCs use the rotor models

```math
\dot\omega_h=a_{11}(\rho)\omega_h+b_{11}u_h,
\qquad
\dot\omega_v=a_{44}(\rho)\omega_v+b_{42}u_v.
```

Here, $\omega_h$ and $\omega_v$ are rotor speeds in rad/s, while $u_h$ and
$u_v$ are motor voltages in V. The gains $b_{11}$ and $b_{42}$ are constant;
the scheduling-dependent rotor coefficients are evaluated at the measured
TRMS state before each set of solves.

## MPC definitions

### Outer MIMO MPC

The outer MPC tracks
$r_o=[\Omega_h^{ref},\theta_h^{ref},\Omega_v^{ref},\theta_v^{ref}-\theta_{v0}]^T$
and computes the rotor-speed reference input
$u_o=[\omega_h^{ref},\omega_v^{ref}]^T$. It solves
```math
\begin{aligned}
\min_{x_o,u_o}\quad
J_o={}&(r_{o,N}-x_{o,N})^TP_o(r_{o,N}-x_{o,N})\\
&+\frac{1}{2}\sum_{k=1}^{N}(r_{o,k}-x_{o,k})^TQ_{e,o}(r_{o,k}-x_{o,k})\\
&+\frac{1}{2}\sum_{k=0}^{N-1}\Delta u_{o,k}^TR_{\Delta u,o}\Delta u_{o,k}\\
\text{subject to}\quad
&x_{o,k+1}=A_{o,d}(\rho)x_{o,k}+B_{o,d}(\rho)u_{o,k}
+B_{d,o}(\rho)d_k, && k=0,\ldots,N-1,\\
&x_{o,\min}\leq x_{o,k}\leq x_{o,\max}, && k=1,\ldots,N,\\
&u_{o,\min}\leq u_{o,k}\leq u_{o,\max}, && k=0,\ldots,N-1.
\end{aligned}
```


### Tail-rotor MPC

The tail-rotor MPC tracks the complete rotor-speed $u_h$ sequence from the outer MPC by solving:
```math
\begin{aligned}
\min_{\omega_h,u_h}\quad
J_h={}&(\omega_{h,N}^{ref}-\omega_{h,N})^TP_h(\omega_{h,N}^{ref}-\omega_{h,N})\\
&+\frac{1}{2}\sum_{k=1}^{N}(\omega_{h,k}^{ref}-\omega_{h,k})^TQ_{e,h}
(\omega_{h,k}^{ref}-\omega_{h,k})\\
&+\frac{1}{2}\sum_{k=0}^{N-1}\Delta u_{h,k}^TR_{\Delta u,h}\Delta u_{h,k}\\
\text{subject to}\quad
&\omega_{h,k+1}=A_{h,d}(\rho)\omega_{h,k}+B_{h,d}u_{h,k}, && k=0,\ldots,N-1,\\
&-2.9\leq\omega_{h,k}\leq2.9, && k=1,\ldots,N,\\
&-2.5\leq u_{h,k}\leq2.5, && k=0,\ldots,N-1.
\end{aligned}
```

### Main-rotor MPC

The main-rotor MPC tracks the complete rotor-speed $u_v$ sequence from the outer MPC by solving:
```math
\begin{aligned}
\min_{\omega_v,u_v}\quad
J_v={}&(\omega_{v,N}^{ref}-\omega_{v,N})^TP_v(\omega_{v,N}^{ref}-\omega_{v,N})\\
&+\frac{1}{2}\sum_{k=1}^{N}(\omega_{v,k}^{ref}-\omega_{v,k})^TQ_{e,v}
(\omega_{v,k}^{ref}-\omega_{v,k})\\
&+\frac{1}{2}\sum_{k=0}^{N-1}\Delta u_{v,k}^TR_{\Delta u,v}\Delta u_{v,k}\\
\text{subject to}\quad
&\omega_{v,k+1}=A_{v,d}(\rho)\omega_{v,k}+B_{v,d}u_{v,k}, && k=0,\ldots,N-1,\\
&-1.6\leq\omega_{v,k}\leq1.6, && k=1,\ldots,N,\\
&-2\leq u_{v,k}\leq2, && k=0,\ldots,N-1.
\end{aligned}
```

## Online use

At each control sample:

1. The angle errors generate references for the TRMS body rates.
2. The outer and inner LPV models are evaluated at the measured state.
3. The outer MIMO MPC computes both predicted rotor-speed sequences, using the
   main-rotor voltage as a disturbance input.
4. Each inner MPC tracks its complete rotor-speed sequence and computes a motor
   voltage.
5. The two voltages are applied to the nonlinear TRMS plant.

## Reference

[1] Rotondo, D., Nejjari, F., & Puig, V. (2013). Quasi-LPV modeling,
identification and control of a Twin Rotor MIMO System. *Control Engineering
Practice, 21*(6), 829-846.
