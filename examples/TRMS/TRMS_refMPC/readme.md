# TRMS Virtual-Reference MPC Example

## What this example demonstrates

This example shows how a CHRONOS custom cost can connect quantities that have
different roles in the prediction model. The controller simultaneously chooses

- the physical motor voltages that change the real rotor speeds, and
- virtual rotor-speed references that drive the predicted body motion.

The custom cost softly ties each virtual reference to its corresponding
physical rotor speed. This lets the MPC find useful rotor-speed references as
part of the optimization instead of computing them beforehand from nonlinear
equilibrium equations.

The key idea is not specific to the TRMS. CHRONOS can penalize a user-defined
linear performance signal containing states, current or previous control
actions, and known online signals. Application-specific objectives can therefore
be added without changing the solver.

## Folder contents

- *TRMS_refMPC_init.m*: defines and builds the CHRONOS controller.
- *TRMS_refMPC_sim_lpv.m*: runs the closed-loop nonlinear-plant simulation.
- *qLPV_TRMS_refMPC_SS.m*: evaluates the modified LPV model at the measured
  state.

## From external reference generation to joint optimization

The Basic TRMS example and this example control the same six-state plant, but
they generate the rotor-speed references differently.

| | Basic TRMS | Virtual-reference TRMS |
| --- | --- | --- |
| MPC inputs | Two motor voltages | Two motor voltages and two virtual rotor-speed references |
| Tracked outputs | All six states | Body rates and angles only |
| Rotor-speed references | Computed before the solve by *compute_ref.m* | Chosen inside the MPC optimization |
| Rotor-to-body coupling | Physical rotor states appear directly in the body model | Virtual rotor-speed inputs appear in the body model |
| Link to physical rotor speeds | Directly through the original state equations | Softly enforced by the custom cost |

In the Basic example, `compute_ref.m` evaluates nonlinear equilibrium
relationships to obtain $\omega_h^{ref}$ and $\omega_v^{ref}$ from the angle
commands. Those values become entries of the six-state tracking reference.

Here, the standard tracking reference contains only body rates and angles. The
two rotor-speed references become additional MPC decision variables, so the
prediction can choose the profiles that best serve angle tracking while
respecting the rotor, voltage, and rate limits.

## Prediction model and decision variables

The physical plant state is

$$
x=\begin{bmatrix}
\omega_h&\Omega_h&\theta_h&\omega_v&\Omega_v&\widetilde\theta_v
\end{bmatrix}^T,
\qquad
\widetilde\theta_v=\theta_v-\theta_{v0}.
$$

The MPC input is extended to

$$
u_{MPC}=\begin{bmatrix}
u_h&u_v&\omega_h^{ref}&\omega_v^{ref}
\end{bmatrix}^T.
$$

| Variable | Role | Unit |
| --- | --- | --- |
| $u_h$, $u_v$ | Physical tail- and main-rotor motor voltages | V |
| $\omega_h$, $\omega_v$ | Predicted physical rotor states | rad/s |
| $\omega_h^{ref}$, $\omega_v^{ref}$ | Virtual inputs used by the predicted body dynamics | rad/s |

Only $u_h$ and $u_v$ are applied to the nonlinear plant. The virtual inputs
exist only inside the controller.

The modified model separates the rotor dynamics from their effect on the body:

$$
\dot x=
\begin{bmatrix}
a_{11}(\rho) & 0 & 0 & 0 & 0 & 0\\
0 & a_{22}(\rho) & a_{23}(\rho) & 0 & a_{25}(\rho) & a_{26}(\rho)\\
0 & a_{32} & 0 & 0 & 0 & 0\\
0 & 0 & 0 & a_{44}(\rho) & 0 & 0\\
0 & a_{52}(\rho) & 0 & 0 & a_{55} & a_{56}(\rho)\\
0 & 0 & 0 & 0 & a_{65} & 0
\end{bmatrix}x+
\begin{bmatrix}
b_{11} & 0 & 0 & 0\\
0 & b_{22}(\rho) & a_{21}(\rho) & a_{24}(\rho)\\
0 & 0 & 0 & 0\\
0 & b_{42} & 0 & 0\\
0 & 0 & 0 & a_{54}(\rho)\\
0 & 0 & 0 & 0
\end{bmatrix}u_{MPC}.
$$

The scheduling vector is
$\rho=[\omega_h,\Omega_h,\theta_h,\omega_v,\theta_v]^T$.
*qLPV_TRMS_refMPC_SS.m* evaluates the scheduling-dependent coefficients at
the measured state.

This separation creates a useful optimization structure:

1. Body tracking determines which virtual rotor-speed profiles would produce
   the requested motion.
2. The motor voltages determine which physical rotor-speed profiles the rotor
   dynamics can produce.
3. The custom cost penalizes disagreement between the two.

## Why the custom cost is necessary

If the virtual references were free, the optimizer could predict body motion
using rotor speeds that the voltage-driven physical rotors do not achieve. The
body prediction would then benefit from fictitious actuation.

The example defines the two-component performance signal

$$
z_k=
\begin{bmatrix}
\omega_{h,k}^{ref}-\omega_{h,k}\\
\omega_{v,k}^{ref}-\omega_{v,k}
\end{bmatrix}
=C_zx_k+D_zu_k
$$

and adds the following term to the cost function:

$$
J_z=\frac{1}{2}\sum_{k=0}^{N-1}z_k^TQ_z z_k
$$

This is a soft consistency condition. A larger $Q_z$ forces the virtual and
physical speeds to agree more closely; a smaller $Q_z$ gives the virtual
references more freedom to improve body tracking. The mismatch need not be
zero because the optimizer trades this cost against the tracking, input-rate,
terminal, and constraint terms.

The initializer constructs the signal directly:

```matlab
% z1 = omega_h_ref - omega_h
% z2 = omega_v_ref - omega_v
Cz = [-1 0 0  0 0 0;
       0 0 0 -1 0 0];

Dz = [0 0 1 0;
      0 0 0 1];

Dsuz = [];
Ddz  = [];
Qz   = diag([100 100]);

mpc = init_mpc_Custom_cost(mpc,Cz,Dz,Dsuz,Ddz,Qz);
```

The first row of `Cz` selects $-\omega_h$ and the first row of `Dz` selects
$\omega_h^{ref}$. The second rows do the same for the main rotor. Therefore,
`Cz*x + Dz*u` is exactly the virtual-to-physical speed error.

## Complete MPC objective

The standard tracking output is

$$
y=\begin{bmatrix}
\Omega_h&\theta_h&\Omega_v&\widetilde\theta_v
\end{bmatrix}^T.
$$

With the configuration in *TRMS_refMPC_init.m*, CHRONOS solves

$$
\begin{aligned}
\min_{x,u}\quad
J={}&(x_N^{ref}-x_N)^TP(x_N^{ref}-x_N)
+\frac{1}{2}\sum_{k=1}^{N}(r_k-y_k)^TQ_e(r_k-y_k)\\
&+\frac{1}{2}\sum_{k=0}^{N-1}\Delta u_k^TR_{du}\Delta u_k
+\frac{1}{2}\sum_{k=0}^{N-1}z_k^TQ_z z_k\\
\text{subject to}\quad
&x_{k+1}=A_k(\rho)x_k+B_k(\rho)u_k,
&&k=0,\ldots,N-1,\\
&x_{\min}\leq x_k\leq x_{\max}, &&k=1,\ldots,N,\\
&u_{\min}\leq u_k\leq u_{\max}, &&k=0,\ldots,N-1,\\
&\Delta u_{\min}\leq\Delta u_k\leq\Delta u_{\max},
&&k=0,\ldots,N-1.
\end{aligned}
$$

The bounds on the four MPC inputs have two meanings: the first two limit motor
voltages, while the last two keep the virtual references within the physical
rotor-speed range. The rate cost and rate bounds also act on all four entries,
so both voltage commands and virtual-reference profiles are kept smooth.

## Online sequence and expected result

At each control sample:

1. The angle errors generate references for $\Omega_h$ and $\Omega_v$.
2. The modified LPV matrices are evaluated at the measured TRMS state.
3. `mpc_solve` jointly optimizes the two motor voltages and two virtual
   rotor-speed references.
4. Only the two voltages are applied to the nonlinear plant.


## Reference

[1] Rotondo, D., Nejjari, F., & Puig, V. (2013). Quasi-LPV modeling,
identification and control of a Twin Rotor MIMO System. *Control Engineering
Practice, 21*(6), 829-846.
