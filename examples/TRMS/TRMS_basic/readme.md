# Twin Rotor MIMO System (TRMS) Basic MPC Example

## Folder contents

In this folder you will find the following files:

* *TRMS_init.m*: defines and builds the CHRONOS controller.
* *TRMS_sim_lpv.m*: runs the closed-loop nonlinear-plant simulation.
* *qLPV_TRMS_SS.m*: evaluates the TRMS LPV model at the measured state.
* *compute_ref.m*: computes the remaining state references from the requested angles.

## Example overview

This is the classical centralized formulation: one six-state MIMO MPC computes
the voltages applied to the tail- and main-rotor motors. The controller follows
horizontal and vertical angle commands while accounting for the coupled,
nonlinear TRMS dynamics described in [1].

The corresponding LPV model has the structure

```math
\dot x=
\begin{bmatrix}
a_{11}(\rho) & 0 & 0 & 0 & 0 & 0 \\
a_{21}(\rho) & a_{22}(\rho) & a_{23}(\rho) & a_{24}(\rho) & a_{25}(\rho) & a_{26}(\rho) \\
0 & a_{32} & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & a_{44}(\rho) & 0 & 0 \\
0 & a_{52}(\rho) & 0 & a_{54}(\rho) & a_{55} & a_{56}(\rho) \\
0 & 0 & 0 & 0 & a_{65} & 0
\end{bmatrix}x
+\begin{bmatrix}
b_{11} & 0 \\
0 & b_{22}(\rho) \\
0 & 0 \\
0 & b_{42} \\
0 & 0 \\
0 & 0
\end{bmatrix}u
```

The MPC state uses the vertical-angle deviation
$\widetilde\theta_v=\theta_v-\theta_{v0}$:

```math
x=\begin{bmatrix}
\omega_h&\Omega_h&\theta_h&\omega_v&\Omega_v&\widetilde\theta_v
\end{bmatrix}^T,
\qquad
u=\begin{bmatrix}u_h&u_v\end{bmatrix}^T.
```

| Symbol | Meaning | Unit |
| --- | --- | --- |
| $\omega_h$, $\omega_v$ | Tail- and main-rotor angular speeds | rad/s |
| $\Omega_h$, $\Omega_v$ | Horizontal and vertical body rates | rad/s |
| $\theta_h$, $\widetilde\theta_v$ | Horizontal angle and vertical-angle deviation | rad |
| $u_h$, $u_v$ | Tail- and main-rotor motor voltages | V |

The scheduling vector is
$\rho=[\omega_h,\Omega_h,\theta_h,\omega_v,\theta_v]^T$.
*qLPV_TRMS_SS.m* evaluates $A(\rho)$ and $B(\rho)$ from the measured state.

## MPC definition

At each sample, CHRONOS solves

```math
\begin{aligned}
\min_{x,u}\quad
J={}&(r_N-x_N)^TP(r_N-x_N)
+\frac{1}{2}\sum_{k=1}^{N}(r_k-x_k)^TQ_e(r_k-x_k)\\
&+\frac{1}{2}\sum_{k=0}^{N-1}u_k^TR_u u_k
+\frac{1}{2}\sum_{k=0}^{N-1}\Delta u_k^TR_{du}\Delta u_k\\
\text{subject to}\quad
&x_{k+1}=A_k(\rho)x_k+B_k(\rho)u_k,
&&k=0,\ldots,N-1.
\end{aligned}
```

The configured bounds are
```math
\begin{bmatrix}
-2.9 \\
-1.0 \\
-1.7 \\
-1.6 \\
-0.6 \\
-0.5
\end{bmatrix}
\leq x \leq
\begin{bmatrix}
2.9 \\
1.0 \\
1.2 \\
1.6 \\
0.6 \\
1.0
\end{bmatrix}
```

```math
\begin{bmatrix}
-2.5 \\
-2.0
\end{bmatrix}
\leq u \leq
\begin{bmatrix}
2.5 \\
2.0
\end{bmatrix}
```

The input-rate bounds are

```math
\begin{bmatrix}-0.5\\-0.4\end{bmatrix}
\leq\Delta u\leq
\begin{bmatrix}0.5\\0.4\end{bmatrix}.
```

## Reference generation

The angle commands alone do not define appropriate references for the remaining
TRMS states. The six-state reference is assembled as follows:

- $\theta_h^{ref}$ and $\theta_v^{ref}$ are the requested angle trajectories.
- $\Omega_h^{ref}$ and $\Omega_v^{ref}$ use the angle error with $\tau=0.5$ s:

  $$\Omega_i^{ref}=\frac{\theta_i^{ref}-\theta_i}{\tau}.$$

- $\omega_h^{ref}$ and $\omega_v^{ref}$ are computed by *compute_ref.m* from
  the nonlinear TRMS equilibrium equations.

The reference passed to `mpc_solve` is

```math
r=\begin{bmatrix}
\omega_h^{ref}&\Omega_h^{ref}&\theta_h^{ref}&
\omega_v^{ref}&\Omega_v^{ref}&\theta_v^{ref}-\theta_{v0}
\end{bmatrix}^T.
```

This nonlinear feedforward calculation supplies the rotor-speed references. The
next TRMS example shows how CHRONOS custom costs can compute them inside the MPC
instead.


## Reference

[1] Rotondo, D., Nejjari, F., & Puig, V. (2013). Quasi-LPV modeling, identification and control of a Twin Rotor MIMO System. Control Engineering Practice, 21(6), 829-846.
