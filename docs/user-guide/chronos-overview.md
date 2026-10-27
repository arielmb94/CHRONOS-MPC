# CHRONOS overview

**Prediction models, control objectives, and operating constraints in one formulation.**

CHRONOS solves Model Predictive Control (MPC) problems using LPV/LTV
prediction models. Given a model, objectives, and operating limits, CHRONOS
computes a sequence of control actions over a prediction horizon of $N$
steps. In receding-horizon operation, the first action is applied and the
problem is updated and solved again at the next control sample.

## Objective function

The objective combines the cost terms selected for the controller:

$$
\min_{x,u,v}\; J = J_{\mathrm{tracking}} + J_{\mathrm{effort}}
+ J_{\mathrm{rate}} + J_{\mathrm{custom}} + J_{\mathrm{terminal}}
+ J_{\mathrm{slack}}.
$$

| Ingredient | Purpose |
| --- | --- |
| Reference tracking | Penalize the difference between predicted outputs $y_k$ and references $r_k$. |
| Control effort | Add quadratic and/or linear costs on the control inputs $u_k$. |
| Control rate | Penalize input changes $\Delta u_k=u_k-u_{k-1}$. At the first step, use the previous applied input. |
| Custom objectives | Add linear and/or quadratic costs on user-defined affine signals. |
| Terminal objective | Penalize the final predicted state's deviation from its reference. |
| Slack penalties | Penalize relaxation of state, output, and custom bounds: $J_{\mathrm{slack}}=q_v^T v$. Here $v$ collects the feasibility slack variables and $q_v$ their penalty weights. |



## Prediction model

The discrete-time prediction model connects successive states:

$$
x_{k+1}=A_kx_k+B_ku_k+B_{d,k}d_k,
\qquad k=0,\ldots,N-1.
$$

Here $x_k$ is the state, $u_k$ the control input, and $d_k$ a known
disturbance or input. The current state $x_0$ is supplied to each solve.
Model matrices may vary along the horizon, but their values are fixed during
that solve. This represents changing operating conditions while preserving
a convex optimization problem.

Outputs are affine signals,
$y_k=C_kx_k+D_ku_k+D_{d,k}d_k$. Custom cost and constraint signals can also
combine the state, current and previous control inputs, and known inputs.

## Constraints

| Ingredient | Desired bounds | Treatment |
| --- | --- | --- |
| States | $x_{\min,k}\le x_k\le x_{\max,k}$ | Soft |
| Control inputs | $u_{\min,k}\le u_k\le u_{\max,k}$ | Hard |
| Control rates | $\Delta u_{\min,k}\le\Delta u_k\le\Delta u_{\max,k}$ | Hard |
| Outputs | $y_{\min,k}\le y_k\le y_{\max,k}$ | Soft |
| Custom affine signals | $h_{\min,k}\le h_k\le h_{\max,k}$ | Soft |

Soft bounds include penalized slack variables, allowing violations when
necessary. Their penalties control the cost of relaxing each bound. Control
and control-rate bounds have no feasibility relaxation.


## Time-varying problem data

Model and signal matrices, stage cost weights, bounds, and slack penalties
can vary across the prediction horizon. For example, a controller can account
for upcoming operating limits, adjust the balance between tracking and
effort, or use different slack penalties at different prediction steps.
Online updates change numerical data within the dimensions and ingredients
fixed when the controller is built.

The next chapter, [From nonlinear to LPV models](from-nonlinear-problems-to-lpv-mpc.md),
describes how to prepare a prediction model. The
[initialization API](../api/init/index.md) and [online API](../api/online/index.md)
document the functions used to define, update, and solve the controller.
