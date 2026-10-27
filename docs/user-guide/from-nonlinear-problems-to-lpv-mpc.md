# From nonlinear to LPV-MPC

**Objective:** learn how to express nonlinear dynamics and signals in the
linear parameter-varying (LPV) form used by CHRONOS.

**Prerequisites:** familiarity with state-space models, MPC, and MATLAB.

**Expected outcome:** choose an exact embedding or a local approximation,
identify its scheduling variables, and obtain the discrete matrices and
signals needed to define a CHRONOS controller.

```{contents}
:local:
:depth: 2
```

## CHRONOS model and notation

Nonlinear continuous-time dynamics are commonly written as

$$
\dot{x}=f(x,u,d),
$$

where $x\in\mathbb{R}^{n_x}$ is the state, $u\in\mathbb{R}^{n_u}$ is the
control input, and $d\in\mathbb{R}^{n_d}$ is a known disturbance or other
known input. To use this model with CHRONOS, express its dynamics in LPV form:

$$
\dot{x}=A(\rho)x+B(\rho)u+B_d(\rho)d.
$$

Here $A(\rho)$, $B(\rho)$, and $B_d(\rho)$ are coefficient matrices selected
by the scheduling variable or vector $\rho$. For a fixed $\rho$, the model is
linear in $x$, $u$, and $d$.

The same idea applies to a nonlinear output or signal,

$$
y=f_y(x,u,d),
$$

which CHRONOS can represent in the form

$$
y=C(\rho)x+D(\rho)u+D_d(\rho)d.
$$

CHRONOS represents custom performance and constraint signals with the
following affine forms:

$$
\begin{aligned}
y_k &= C_k x_k+D_k u_k+D_{d,k}d_k,\\
z_k &= C_{z,k}x_k+D_{z,k}u_k+D_{xu,z,k}xu_k+D_{dz,k}d_{z,k},\\
h_k &= C_{h,k}x_k+D_{h,k}u_k+D_{xu,h,k}xu_k+D_{dh,k}d_{h,k}.
\end{aligned}
$$

The dimensions of $y_k$, $z_k$, and $h_k$ are $n_y$, $n_z$, and $n_h$,
respectively. The known inputs satisfy $d_z\in\mathbb{R}^{n_{dz}}$ and
$d_h\in\mathbb{R}^{n_{dh}}$; they are specific to the custom-cost and
custom-constraint signals. When a control-rate or custom-signal term uses
$xu_k$, it denotes the control action preceding $u_k$. Once the scheduling
values are fixed, each supported signal is affine in its decision variables.

| Signal | Coefficient dimensions |
| --- | --- |
| $y$ | $C:n_y\times n_x$; $D:n_y\times n_u$; $D_d:n_y\times n_d$ |
| $z$ | $C_z:n_z\times n_x$; $D_z,D_{xu,z}:n_z\times n_u$; $D_{dz}:n_z\times n_{dz}$ |
| $h$ | $C_h:n_h\times n_x$; $D_h,D_{xu,h}:n_h\times n_u$; $D_{dh}:n_h\times n_{dh}$ |

During a solve, scheduling values are fixed for each stage. One value can be
shared across the prediction horizon, or the values can vary by stage. If
$\rho$ depends on the state, its values are available from CHRONOS's internal
state trajectory; an external forecast is another option. Matrix coefficients
are evaluated at the corresponding scheduling values. Later guide pages cover
scheduling and update details.

## Three ways to represent nonlinear equations

The following sections show three ways to derive LPV coefficients for
nonlinear dynamics and signals: direct factorization, algebraic manipulation,
and Taylor linearization.

### 1. Direct factorization

If the equations already collect the states, controls, and known inputs
linearly, write

$$
\dot{x}=A_c(\rho)x+B_c(\rho)u+E_c(\rho)d.
$$

This is an exact factorization: the nonlinearity is in the coefficients as
functions of the schedule $\rho$, rather than in products among the modeled
state and control variables after the schedule is fixed.

For a scalar example,

$$
\dot{x}=-x^2+u=(-x)x+u=A(\rho)x+Bu,
\qquad \rho=x,\quad A(\rho)=-\rho,\quad B=1.
$$

For a fixed scheduling value $\rho$, $A(\rho)$ is fixed and the model is
linear in the state and input. Since $\rho=x$ in this example, CHRONOS can
update the system matrices at each time instance with
[`update_mpc_dynamics`](../api/online/update_mpc_dynamics.md).

The [scaled autonomous vehicle example](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Scaled%20Autonomous%20Vehicle/readme.md)
uses the lateral bicycle model

$$
\begin{bmatrix}\dot v_y\\\ddot\psi\end{bmatrix}
=A_c(v_x)
\begin{bmatrix}v_y\\\dot\psi\end{bmatrix}
+B_c(v_x)\delta,
$$

where longitudinal speed $v_x$ schedules the matrices. The example evaluates
them in [`update_BM.m`](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Scaled%20Autonomous%20Vehicle/update_BM.m)
and initializes the prediction model in
[`sav_lat_controller_init.m`](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Scaled%20Autonomous%20Vehicle/sav_lat_controller_init.m).

### 2. Algebraic manipulation

If a term is not directly factorizable, an exact algebraic identity may expose
the state or input factor needed for a matrix representation. For the
two-tank model,

$$
\begin{aligned}
\dot h_1 &= u/A_b-\sqrt{2gh_1}/A_b,\\
\dot h_2 &= \sqrt{2gh_1}/A_b-\sqrt{2gh_2}/A_b.
\end{aligned}
$$

For positive tank heights,

$$
\sqrt{2gh_i}=\frac{\sqrt{2gh_i}}{h_i}h_i,
\qquad h_i>0.
$$

Thus, with $x=[h_1\;h_2]^T$ and $\rho=[h_1\;h_2]^T$, the same equations
can be written exactly on that domain as

$$
\dot{x}=
\begin{bmatrix}
-\dfrac{\sqrt{2gh_1}}{A_bh_1} & 0\\
\dfrac{\sqrt{2gh_1}}{A_bh_1} & -\dfrac{\sqrt{2gh_2}}{A_bh_2}
\end{bmatrix}x
+\begin{bmatrix}1/A_b\\0\end{bmatrix}u.
$$

The coefficient expression is undefined at $h_i=0$, so this representation
requires strictly positive heights. The [two-tank example](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Two%20Tank/readme.md)
uses this embedding; its initializer and LPV simulation are
[`two_tank_init.m`](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Two%20Tank/two_tank_init.m)
and [`sim_two_tank_lpv.m`](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Two%20Tank/sim_two_tank_lpv.m).

### 3. Taylor linearization

When exact factorization is not useful, approximate only the nonfactorable
terms around a local operating point. For example, the stirring-tank model
approximates the reaction term $-kce^{-M/v}$ around $(c^o,v^o)$:

$$
q^o=\frac{k c^o M e^{-M/v^o}}{(v^o)^2},
\qquad
-kce^{-M/v}\approx-k e^{-M/v^o}c-q^o(v-v^o).
$$

The rest of the model can retain its exact LPV embedding. This termwise
approximation introduces a coefficient multiplying $v$ while keeping the
physical states $c$ and $v$. See the
[stirring-tank example](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Stirring%20Tank%20Reactor/readme.md)
and its [`stirring_tank_init.m`](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Stirring%20Tank%20Reactor/stirring_tank_init.m).

Alternatively, apply a first-order Taylor expansion to the full vector field
at $(x^o,u^o)$:

$$
f(x,u)\approx f(x^o,u^o)+A_o(x-x^o)+B_o(u-u^o)
=A_o x+B_o u+c_o,
$$

where $A_o=\partial f/\partial x$ and $B_o=\partial f/\partial u$ at the
operating point, and

$$
c_o=f(x^o,u^o)-A_o x^o-B_o u^o.
$$

This keeps the original state and control coordinates; it does not require
the controller to use deviation variables $x-x^o$ and $u-u^o$. Put the
constant offset in CHRONOS's known-input term, for example

$$
\dot{x}\approx A_o x+B_o u+B_{d,o}d,
\qquad d=1,\quad B_{d,o}=c_o.
$$

For an offset decomposed into several known quantities, use a matching vector
$d$ and matrix $B_d$. Taylor expansion of selected terms and Taylor expansion
of the full vector field are both local approximations: their accuracy is
expected near the chosen operating point, not over the entire state space.

## Discretize the model

CHRONOS uses a discrete-time prediction model. For a continuous LPV model

$$
\dot{x}=A_c(\rho)x+B_c(\rho)u+E_c(\rho)d,
$$

the discrete prediction model is

$$
x_{k+1}=A_kx_k+B_ku_k+B_{d,k}d_k,
\qquad k=0,\ldots,N-1.
$$

Under forward Euler with sample time $T_s$, the coefficients are

$$
A_k=I+T_sA_{c,k},\qquad
B_k=T_sB_{c,k},\qquad
B_{d,k}=T_sE_{c,k}.
$$

This assumes the coefficients and known input are held over the sample
interval. Discretize the state, control, and known-input terms consistently.
Other discretization methods can be used when appropriate to the plant and
sampling assumptions.

## Map the model to CHRONOS

Use [`init_mpc_dynamics`](../api/init/init_mpc_dynamics.md) to define the
discrete $A$, $B$, and optional $B_d$ matrices. The known input $d_k$ is
supplied when building and solving the controller. For changing coefficients,
[`update_mpc_dynamics`](../api/online/update_mpc_dynamics.md) accepts new
constant or stage-varying matrices without changing their initialized
dimensions.

For a nonlinear output $y=g(x,u,d)$, first factor or approximate it in the
form $Cx+Du+D_dd$, then use [`init_mpc_output`](../api/init/init_mpc_output.md).
After building, update staged output coefficients with
[`update_mpc_output_vector`](../api/online/update_mpc_output_vector.md). The
resulting output can be tracked or bounded. A custom performance signal $z$
can be penalized by a quadratic or linear custom cost, and a custom signal $h$
can be bounded with a custom constraint; see
[`init_mpc_Custom_cost`](../api/init/init_mpc_Custom_cost.md) and
[`init_mpc_custom_cnstr`](../api/init/init_mpc_custom_cnstr.md).

CHRONOS supports quadratic tracking, control, control-rate, custom, and
terminal costs, along with linear bounds on state, control, control rate,
output, and custom signals. An arbitrary nonlinear cost or constraint
function is not passed directly to the solver: express it using these
supported signals and terms, and use an exact representation or a stated
approximation for any nonlinear mapping. See the API pages for
[tracking](../api/init/init_mpc_Tracking_cost.md),
[control](../api/init/init_mpc_Control_cost.md),
[control-rate](../api/init/init_mpc_ControlRate_cost.md),
[state bounds](../api/init/init_mpc_state_cnstr.md),
[control bounds](../api/init/init_mpc_control_cnstr.md),
[control-rate bounds](../api/init/init_mpc_control_rate_cnstr.md), and
[output bounds](../api/init/init_mpc_output_cnstr.md).

Return to the [User Guide index](index.md).
