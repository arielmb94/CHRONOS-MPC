# From nonlinear to LPV models

CHRONOS does not attempt to work directly with nonlinear systems. Instead, it relies on the user providing the MPC problem ingredients (dynamics, costs, constraints, etc.) using Linear Parameter Varying (LPV) or Linear Time Varying (LTV) formulations. In exchange for this effort, the tradeoff for users can be summarized in the following *contract*:

```{note}
Give CHRONOS your dynamics and the signals used to define costs and constraints in its supported LPV/LTV forms. In return, CHRONOS turns them into a structured MPC problem and takes care of the optimization work in the background using implementations designed for maximum computational efficiency.

For a fixed LPV / LTV schedule, every optimum is global. Analytical derivatives, structured constraint reductions, and a custom KKT recursion let CHRONOS exploit the model’s structure, solving the problem without requiring generic horizon-wide dense KKT linear solvers. Stage-varying models and custom cost/constraint signals still let you capture rich model and performance behavior.
```

This guide aims to explain how to translate your nonlinear models and signals into the LPV format.


```{contents}
:local:
:depth: 2
```

## LPV models

Most system dynamics can be expressed in the form of a nonlinear ODE, expressed generally as:

$$
\dot{x}=f(x,u,d),
$$

where $x\in\mathbb{R}^{n_x}$ is the state, $u\in\mathbb{R}^{n_u}$ is the
control input, and $d\in\mathbb{R}^{n_d}$ is a known disturbance or other
known input. 

CHRONOS does not work directly with nonlinear systems, instead, it expects the dynamics to be expressed in the standard LPV form:

$$
\dot{x}=A(\rho)x+B(\rho)u+B_d(\rho)d.
$$

where $A(\rho)$, $B(\rho)$, and $B_d(\rho)$ are coefficient matrices whose elements are not fixed for ever, but can evolve over time based on a scheduling vector $\rho$. Through the process know as *nonlinear embedding*, the varying parameters $\rho$ absorb the nonlinearities and time-varying behaviour of the original nonlinear ODE model. Through this process we then obtain a linear-like representation of the original nonlinear system. 

The advantage of the LPV form is that it is input-to-output linear in the states $x$ and control inputs $u$. This results in a highly structured modeling framework, which allows us to use linear algebra to work with our nonlinear system despite the nonlinearities embedded within the $A(\rho)$, $B(\rho)$, and $B_d(\rho)$ coefficient matrices.

The same idea applies to a nonlinear signal,

$$
y=f_y(x,u,d),
$$

which could be represented in the LPV form

$$
y=C(\rho)x+D(\rho)u+D_d(\rho)d.
$$



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
update the system matrices at each time instant with
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

An exact LPV factorization can still produce a fixed-schedule model that
misses a state coupling needed to represent how the control affects the
variable of interest. In that case, locally linearize only the nonlinear
term that carries the missing coupling. The remaining dynamics can keep their
exact factorization, and the physical state and control coordinates do not
need to change.

The stirring-tank example shows this situation. Its nonlinear model is

$$
\begin{aligned}
\dot{c} &= (1-c)/\theta_f-kce^{-M/v},\\
\dot{v} &= (x_f-v)/\theta_f+kce^{-M/v}-\alpha u(v-x_c),
\end{aligned}
$$

where the coolant input $u$ acts on temperature $v$, and the concentration
$c$ depends on temperature through the reaction term $-kce^{-M/v}$. A direct
factorization gives

$$
\begin{bmatrix}\dot{c}\\\dot{v}\end{bmatrix}
=
\begin{bmatrix}
-1/\theta_f-ke^{-M/v} & 0\\
ke^{-M/v} & -1/\theta_f
\end{bmatrix}
\begin{bmatrix}c\\v\end{bmatrix}
+\begin{bmatrix}0\\-\alpha(v-x_c)\end{bmatrix}u
+\begin{bmatrix}1/\theta_f\\x_f/\theta_f\end{bmatrix}d,
\qquad d=1.
$$

This factorization is exact when the coefficient schedule uses the current
temperature. But when the schedule is fixed, the first row has no term
multiplying $v$, and $u$ enters only the temperature equation. The frozen
linear model therefore has no path from $u$ through $v$ to $c$ and is not
controllable.

The nonlinear model does contain this coupling through $-kce^{-M/v}$. To
expose it, approximate only that term with a first-order Taylor expansion
around $(c^o,v^o)$:

$$
q^o=\frac{k c^o M e^{-M/v^o}}{(v^o)^2},
\qquad
-kce^{-M/v}\approx-k e^{-M/v^o}c-q^o(v-v^o).
$$

The term $-q^o v$ adds a temperature coefficient to the concentration
equation; the offset $q^o v^o$ is carried through the known-input term. With
$x=[c\;v]^T$ and $d=[1\;v^o]^T$, the resulting hybrid LPV model is

$$
\dot{x}=
\begin{bmatrix}
-1/\theta_f-ke^{-M/v^o} & -q^o\\
ke^{-M/v^o} & -1/\theta_f
\end{bmatrix}x
+\begin{bmatrix}0\\-\alpha(v^o-x_c)\end{bmatrix}u
+\begin{bmatrix}
1/\theta_f & q^o\\
x_f/\theta_f & 0
\end{bmatrix}d.
$$

The temperature-to-concentration coefficient now connects the state affected
by the control to the concentration equation. Only the term needed to expose
this path is approximated; the other terms retain their factorized form. In
the repository example, $(c^o,v^o)$ is the measured reactor state at the
current control sample. See the
[stirring-tank example README](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Stirring%20Tank%20Reactor/readme.md)
and [`stirring_tank_init.m`](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Stirring%20Tank%20Reactor/stirring_tank_init.m).

More generally, a first-order expansion of the full vector field at
$(x^o,u^o)$ gives

$$
f(x,u)\approx f(x^o,u^o)+A_o(x-x^o)+B_o(u-u^o)
=A_o x+B_o u+c_o,
$$

where $A_o=\partial f/\partial x$, $B_o=\partial f/\partial u$, and

$$
c_o=f(x^o,u^o)-A_o x^o-B_o u^o.
$$

At an equilibrium, $f(x^o,u^o)=0$, so
$c_o=-A_o x^o-B_o u^o$. This affine offset can be represented by
CHRONOS's known-input term, for example $d=1$ and $B_{d,o}=c_o$. This form
keeps the original state and control coordinates rather than introducing
deviation variables. Partial and full Taylor expansions are local
approximations whose accuracy is expected near the chosen operating point.

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
