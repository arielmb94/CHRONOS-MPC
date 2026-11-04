# Configuring the MPC controller

[From nonlinear to LPV models](from-nonlinear-problems-to-lpv-mpc.md) explains
how to obtain a discrete prediction model and additional signal models in
CHRONOS's supported form. This chapter shows how to register those models,
choose objectives and limits, and build the controller through the
initialization API. The [API reference](../api/init/index.md) provides the
full argument descriptions.

Define and build the controller once, then reuse it at each sampling instant.
The next chapter, [Solving and updating the MPC](solving-and-updating-the-mpc.md), explains
how to update configured data and solve with the latest state, references,
and forecasts.

```{contents}
:local:
:depth: 2
```

## Time-varying by design

CHRONOS was designed from the outset for time-varying MPC problems. It supports [linear parameter-varying / linear time-varying (LPV/LTV)](from-nonlinear-problems-to-lpv-mpc.md) prediction models, with model coefficients, cost weights, and constraint bounds that can vary from stage to stage along the prediction horizon.

For example, evaluating an LPV model at a sequence of operating points gives an LTV prediction model with a local linear model at each stage. This lets CHRONOS approximate nonlinear plant behavior along the prediction horizon.

The MPC data are held fixed during each solve, so the prediction dynamics are linear in the state and control decision variables. With affine constraints and a convex quadratic objective, CHRONOS can use very efficient convex optimization algorithms to solve the resulting MPC problem.

The following data can vary across stages of the prediction horizon:

*   **Model and signal coefficients** for prediction dynamics, outputs, custom costs, and custom constraints.
*   **Known exogenous inputs** to the dynamics and outputs, custom costs, and custom constraints.
*   **Cost weights**, including tracking, control effort, control rate, and custom cost weights.
*   **Tracking reference trajectory**.
*   **Lower and upper constraint bounds**.
*   **Soft-constraint slack penalties** for constraints on states, outputs, and custom constraint signals.

### Supplying time-varying data

CHRONOS time-varying data can be vectors, such as constraint-bound vectors, or matrices, such as prediction-model coefficients. CHRONOS handles time variation differently depending on the data type. A vector uses its **second dimension** for time: each column contains one stage's values. A matrix uses its **third dimension**: each page contains one stage's coefficients.

| Data | Constant | Time-varying  |
| --- | --- | --- |
| Vector (`n` components) | `n`-by-1 | `n`-by-`L` |
| Matrix (`nr` rows by `nc` columns) | `nr`-by-`nc` | `nr`-by-`nc`-by-`L` |

Here `L` is the number of supplied stages. A single column or matrix is reused across the horizon. If $L$ is shorter than the prediction horizon $N$, CHRONOS automatically repeats the final column or page; excess stages are ignored.

### Prediction stages

CHRONOS separates in the background the prediction horizon into three stage regions depending on which optimization variables are accessible at each one.

| Stage | Interpretation |
| --- | --- |
| $k=0$ | $u_0$ is optimized, the current state `x_prev` and previous input `u_prev` are known. |
| $k=1,\ldots,N-1$ | Predicted states and controls are optimized. |
| $k=N$ | Only the terminal state is optimized. |


## MPC formulation

CHRONOS solves the following MPC problem:

$$
\begin{aligned}
\min_{x,u}\; J ={}&
\frac12\sum_{k\in\mathcal K_y}(r_k-y_k)^TQ_{e,k}(r_k-y_k)\\
&+\sum_{k=0}^{N-1}\left(\frac12u_k^TR_{u,k}u_k+r_{u,k}^Tu_k\right)\\
&+\frac12\sum_{k=0}^{N-1}\Delta u_k^TR_{du,k}\Delta u_k\\
&+\sum_{k\in\mathcal K_z}\left(\frac12z_k^TQ_{z,k}z_k+q_{z,k}^Tz_k\right)\\
&+(x_{N,\mathrm{ref}}-x_N)^TP(x_{N,\mathrm{ref}}-x_N)
\end{aligned}
$$

with prediction dynamics and desired operating bounds:

$$
\begin{aligned}
x_{k+1}&=A_kx_k+B_ku_k+B_{d,k}d_k, && k=0,\ldots,N-1,\\
x_{\min,k}&\le x_k\le x_{\max,k}, && k=1,\ldots,N,\\
u_{\min,k}&\le u_k\le u_{\max,k}, && k=0,\ldots,N-1,\\
\Delta u_{\min,k}&\le\Delta u_k\le\Delta u_{\max,k}, && k=0,\ldots,N-1,\\
y_{\min,k}&\le y_k\le y_{\max,k}, && k\in\mathcal K_y,\\
h_{\min,k}&\le h_k\le h_{\max,k}, && k\in\mathcal K_h.
\end{aligned}
$$

| Symbol | Meaning |
| --- | --- |
| $N$ |Prediction horizon length |
| $x_k\in\mathbb R^{n_x}$ |  State vector|
| $u_k\in\mathbb R^{n_u}$ | Computed control input |
| $d_k\in\mathbb R^{n_d}$ | Known inputs |
| $y_k\in\mathbb R^{n_y}$ | Output signals |
| $r_k\in\mathbb R^{n_y}$ | Tracking reference |
| $z_k\in\mathbb R^{n_z}$ | Custom cost signals |
| $h_k\in\mathbb R^{n_h}$ | Custom constraint signals |


## CHRONOS MPC definition workflow

To define your MPC problem in CHRONOS, follow these steps:


1.  Initialize the CHRONOS MPC structure
2.  Configure the ingredients your MPC requires
3.  Call the builder


Always start by initializing a CHRONOS MPC structure:

```matlab
mpc = init_chronos_mpc(N);
```

[`init_chronos_mpc(N)`](../api/init/init_chronos_mpc.md) selects the horizon and creates the structure to be populated by subsequent calls. The supported horizon is an
integer $N\ge1$; omitting the argument selects $N=10$. For a model sampled
every $T_s$ seconds, $N$ prediction steps cover $NT_s$ seconds.

### The MPC structure

The `mpc` structure holds the controller's key data, including its prediction
model, configured costs and constraints, and the solver data needed during a
solve. During definition, ingredient initialization calls populate and enable
the fields required by your application. The call to `build_chronos_mpc` then
computes and allocates the fixed-size solver data and workspace for the
configured ingredients and initializes the predicted trajectories.

During online operation, you can update numerical data between
solves. However, after the build, the horizon, state and input dimensions, and enabled ingredients are fixed. They cannot be changed online because CHRONOS prepares its solver data and workspace for that structure and those dimensions; redefine and rebuild the controller to make structural changes.

### Ingredient initialization order

Initialize the MPC ingredients in the following order:

1.  Dynamics model
2.  Output model (optional)
3.  Constraints and costs

The following sections describe each of the ingredients that help you define your MPC problem.

## Prediction model

### Set up the dynamics

Pass the discrete-time model at your controller's sample time to
[`init_mpc_dynamics`](../api/init/init_mpc_dynamics.md).

A single matrix is enough for a coefficient that stays constant across the
horizon. If a known input such as a disturbance forecast affects the state
prediction, include its coefficient `Bd` and supply its values as `d_in` when building and solving. CHRONOS predicts the next state at each stage using

$$
x_{k+1}=A_kx_k+B_ku_k+B_{d,k}d_k.
$$

In the coefficient tables below, **Size** is the stored field size. You can
supply a constant or a profile with $L$ stages.

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.A` | State transition matrix | $n_x\times n_x\times N$ | $k=0,\ldots,N-1$ |
| `mpc.B` | Control-input matrix | $n_x\times n_u\times N$ | $k=0,\ldots,N-1$ |
| `mpc.Bd` | Known-input matrix | $n_x\times n_d\times N$ | $k=0,\ldots,N-1$ |

**Associated signals**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.x` | Predicted state | $n_x\times N$ | $k=1,\ldots,N$ |
| `mpc.u` | Optimized control input | $n_u\times N$ | $k=0,\ldots,N-1$ |
| `mpc.d` | Known input supplied through `d_in` | $n_d\times N$ | $k=0,\ldots,N-1$ |

### Choose a custom output model

By default, the call to [`init_mpc_dynamics`](../api/init/init_mpc_dynamics.md) sets the full state as the output, $y=x$. To track a subset of states or a different quantity, define an output with
[`init_mpc_output`](../api/init/init_mpc_output.md).

Supply at least one nonzero coefficient among `C` and `D`. Use `C` for state
dependence and `D` for direct control-input dependence; a control-only output
can use `C=[]`. If the output also
depends on a known input, `Dd` uses the same $d_k$ vector as the dynamics. The output at each stage is

$$
y_k=C_kx_k+D_ku_k+D_{d,k}d_k.
$$

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.C` | State-to-output matrix | $n_y\times n_x\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.D` | Control-to-output matrix | $n_y\times n_u\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.Dd` | Known-input-to-output matrix | $n_y\times n_d\times (N-1)$ | $k=1,\ldots,N-1$ |

**Associated signals**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.y` | Output computed from the configured output model | $n_y\times (N-1)$ | $k=1,\ldots,N-1$ |

If the output model uses an active `mpc.D` field, CHRONOS will also compute an output signal for $k=0$ stored in `mpc.y_0`. Similarly, if the output model has rows that depend exclusively on the states, CHRONOS will compute an output signal for the stage $k=N$ stored in `mpc.y_ter`.

## Cost terms

Combine the optional costs that reflect what you want the controller to achieve. Use symmetric positive semidefinite matrices for the stage quadratic weights.

| Cost type | When to use it |
| --- | --- |
| {ref}`Tracking <mpc-cost-tracking>` | Keep predicted outputs close to a target or reference trajectory. |
| {ref}`Control effort <mpc-cost-control-effort>` | Penalize large control inputs or account for a linear cost per unit of input. |
| {ref}`Control rate <mpc-cost-control-rate>` | Discourage large changes in the control input between samples. |
| {ref}`Terminal state <mpc-cost-terminal-state>` | Bring the final predicted state close to a full-state target. |
| {ref}`Custom cost <mpc-cost-custom>` | Penalize a signal you define from the state, current or previous input, or known signals. |

(mpc-cost-tracking)=
### Track a reference

To keep the configured output close to a target, set its error weight `Qe`
with [`init_mpc_Tracking_cost`](../api/init/init_mpc_Tracking_cost.md).
Supply one `ny`-by-`ny` matrix for a constant weight or a profile with one
matrix per supplied stage. For example, increasing the weight in later
pages makes errors further ahead count more.

Pass the target as `r_in` when calling `mpc_solve`. Use an `ny`-by-1 column
for a constant target or an `ny`-by-`L` reference trajectory. The output
model determines which terms contribute at the current and terminal stages.

$$
J_{\mathrm{tracking}}=\frac12\sum_{k\in\mathcal K_y}
(r_k-y_k)^TQ_{e,k}(r_k-y_k).
$$

**Cost weights**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.Qe` | Quadratic weight on the tracking error | $n_y\times n_y\times (N-1)$ | $k=1,\ldots,N-1$ |

**Associated signals**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.r` | Tracking reference supplied through `r_in` | $n_y\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.y` | Output computed from the configured output model | $n_y\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.err` | Tracking error, `mpc.r - mpc.y` | $n_y\times (N-1)$ | $k=1,\ldots,N-1$ |

CHRONOS computes `mpc.y` from the output model and stores `mpc.r - mpc.y`
in `mpc.err`.

If the output model uses an active `mpc.D` field, CHRONOS will automatically populate the data required to use the output at $k=0$. Similarly, if the output model has rows that depend exclusively on the states, CHRONOS will automatically populate the data required to use the output at $k=N$.

(mpc-cost-control-effort)=
### Penalize control effort

You can use quadratic or linear costs on the control action, or both. Use `Ru` to penalize control effort quadratically. As the control signal
grows, the penalty grows with the square of its
magnitude. Use `ru` for a linear cost. Each unit increase in a control input changes the cost by the same amount regardless of the current input value. Supply `Ru`, `ru`, or both to [`init_mpc_Control_cost`](../api/init/init_mpc_Control_cost.md).

A (positive) entry of `ru` favors smaller values of the corresponding input. When using a linear control cost, use nonnegative input bounds when modeling a quantity that cannot be negative, such as a one-way flow. Otherwise, the linear cost `ru` will push the control action toward negative values, although the resulting control balances this linear term against the other costs and constraints.

$$
J_{\mathrm{effort}}=\sum_{k=0}^{N-1}
\left(\frac12u_k^TR_{u,k}u_k+r_{u,k}^Tu_k\right).
$$

Supply `Ru` as a matrix or a profile of matrix pages, and `ru` as a vector
or a profile of columns. The tables show the stored sizes.

**Cost weights**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.Ru` | Quadratic weight on the control input | $n_u\times n_u\times N$ | $k=0,\ldots,N-1$ |
| `mpc.ru` | Linear weight on the control input | $n_u\times N$ | $k=0,\ldots,N-1$ |

**Associated signals**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.u` | Control inputs penalized by the cost | $n_u\times N$ | $k=0,\ldots,N-1$ |

(mpc-cost-control-rate)=
### Smooth changes in control

To discourage abrupt changes in the control input, set `Rdu` with
[`init_mpc_ControlRate_cost`](../api/init/init_mpc_ControlRate_cost.md).
Use a single `nu`-by-`nu` matrix or an `nu`-by-`nu`-by-`L` profile.
The weight applies to the change per sample, with the first change measured
from the previously applied input.

$$
J_{\mathrm{rate}}=\frac12\sum_{k=0}^{N-1}
\Delta u_k^TR_{du,k}\Delta u_k,
\qquad \Delta u_0=u_0-u_{\mathrm{prev}}.
$$

**Cost weights**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.Rdu` | Quadratic weight on changes in the control input | $n_u\times n_u\times N$ | $k=0,\ldots,N-1$ |

**Associated signals**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.du` | Changes in the control input per sample | $n_u\times N$ | $k=0,\ldots,N-1$ |
| `mpc.xu` | Preceding control input at each predicted state stage | $n_u\times N$ | $k=1,\ldots,N$ |

`mpc.du` stores the input changes. `mpc.xu` stores the preceding control
at each predicted state stage. At $k=0$, CHRONOS uses `u_prev` directly.

(mpc-cost-terminal-state)=
### Set a terminal state objective

The terminal cost is an important part of designing for closed-loop
stability. It lets the controller account for how costly it will be to bring
the final predicted state to the target after the horizon ends. If the system
oscillates with a short horizon, a well-chosen terminal cost can help reduce
those oscillations.

Choose one of two ways to set up the terminal cost.

- **Supply your own terminal weight.** Pass `P` as a single symmetric
  positive definite `nx`-by-`nx` matrix to
  [`init_mpc_ter_ingredients`](../api/init/init_mpc_ter_ingredients.md).
- **Compute the terminal weight with DLQR.** Supply one state-weight matrix
  `Qx` and one control-weight matrix `Ru` to
  [`init_mpc_ter_ingredients_dlqr`](../api/init/init_mpc_ter_ingredients_dlqr.md). This initializer computes `P` using the supplied weights and the `A` and `B` matrices stored during dynamics initialization. Choose a representative operating point for those matrices.

At solve time, pass an `nx`-by-1 full-state reference as `xN_ref_in`. If the
tracking output is the full state in the same order, set `xN_ref_is_y=true`
in the initializer and pass `xN_ref_in=[]` to reuse the terminal tracking
reference. An explicit `xN_ref_in` takes precedence.

$$
J_{\mathrm{terminal}}=(x_{N,\mathrm{ref}}-x_N)^TP
(x_{N,\mathrm{ref}}-x_N).
$$

**Cost weights**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.P` | Quadratic weight on the terminal state error | $n_x\times n_x$ | $k=N$ |

**Associated signals**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.xN_ref` | Full-state reference for the terminal state | $n_x\times 1$ | $k=N$ |

(mpc-cost-custom)=
### Customize the controller's behavior

If the standard MPC costs do not capture an objective you want the
controller to pursue, CHRONOS lets you define a custom signal $z$ that describes your objective and add a cost on it. You can combine the state, current and previous control inputs, and a known input in $z$ to form the signal model. It is defined at each stage as

$$
z_k=C_{z,k}x_k+D_{z,k}u_k+D_{xu,z,k}u_{k-1}+D_{dz,k}d_{z,k}
$$

You can add a quadratic weight `Qz`, a linear weight `qz`, or both to penalize the signal $z$. Supply the signal coefficients and weights to [`init_mpc_Custom_cost`](../api/init/init_mpc_Custom_cost.md).

Supply at least one active coefficient among `Cz`, `Dz`, and `Dxuz` and at least one active weight. If the signal uses a known profile $d_z$, pass it as `dz_in` when solving with `mpc_solve`. Note that this profile is independent of `d_in` from the prediction model known inputs. 

$$
J_{\mathrm{custom}}=\sum_{k\in\mathcal K_z}
\left(\frac12z_k^TQ_{z,k}z_k+q_{z,k}^Tz_k\right).
$$

**Signal model coefficients**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.Cz` | Coefficient on the state | $n_z\times n_x\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.Dz` | Coefficient on the current control input | $n_z\times n_u\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.Dxuz` | Coefficient on the previous control input | $n_z\times n_u\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.Ddz` | Coefficient on the known input | $n_z\times n_{dz}\times (N-1)$ | $k=1,\ldots,N-1$ |

**Cost weights**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.Qz` | Quadratic weight on the custom signal | $n_z\times n_z\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.qz` | Linear weight on the custom signal | $n_z\times (N-1)$ | $k=1,\ldots,N-1$ |

**Associated signals**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.z` | Computed custom cost signal | $n_z\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.dz` | Known input supplied through `dz_in` | $n_{dz}\times N$ | $k=0,\ldots,N-1$ |
| `mpc.xu` | Preceding control input at each predicted state stage | $n_u\times N$ | $k=1,\ldots,N$ |
| `mpc.du` | Changes in the control input per sample | $n_u\times N$ | $k=0,\ldots,N-1$ |

`mpc.z` stores the computed custom signal. An active `Ddz` enables
`mpc.dz` for the known input.
An active `Dxuz` enables the same `mpc.xu` and `mpc.du`
fields used for control-rate terms.

If the signal model uses an active `mpc.Dz` field, CHRONOS will automatically populate the data required to use the custom cost at $k=0$. Similarly, if the signal model has rows that depend exclusively on the states, CHRONOS will automatically populate the data required to use the custom cost at $k=N$.

#### Example: Penalize the control rate as a custom cost

Let's consider adding a custom cost to mimic the behavior of the control-rate cost. Our $z$ signal would be the control-rate change:

$$z_k = u_k-u_{k-1}$$

Since the custom cost model exposes both the current control action and the previous-stage action, we can set up the model with the following coefficients.

| $C_{z,k}$ | $D_{z,k}$ | $D_{xu,z,k}$ | $D_{dz,k}$ |
| --- | --- | --- | --- | 
| [] | 1 | -1 | [] | 

#### Example: Penalize the control rate and the tracking error as a custom cost

Let's now consider recreating both the control-rate cost and the tracking-error cost using the custom cost formulation. Our custom signal would have two rows, one for each cost:

$$z_k = 
\begin{bmatrix}
u_k-u_{k-1}\\
r_k-y_k
\end{bmatrix}
$$

For simplicity, let's assume that the tracking signal is the full state vector, $y=x$. We can set up the model with the following coefficients.

| $C_{z,k}$ | $D_{z,k}$ | $D_{xu,z,k}$ | $D_{dz,k}$ |
| --- | --- | --- | --- | 
| $\begin{bmatrix} 0\\-I\end{bmatrix}$ | $\begin{bmatrix} 1\\0\end{bmatrix}$ | $\begin{bmatrix} -1\\0\end{bmatrix}$  | $\begin{bmatrix} 0\\I\end{bmatrix}$ | 

In this example, we need to pass the tracking reference to `mpc_solve` through the `dz_in` input field.

#### Explore more custom-cost examples

Explore more uses of custom costs in these repository examples.

- [Virtual-reference MPC for the twin-rotor system](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/TRMS/TRMS_refMPC/TRMS_refMPC_init.m)
- [First Step tutorial for the two-tank system](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/tutorials/01%20-%20First%20Step/beginning_tutorial.mlx)

(defining-mpc-constraints)=
## Constraints

Decide which limits the controller must respect and where it can trade a
violation against its other objectives.

| Constraint type | When to use it | Enforcement |
| --- | --- | --- |
| {ref}`State <mpc-constraint-state>` | Keep predicted states within an operating range. | Soft |
| {ref}`Control input <mpc-constraint-control>` | Keep actuator inputs within their operating limits. | Hard |
| {ref}`Control rate <mpc-constraint-control-rate>` | Limit how much the control input can change per sample. | Hard |
| {ref}`Output <mpc-constraint-output>` | Bound a quantity defined by the output model. | Soft |
| {ref}`Custom signal <mpc-constraint-custom>` | Bound a signal you define from the state, current or previous input, or known signals. | Soft |


For soft constraints, CHRONOS adds slack variables to prevent the optimization problem from becoming infeasible when working near the limits. The slack values can be penalized individually for each lower and upper bound constraint and for each stage of the prediction horizon.

The slack penalty is a linear cost term, $J_{\mathrm{slack}}=q_vv$, where $q_v$ is the slack cost and $v$ is the slack value. The slack cost determines how the constraint behaves near its boundaries.

* A **small** slack cost $q_v$ leads to a relaxed constraint that can be violated.
* A **large** slack cost $q_v$ leads to a constraint that behaves similarly to a hard constraint, although it can still be violated if the solver has no room to maneuver.

If you do not set a value for the slack cost on soft constraints, CHRONOS automatically selects the value stored in the field `mpc.qv`.


(mpc-constraint-state)=
### Limit the predicted states

To keep states within an operating range, supply
`x_min` and `x_max` to
[`init_mpc_state_cnstr`](../api/init/init_mpc_state_cnstr.md).
Use `nx` rows for a component column or an `nx`-by-`L` profile. Use
`qv_min` and `qv_max` with the same layout to set the cost of relaxing each
side. The first column applies to the first predicted state $x_1$. The
bounds continue through the terminal state $x_N$, with the terminal
penalties stored separately.

**`mpc.x_cnstr`**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.x_cnstr.min` | Lower bound | $n_x\times N$ | $k=1,\ldots,N$ |
| `mpc.x_cnstr.max` | Upper bound | $n_x\times N$ | $k=1,\ldots,N$ |
| `mpc.x_cnstr.qv_min` | Penalty for relaxing the lower bound | $n_x\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.x_cnstr.qv_max` | Penalty for relaxing the upper bound | $n_x\times (N-1)$ | $k=1,\ldots,N-1$ |

**Associated signals**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.x` | Predicted states to which the bounds apply | $n_x\times N$ | $k=1,\ldots,N$ |

(mpc-constraint-control)=
### Limit the control input

To keep each actuator within its operating limits, supply `u_min` and
`u_max` to [`init_mpc_control_cnstr`](../api/init/init_mpc_control_cnstr.md).
Use `nu` rows for a component column or an `nu`-by-`L` profile. The first
column bounds the next control action $u_0$. For example, use a lower bound
of zero for an actuator that only operates in one direction.

$$
u_{\min,k}\le u_k\le u_{\max,k},\qquad k=0,\ldots,N-1.
$$

**`mpc.u_cnstr`**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.u_cnstr.min` | Lower bound | $n_u\times N$ | $k=0,\ldots,N-1$ |
| `mpc.u_cnstr.max` | Upper bound | $n_u\times N$ | $k=0,\ldots,N-1$ |

**Associated signals**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.u` | Control inputs to which the bounds apply | $n_u\times N$ | $k=0,\ldots,N-1$ |

(mpc-constraint-control-rate)=
### Limit changes in control

To cap how far an actuator can move in one sample, supply `du_min` and
`du_max` to
[`init_mpc_control_rate_cnstr`](../api/init/init_mpc_control_rate_cnstr.md).
Use `nu` rows for a component column or an `nu`-by-`L` profile. Set the
bounds in units of input increment per sample. For example, `-0.1` and
`0.1` allow each input to change by up to 0.1 in either direction per sample.

The first column bounds the move from the last applied input to $u_0$.
Pass that last applied input as `u_prev` when calling
[`mpc_solve`](../api/online/mpc_solve.md).

$$
\begin{aligned}
\Delta u_0 &= u_0-u_{\mathrm{prev}}, \\
\Delta u_k &= u_k-xu_k, && k=1,\ldots,N-1, \\
\Delta u_{\min,k} &\le \Delta u_k \le \Delta u_{\max,k}, && k=0,\ldots,N-1.
\end{aligned}
$$

**`mpc.du_cnstr`**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.du_cnstr.min` | Lower bound | $n_u\times N$ | $k=0,\ldots,N-1$ |
| `mpc.du_cnstr.max` | Upper bound | $n_u\times N$ | $k=0,\ldots,N-1$ |

**Associated signals**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.du` | Changes in the control input per sample | $n_u\times N$ | $k=0,\ldots,N-1$ |
| `mpc.xu` | Preceding control input at each predicted state stage | $n_u\times N$ | $k=1,\ldots,N$ |

Rate bounds enable `mpc.du` for the input changes and `mpc.xu` for the
preceding-input term at predicted state stages. These fields are shared
with the rate cost when it is active.

(mpc-constraint-output)=
### Bound the predicted output

To keep a quantity defined by your output model within a range, supply
`y_min` and `y_max` to
[`init_mpc_output_cnstr`](../api/init/init_mpc_output_cnstr.md).
Use `ny` rows for a component column or an `ny`-by-`L` profile. Use
`qv_min` and `qv_max` to set the relaxation penalties. 

CHRONOS bounds the output signal at $k=0$ for rows with a nonzero control coefficient in the first `D` page. At $k=N$, only rows that depend on the state alone at that stage are bounded.

If the output model uses a known input, supply it as `d_in` when calling `mpc_solve`. The dynamics and output model share this input.


**`mpc.y_cnstr`**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.y_cnstr.min` | Lower bound | $n_y\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.y_cnstr.max` | Upper bound | $n_y\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.y_cnstr.qv_min` | Penalty for relaxing the lower bound | $n_y\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.y_cnstr.qv_max` | Penalty for relaxing the upper bound | $n_y\times (N-1)$ | $k=1,\ldots,N-1$ |

**Associated signals**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.y` | Output computed from the configured output model | $n_y\times (N-1)$ | $k=1,\ldots,N-1$ |

CHRONOS computes `mpc.y` from the configured output model. Tracking uses
the same signal when a tracking cost is active.

(mpc-constraint-custom)=
### Constrain a custom quantity

If the standard bounds do not describe the limit you need, define a
separate signal $h_k$. It can combine the state, current input, preceding
input, and a dedicated known input $d_{h,k}$.

Define the signal and its bounds with
[`init_mpc_custom_cnstr`](../api/init/init_mpc_custom_cnstr.md).
Supply `Ch`, `Dh`, `Dxuh`, and `Ddh` as constant matrices or profiles
with `L` pages. Supply at least one nonzero coefficient among `Ch`, `Dh`,
and `Dxuh`. Use `[]` for unused coefficients. Set `h_min` and `h_max` and
optionally `qv_min` and `qv_max`, with `nh` rows in each component column
or `nh`-by-`L` profile.

CHRONOS bounds the signal $h_k$ at $k=0$ for rows with a nonzero control coefficient in the first `Dh` page. Meanwhile, at $k=N$, only rows that depend on the state alone at that stage are bounded.

If you set a nonzero `Ddh`, pass the known input as `dh_in` when building
and when supplying its values for a solve. Use an `ndh`-by-1 column or an
`ndh`-by-`L` profile starting at $k=0$. This input is separate from `d_in`
for dynamics and outputs and `dz_in` for custom costs.

$$
h_k=C_{h,k}x_k+D_{h,k}u_k+D_{xu,h,k}u_{k-1}+D_{dh,k}d_{h,k}
$$

**Signal model coefficients**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.Ch` | Coefficient on the state | $n_h\times n_x\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.Dh` | Coefficient on the current control input | $n_h\times n_u\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.Dxuh` | Coefficient on the previous control input | $n_h\times n_u\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.Ddh` | Coefficient on the known input | $n_h\times n_{dh}\times (N-1)$ | $k=1,\ldots,N-1$ |

**`mpc.h_cnstr`**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.h_cnstr.min` | Lower bound | $n_h\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.h_cnstr.max` | Upper bound | $n_h\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.h_cnstr.qv_min` | Penalty for relaxing the lower bound | $n_h\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.h_cnstr.qv_max` | Penalty for relaxing the upper bound | $n_h\times (N-1)$ | $k=1,\ldots,N-1$ |

**Associated signals**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.h` | Computed custom constraint signal | $n_h\times (N-1)$ | $k=1,\ldots,N-1$ |
| `mpc.dh` | Known input supplied through `dh_in` | $n_{dh}\times N$ | $k=0,\ldots,N-1$ |
| `mpc.xu` | Preceding control input at each predicted state stage | $n_u\times N$ | $k=1,\ldots,N$ |
| `mpc.du` | Changes in the control input per sample | $n_u\times N$ | $k=0,\ldots,N-1$ |

CHRONOS computes `mpc.h` from the supplied coefficients. An active `Ddh`
enables `mpc.dh`. An active `Dxuh` enables the same `mpc.xu` and `mpc.du`
fields used for control-rate terms.

#### Example: Keep two states close to each other

Let's assume a two-state system in which the two states must remain close. The signal to constrain is the difference between them:

$$
h_k = x_{1,k}-x_{2,k}.
$$

The signal model can then be defined with the following choice of coefficients:

| $C_{h,k}$ | $D_{h,k}$ | $D_{xu,h,k}$ | $D_{dh,k}$ |
| --- | --- | --- | --- | 
| [1 -1] | [] | [] | [] | 

Now, we can set `h_max` and `h_min` to the maximum and minimum deviation allowed between the states. The slack costs `qv_min` and `qv_max` determine how strongly this constraint is enforced.

#### Example: Keep the input near a supplied target

To keep the control input close to a supplied target, define the constrained signal as the difference between the input and that target.

$$
h_k=u_k-u_{\mathrm{ref},k},
$$

where $u_{\mathrm{ref},k}$ is the reference control command, which is treated as a known input. Define the signal model with the following coefficients.

| $C_{h,k}$ | $D_{h,k}$ | $D_{xu,h,k}$ | $D_{dh,k}$ |
| --- | --- | --- | --- |
| [] | $I$ | [] | $-I$ |

Set `h_min` and `h_max` to the allowed deviation from the target, such as
`-0.2` and `0.2`. Use `qv_min` and `qv_max` to set the cost of relaxing
those bounds. Pass the input target as `dh_in`. The bound includes the
next control action because the signal depends on $u_0$.

#### Example: Keep a minimum distance from a reference point

For a two-state position model, suppose the controller should keep the
position $(x_k,y_k)$ farther than $d_{\mathrm{ref}}\ge0$ from a fixed point
$(x_0,y_0)$. Squaring the distance gives the equivalent condition

$$
d_k=\sqrt{(x_k-x_0)^2+(y_k-y_0)^2},
\qquad d_k>d_{\mathrm{ref}}
\iff d_k^2>d_{\mathrm{ref}}^2.
$$

Use the squared-distance factorization with scheduled positions
$(\hat x_k,\hat y_k)$, taken from the previous MPC predictions and held
fixed during the solve. The resulting custom signal approximates $d_k^2$.
See the {ref}`distance factorization example <squared-distance-to-a-reference-point>`.

$$
h_k=
\begin{bmatrix}\hat x_k-2x_0 & \hat y_k-2y_0\end{bmatrix}
\begin{bmatrix}x_k\\y_k\end{bmatrix}
+x_0^2+y_0^2.
$$

Define the signal model with the following choice of coefficients.

| $C_{h,k}$ | $D_{h,k}$ | $D_{xu,h,k}$ | $D_{dh,k}$ |
| --- | --- | --- | --- |
| $\begin{bmatrix}\hat x_k-2x_0 & \hat y_k-2y_0\end{bmatrix}$ | [] | [] | $x_0^2+y_0^2$ |

Pass `dh_in = 1` to supply the constant offset. To request a margin on the
approximated squared distance, choose $\epsilon>0$ in squared-distance units
and set `h_min` to $d_{\mathrm{ref}}^2+\epsilon$. Leave `h_max` empty and
use `qv_min` to set the cost of relaxing the lower bound. This is a soft
bound on a scheduled approximation: neither the margin nor a large slack
penalty guarantees strict separation in the original nonlinear distance.

#### Explore more custom-constraint examples

Explore custom constraints in these repository examples.

- [First Step tutorial for the two-tank system](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/tutorials/01%20-%20First%20Step/beginning_tutorial.mlx)

(defining-mpc-build)=
## Building the controller

After defining the prediction model, costs, and constraints, call
[`build_chronos_mpc`](../api/init/build_chronos_mpc.md) to prepare the
controller for its first solve. The builder allocates the solver workspace,
prepares the configured costs and constraints, resolves default slack
penalties, and initializes the predicted trajectories.

### Start from the current operating point

Supply the current measured or estimated state as `x_prev` and the input
applied at the preceding sample as `u_prev`. Pass the configured `mpc`
structure and retain the returned structure for subsequent solver calls.

```matlab
mpc = build_chronos_mpc(mpc, x_prev, u_prev);
```

**Required inputs**

| Input | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `x_prev` | Current measured or estimated state | $n_x\times 1$ | $k=0$ |
| `u_prev` | Last applied control input | $n_u\times 1$ | Previous sample |

The builder uses these values to start a prediction rollout. Choose hard
control and control-rate bounds that are mutually compatible with the
rollout and the previously applied input.

### Supply known input forecasts

If a known input affects the dynamics or output model, pass its values as
`d_in`. If the custom constraint uses `Ddh`, pass its dedicated known
input as `dh_in`. Supply these values only for configured model terms.

Use an `nd`-by-1 column for a constant `d_in` or an `nd`-by-`L` forecast.
For `dh_in`, use `ndh` rows with the same column convention. Here `L` is
the number of supplied columns. The first column applies at $k=0$.
CHRONOS repeats the final column of a short forecast and ignores columns
beyond the horizon.

```matlab
mpc = build_chronos_mpc(mpc, x_prev, u_prev, d_in, dh_in);
```

**Known-input storage**

| Field | Meaning | Size | Physical stages |
| --- | --- | --- | --- |
| `mpc.d` | Known input shared by the dynamics and output model | $n_d\times N$ | $k=0,\ldots,N-1$ |
| `mpc.dh` | Dedicated known input for custom constraints | $n_{dh}\times N$ | $k=0,\ldots,N-1$ |

These fields are active when their corresponding known-input terms are
configured. Newly initialized storage starts at zero. Passing `[]` retains
the stored values.

### Guide the initial prediction

To guide the initial rollout toward a state target, supply `x_ref` as an
`nx`-by-1 column when a rollout gain `mpc.K` is available. For example,
[`init_mpc_ter_ingredients_dlqr`](../api/init/init_mpc_ter_ingredients_dlqr.md)
computes this gain when setting up terminal ingredients. At each rollout
stage, the builder uses

$$
u_{\mathrm{raw},k}=K(x_{\mathrm{ref}}-x_k).
$$

Without both the gain and the reference, the rollout starts each stage
from the preceding input. The builder adjusts the trial input using the
configured control-rate and control bounds.

The optional arguments default to `[]`. Omit trailing optional arguments
or use `[]` to skip an argument before a later one. For example, to supply
only the rollout reference, call

```matlab
mpc = build_chronos_mpc(mpc, x_prev, u_prev, [], [], x_ref);
```

## Next step

Keep the built `mpc` structure for the online control loop. At each sampling
instant, update any configured data that have changed, solve using the current
state, last applied input, references, and known-input forecasts, and apply
the first control action.

Continue with [Solving and updating the MPC](solving-and-updating-the-mpc.md) for the update
and solve calls. Changes to the controller's structure require redefining
and rebuilding it.
