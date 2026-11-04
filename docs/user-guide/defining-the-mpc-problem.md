# Defining the MPC problem

This chapter describes the model, costs, and constraints supported by CHRONOS
and how to configure them through the initialization API. The
[API reference](../api/init/index.md) provides the full argument descriptions.

```{contents}
:local:
:depth: 2
```

## MPC formulation

For a supplied LPV/LTV model sequence, CHRONOS solves the following MPC problem

$$
\begin{aligned}
\min_{x,u}\; J ={}&
\frac12\sum_{k\in\mathcal K_y}(r_k-y_k)^TQ_{e,k}(r_k-y_k)\\
&+\sum_{k=0}^{N-1}\left(\frac12u_k^TR_{u,k}u_k+r_{u,k}^Tu_k\right)\\
&+\frac12\sum_{k=0}^{N-1}\Delta u_k^TR_{du,k}\Delta u_k\\
&+\sum_{k\in\mathcal K_z}\left(\frac12z_k^TQ_{z,k}z_k+q_{z,k}^Tz_k\right)\\
&+(x_{N,\mathrm{ref}}-x_N)^TP(x_{N,\mathrm{ref}}-x_N).
\end{aligned}
$$

where prediction dynamics and desired operating bounds are

$$
\begin{aligned}
x_{k+1}&=A_kx_k+B_ku_k+B_{d,k}d_k, && k=0,\ldots,N-1,\\
u_{\min,k}&\le u_k\le u_{\max,k}, && k=0,\ldots,N-1,\\
\Delta u_{\min,k}&\le\Delta u_k\le\Delta u_{\max,k}, && k=0,\ldots,N-1,\\
x_{\min,k}&\le x_k\le x_{\max,k}, && k=1,\ldots,N,\\
y_{\min,k}&\le y_k\le y_{\max,k}, && k\in\mathcal K_y,\\
h_{\min,k}&\le h_k\le h_{\max,k}, && k\in\mathcal K_h.
\end{aligned}
$$

| Symbol | Meaning |
| --- | --- |
| $N$ | Number of prediction steps |
| $x_k\in\mathbb R^{n_x}$ | Predicted state|
| $u_k\in\mathbb R^{n_u}$ | Predicted control input |
| $d_k\in\mathbb R^{n_d}$ | Known input shared by dynamics and outputs |
| $y_k\in\mathbb R^{n_y}$ | Output signal and its tracking reference |
| $r_k\in\mathbb R^{n_y}$ | Tracking reference |
| $z_k\in\mathbb R^{n_z}$ | Custom cost signal |
| $h_k\in\mathbb R^{n_h}$ | Custom constraint signal |

The model and signal coefficients are fixed during a solve, however, they can change between prediction stages. Allowing to account for nonlinear behaviour in the form of time varying models. With valid convex weights, the problem has linear dynamics, affine constraints, and a convex quadratic/linear objective.

### CHRONOS MPC definition workflow

**Create the MPC structure → Configure the ingredients → Build the controller**

Always start by initializing a CHRONOS MPC structure:

```matlab
mpc = init_chronos_mpc(N);
```

[`init_chronos_mpc(N)`](../api/init/init_chronos_mpc.md) selects the horizon and creates the
structure to be populated by subsequent calls. The supported horizon is an
integer $N\ge1$; omitting the argument selects $N=10$. For a model sampled
every $T_s$ seconds, $N$ prediction steps cover $NT_s$ seconds. The model
matrices must already describe that discrete-time sampling interval.

Define dynamics before ingredients that use their dimensions. Define the
output before tracking costs or output constraints. Finally,
{ref}`build_chronos_mpc <defining-mpc-build>` prepares the configured
problem for solving; its inputs are explained at the end of this chapter.

The following sections describe the CHRONOS ingredients and the data supplied
to their initialization functions.

## Horizon stages and time-varying data

CHRONOS supports variation across the prediction horizon in model matrices,
references, stage cost weights, bounds, and slack penalties. This allows a
controller to account for changing dynamics, upcoming operating limits, or
objectives that change further into the prediction. The terminal weight `P`
is a single matrix rather than a horizon profile.

### Supplying time-varying data

A vector uses its **second dimension** for time: each column contains one
stage's values. A matrix uses its **third dimension**: each page contains
one stage's coefficients.

| Data | Constant input | Time-varying input |
| --- | --- | --- |
| Vector with `n` components | `n`-by-1 | `n`-by-`L` |
| Matrix with `nr` rows and `nc` columns | `nr`-by-`nc` | `nr`-by-`nc`-by-`L` |

Here `L` is the number of supplied stages. A single column or matrix is reused
across the horizon. Short profiles repeat their final column or page; excess
stages are ignored. For example, a state-bound profile with three columns and
`N=5` uses its third column for prediction steps 3, 4, and 5.

Bounds and slack penalties also accept scalars, expanded to every component.
For a multi-component profile, supply all component rows: a `1`-by-`L` row
does not broadcast across multiple components.

### Prediction stages

The decision variables depend on the prediction stage. At the current stage,
the measured state is known and the next control action is selected. Further
ahead, both states and controls are predicted; the horizon ends with a final
state and no additional control action.

| Stage | Interpretation |
| --- | --- |
| $k=0$ | Current state `x_prev` and previous input `u_prev` are known; $u_0$ is optimized. |
| $k=1,\ldots,N-1$ | Predicted states and controls are optimized. |
| $k=N$ | Final predicted state; there is no $u_N$ decision. |

The main trajectory fields reflect these stage ranges:

| Signal | Field | Physical stages |
| --- | --- | --- |
| State | `mpc.x` | $k=1,\ldots,N$ |
| Control | `mpc.u` | $k=0,\ldots,N-1$ |
| Control change | `mpc.du` | $k=0,\ldots,N-1$ |
| Known inputs | `mpc.d`, `mpc.dz`, `mpc.dh` | $k=0,\ldots,N-1$ |
| Interior output and custom signals | `mpc.y`, `mpc.z`, `mpc.h` | $k=1,\ldots,N-1$ |

For example, `mpc.x(:,1)` is the predicted state at $k=1$, whereas
`mpc.u(:,1)` is the control at $k=0$. States and controls are solver results;
references, known inputs, weights, and bounds are data you supply.

For defining profiles, the first dynamics page and first control-related
column/page apply at $k=0$; the first state-bound column applies at $k=1$.
Output/custom profiles use columns or pages 1 through `N` for prediction
steps 1 through `N`. Eligible current-stage terms reuse the first supplied
column/page, so these are not independent `N+1`-stage profiles. CHRONOS
handles endpoint row selection automatically. With `N=1`, there are only the
current and terminal stages.

### Updating between solves

A horizon profile describes conditions anticipated **within one solve**.
The [online update workflow](online-updates.md) replaces numerical data **between
solves**, for example when a new model schedule or forecast is available.
Updates preserve the built dimensions and enabled ingredients. Use these
functions rather than editing coefficient fields directly, so dependent
cost data is refreshed.





## Prediction model

### Dynamics

$$
x_{k+1}=A_kx_k+B_ku_k+B_{d,k}d_k.
$$

| Coefficient | Per-stage size | Field |
| --- | --- | --- |
| $A_k$ | $n_x\times n_x$ | `mpc.A` |
| $B_k$ | $n_x\times n_u$ | `mpc.B` |
| $B_{d,k}$ | $n_x\times n_d$ | `mpc.Bd` when active |

[`init_mpc_dynamics`](../api/init/init_mpc_dynamics.md) establishes the state
and control dimensions and initializes the default output $y=x$. The known
input is supplied separately as `d_in` when building or solving. `A`, `B`,
and active `Bd` may each be constant or staged.

### Outputs

$$
y_k=C_kx_k+D_ku_k+D_{d,k}d_k.
$$

The output selects or combines quantities for tracking and output constraints.

| Coefficient | Per-stage size | Field |
| --- | --- | --- |
| $C_k$ | $n_y\times n_x$ | `mpc.C` |
| $D_k$ | $n_y\times n_u$ | `mpc.D` |
| $D_{d,k}$ | $n_y\times n_d$ | `mpc.Dd` |

[`init_mpc_output`](../api/init/init_mpc_output.md) replaces the default
output. The disturbance dimension is shared with the dynamics. Empty or
all-zero optional coefficient matrices are inactive; at least `C` or `D`
must contain an active term to define an output.

At the current stage, only output rows that depend on the current control
contribute. At the terminal stage, only state-dependent rows without control
or disturbance dependence contribute. CHRONOS selects these rows automatically.

## Cost terms

Stage quadratic weights are assumed symmetric positive semidefinite; their
definiteness is not checked by the initializer. Choose weights that produce
a well-posed optimization problem. The following costs are optional and may
be combined. Slack penalties are discussed with the constraints below.

### Reference tracking

$$
J_{\mathrm{tracking}}=\frac12\sum_{k\in\mathcal K_y}
(r_k-y_k)^TQ_{e,k}(r_k-y_k).
$$

The weight `Qe` has per-stage size `ny`-by-`ny` and can vary along the
horizon. Increasing its later pages gives tracking errors further ahead more
weight. The output model determines the eligible initial and terminal terms.

Supply the reference as `r_in` to `mpc_solve`, using `ny`-by-1 or
`ny`-by-`L` data. The error is $r_k-y_k$. The corresponding fields are
`mpc.Qe`, `mpc.r`, and `mpc.err`.

Reference: [`init_mpc_Tracking_cost`](../api/init/init_mpc_Tracking_cost.md).

### Control effort

$$
J_{\mathrm{effort}}=\sum_{k=0}^{N-1}
\left(\frac12u_k^TR_{u,k}u_k+r_{u,k}^Tu_k\right).
$$

The quadratic term penalizes input magnitude; the linear term allows a
direction-dependent cost. `Ru` is `nu`-by-`nu`-by-`L`, and `ru` is
`nu`-by-`L`. The corresponding fields are `mpc.Ru` and `mpc.ru`.
There is no terminal control-effort term.

Either weight can be omitted independently. Reference:
[`init_mpc_Control_cost`](../api/init/init_mpc_Control_cost.md).

### Control rate

$$
J_{\mathrm{rate}}=\frac12\sum_{k=0}^{N-1}
\Delta u_k^TR_{du,k}\Delta u_k,
\qquad \Delta u_0=u_0-u_{\mathrm{prev}}.
$$

`Rdu` is `nu`-by-`nu`-by-`L`, stored as `mpc.Rdu`.
This cost penalizes changes between successive inputs, rather than input
magnitude. It enables the previous-control signal `mpc.xu`; `mpc.du` stores
the resulting input differences. These differences are not divided by the
sampling period.

Reference: [`init_mpc_ControlRate_cost`](../api/init/init_mpc_ControlRate_cost.md).

### Terminal cost

$$
J_{\mathrm{terminal}}=(x_{N,\mathrm{ref}}-x_N)^TP
(x_{N,\mathrm{ref}}-x_N).
$$

This term has no factor of $1/2$. `mpc.P` is one `nx`-by-`nx` matrix,
with symmetric positive definiteness assumed. `mpc.xN_ref` is an
`nx`-by-1 reference. The terminal cost is additional to any terminal output
tracking cost.

Supply `P` directly with
[`init_mpc_ter_ingredients`](../api/init/init_mpc_ter_ingredients.md), or use
[`init_mpc_ter_ingredients_dlqr`](../api/init/init_mpc_ter_ingredients_dlqr.md)
to compute `P` and an initial-rollout gain from the final horizon transition.
The DLQR inputs `Qx` and `Ru` are single matrices, not horizon profiles.

At solve time, pass a full-state reference as `xN_ref_in`. Alternatively,
the initializer's optional `xN_ref_is_y=true` allows reuse of the terminal
tracking reference when `xN_ref_in=[]`, but only when the output is the
full state in the same order. An explicit `xN_ref_in` takes precedence.

### Custom costs

$$
\begin{aligned}
z_k&=C_{z,k}x_k+D_{z,k}u_k+D_{xu,z,k}xu_k+D_{dz,k}d_{z,k},\\
J_{\mathrm{custom}}&=\sum_{k\in\mathcal K_z}
\left(\frac12z_k^TQ_{z,k}z_k+q_{z,k}^Tz_k\right).
\end{aligned}
$$

| Input | Per-stage size | Field |
| --- | --- | --- |
| `Cz` | `nz`-by-`nx` | `mpc.Cz` |
| `Dz`, `Dxuz` | `nz`-by-`nu` | `mpc.Dz`, `mpc.Dxuz` |
| `Ddz` | `nz`-by-`ndz` | `mpc.Ddz` |
| `Qz` | `nz`-by-`nz` | `mpc.Qz` |
| `qz` | `nz`-by-1 | `mpc.qz` |

Coefficients and quadratic weights use matrix pages; linear weights use
columns. `mpc.z` stores the derived signal. The previous-control signal
$xu_k=u_{k-1}$ uses the last applied input at $k=0$. `dz_in` supplies the
dedicated known input `mpc.dz`, independently of `d_in`.

For example, defining $z=u-d_z$ allows a cost on deviation from a known
control-reference profile. At least one of `Cz`, `Dz`, or `Dxuz` and one
cost weight must be active. Current-stage terms require current-control
dependence; terminal terms must depend on the state only. CHRONOS handles
this selection from the supplied coefficients.

Reference: [`init_mpc_Custom_cost`](../api/init/init_mpc_Custom_cost.md).

(defining-mpc-constraints)=
## Constraints and slack penalties

### Hard and soft bounds

Control and control-rate bounds are hard. State, output, and custom bounds
use separate feasibility slacks for each enabled side:

$$
q_{\min,k}-v^-_k\le q_k\le q_{\max,k}+v^+_k,
\qquad
J_{\mathrm{slack},k}=(q^-_{v,k})^Tv^-_k+(q^+_{v,k})^Tv^+_k.
$$

Here $q$ is the bounded signal and $v^-,v^+\ge0$ are its lower- and
upper-bound relaxations. The penalty inputs are `qv_min` and `qv_max`.
Larger positive penalties make
violations more costly; they do not convert a soft bound into a hard one.
These penalties are added to the controller cost shown in the opening
formulation.

Limits and penalties accept a scalar, a constant component column, or a
component-by-stage array. `[]` omits a bound; zero is an active bound.
Omitted or all-zero penalties request defaults resolved during building
from available tracking weights and `mpc.qv`. Explicit penalties are useful
when the application requires a particular tradeoff.

### State constraints

State bounds apply to $x_1,\ldots,x_N$, not the supplied current state.
Supply `x_min` and `x_max` as `nx`-component columns or `nx`-by-`L` profiles.
The first column bounds the first predicted state. The `qv_min` and `qv_max`
inputs set the penalties for relaxing each bound and can follow the same
horizon profile.

Reference: [`init_mpc_state_cnstr`](../api/init/init_mpc_state_cnstr.md).

### Control constraints

Input bounds apply to $u_0,\ldots,u_{N-1}$. Supply `u_min` and `u_max`
as `nu`-component columns or `nu`-by-`L` profiles. The first column limits
the next control action. No terminal bound or feasibility slack is added.

Reference: [`init_mpc_control_cnstr`](../api/init/init_mpc_control_cnstr.md).

### Control-rate constraints

These bounds limit $\Delta u_k=u_k-u_{k-1}$ for $k=0,\ldots,N-1$,
using the last applied input at $k=0$. Supply `du_min` and `du_max` as
`nu`-component columns or `nu`-by-`L` profiles. They have no feasibility
relaxation.

The limits are input increments per sample, not continuous-time derivatives.
Reference: [`init_mpc_control_rate_cnstr`](../api/init/init_mpc_control_rate_cnstr.md).

### Output constraints

Output bounds act on the previously defined $y_k$. Supply `y_min` and
`y_max` with `ny` rows, and use `qv_min` and `qv_max` to set their relaxation
penalties. Bounds and penalties may vary across columns, following the output
profile convention. CHRONOS applies the output model's endpoint selection.

For example, lower the upper limit in later columns to anticipate a tighter
operating region. A decreasing penalty profile allows more flexibility in
predictions further ahead, while prioritizing near-term constraint satisfaction.

Reference: [`init_mpc_output_cnstr`](../api/init/init_mpc_output_cnstr.md).

### Custom constraints

Define a signal independently of the tracking output:

$$
h_k=C_{h,k}x_k+D_{h,k}u_k+D_{xu,h,k}xu_k+D_{dh,k}d_{h,k}.
$$

| Input | Per-stage size | Field |
| --- | --- | --- |
| `Ch` | `nh`-by-`nx` | `mpc.Ch` |
| `Dh`, `Dxuh` | `nh`-by-`nu` | `mpc.Dh`, `mpc.Dxuh` |
| `Ddh` | `nh`-by-`ndh` | `mpc.Ddh` |
| `h_min`, `h_max` | `nh`-by-1 | `mpc.h_cnstr.min`, `.max` |
| `qv_min`, `qv_max` | `nh`-by-1 | `mpc.h_cnstr.qv_min`, `.qv_max` |

Matrix coefficients use pages and bounds/penalties use columns. `mpc.h`
stores the derived signal. The known input `dh_in` supplies `mpc.dh`
independently of the dynamics and custom-cost inputs.

For example, choose `Ch` to sum selected states and apply an upper bound to
that sum. At least one of `Ch`, `Dh`, or `Dxuh` and one bound must be active.
Current-stage rows require current-control dependence; terminal rows must
depend on the state only. The API handles this selection automatically.

Reference: [`init_mpc_custom_cnstr`](../api/init/init_mpc_custom_cnstr.md).

(defining-mpc-build)=
## Building the controller

Once the ingredients are defined, call
[`build_chronos_mpc`](../api/init/build_chronos_mpc.md) to prepare the controller.

The builder prepares the solver workspace, costs, and constraints, resolves
default slack penalties, and initializes the predicted trajectories.

| Argument | Size / default | Role |
| --- | --- | --- |
| `mpc` | Required structure | Controller with dynamics and selected ingredients defined. |
| `x_prev` | Required `nx`-by-1 | Current measured or estimated state for the initial rollout. |
| `u_prev` | Required `nu`-by-1 | Previously applied input; also the starting input for the rollout. |
| `d_in` | `nd`-by-1 or `nd`-by-`L`; default `[]` | Known inputs for dynamics/output. Nonempty data fills `mpc.d`; `[]` leaves existing stored values unchanged. |
| `dh_in` | `ndh`-by-1 or `ndh`-by-`L`; default `[]` | Known inputs for the custom constraints. Nonempty data fills `mpc.dh`; `[]` leaves existing stored values unchanged. |
| `x_ref` | `nx`-by-1; default `[]` | Reference used to guide the initial rollout when a rollout gain `mpc.K` is available. |

Known-input profiles repeat their last short column to $N$ stages and ignore
excess columns, with column 1 corresponding to $k=0$. Supply them only when
their model terms have been configured. Newly initialized active known-input
storage starts at zero.

When a rollout gain and `x_ref` are available, the builder uses them to guide
the initial prediction. Otherwise it starts from the previous input. Building
prepares an initial guess; the optimized solution is computed by `mpc_solve`.
Hard input and rate bounds must be mutually compatible with the initial
rollout and previous input.

Only `mpc`, `x_prev`, and `u_prev` are required. Omit trailing optional
arguments or use `[]` to skip an argument before a later one. The
[builder reference](../api/init/build_chronos_mpc.md) provides call examples.

The builder's `x_ref` is **not** the tracking reference `r_in` or the terminal
reference `xN_ref_in` supplied to subsequent solver calls. The builder also
has no `dz_in` argument; custom-cost known inputs are supplied during solving.

After building, update numerical data through the supported online API.
Changing the horizon, dimensions, or enabled ingredient structure requires
redefining and rebuilding the controller, and regenerating any compiled
solver. Do not add ingredients by calling initialization functions on an
already built online controller.

Continue with [Solving and updating online](online-updates.md) for the
sample-by-sample workflow and ingredient updates. See
[`build_chronos_mpc`](../api/init/build_chronos_mpc.md) and the
[online API](../api/online/index.md) for full function references.
