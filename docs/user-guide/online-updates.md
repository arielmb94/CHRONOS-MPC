# Solving and updating online

After [defining and building the MPC problem](defining-the-mpc-problem.md),
reuse the controller at each sampling instant: update the ingredients that
have changed, solve with the latest state and forecasts, and apply the first
control action. This chapter explains which data to supply for each ingredient;
the [online API reference](../api/online/index.md) gives full argument details.

```{contents}
:local:
:depth: 2
```

## Solving the MPC problem

At each sample, pass the current measured or estimated state as `x_prev` and
the input actually applied at the preceding sample as `u_prev`:

```matlab
% First call any model, cost, or limit updaters needed for this sample.
[u0, mpc, iter] = mpc_solve(mpc, x_prev, u_prev, r_in, xN_ref_in, ...
                           d_in, dz_in, dh_in);
% Apply u0 to the plant; retain mpc for the next sample.
```

`u0` is the next control action, equal to `mpc.u(:,1)`. Always retain the
returned `mpc`: it contains the latest predictions and solver state, which
provide the starting point for the next solve. CHRONOS does not automatically
shift that state or your supplied profiles forward by one sample. Prepare
new forecasts relative to the current sample whenever needed.

### Signals supplied to the solve

| Input | What to supply | Time variation |
| --- | --- | --- |
| `x_prev` | Current state, `nx`-by-1 | One measurement/estimate for this sample |
| `u_prev` | Last applied control, `nu`-by-1 | One value for this sample; used by input-rate and previous-control terms |
| `r_in` | Output tracking reference, `ny` rows | One column or a horizon profile, following the output convention below |
| `xN_ref_in` | Terminal state reference, `nx`-by-1 | One target for stage `N` |
| `d_in` | Known input shared by dynamics and output, `nd` rows | Columns for stages `0` through `N-1` |
| `dz_in` | Known input for the custom cost, `ndz` rows | Columns for stages `0` through `N-1` |
| `dh_in` | Known input for custom constraints, `ndh` rows | Columns for stages `0` through `N-1` |

Pass all argument positions, using `[]` for unused signals. For `r_in`,
`d_in`, `dz_in`, and `dh_in`, `[]` retains the stored values; it does not
reset them to zero. Supply nonempty data only for configured ingredients.
The three known inputs are forecasts supplied by your application, not
variables optimized by CHRONOS.

With terminal ingredients enabled, supply `xN_ref_in` on each solve. If you
configured `xN_ref_is_y=true`, `[]` instead reuses the terminal tracking
reference; this requires the tracking output to be the full state in the
same order. An explicit terminal state reference takes precedence.

### What CHRONOS does during a solve

CHRONOS loads the references and forecasts, refreshes cost gradients and
Hessians requested by the update functions, and updates the dynamics data
using the current state. It then iterates on the predicted states, controls,
and constraint slacks. Model matrices and forecasts remain fixed throughout
that solve; CHRONOS does not recompute an LPV schedule from the new predictions.

The returned `mpc.x` contains states at stages `1` through `N`, and `mpc.u`
contains controls at stages `0` through `N-1`. Configured derived signals,
such as output, custom signals, and input increments, are recomputed from
the iterate. Output/custom interior fields cover stages `1` through `N-1`;
active endpoints use separate `_0` and `_ter` fields.

`iter` reports the number of Newton iterations, not a success status. The
solve can stop at its tolerance, iteration limit, or a sufficiently small
accepted step. See [`mpc_solve`](../api/online/mpc_solve.md) for its interface.

### Time-varying updates and fixed structure

There are two independent kinds of variation: a profile describes changes
**along the prediction horizon**, and an online call replaces that profile
**between solves**. A constant matrix can therefore be replaced by a staged
profile without rebuilding, provided its component dimensions stay the same.

| Supplied data | Constant over the horizon | Horizon profile |
| --- | --- | --- |
| Vector with `n` components | `n`-by-1 | `n`-by-`L` |
| Matrix with `nr` rows and `nc` columns | `nr`-by-`nc` | `nr`-by-`nc`-by-`L` |

A single column/page is reused throughout the horizon. Short profiles repeat
their final column/page; excess stages are ignored. Use full component sizes
for model matrices and cost weights in online calls. Bounds and slack
penalties additionally accept a scalar for all components; a `1`-by-`L`
row does not broadcast to multiple components.

| Ingredient profile | Meaning of supplied entry `j` |
| --- | --- |
| Dynamics matrices, known inputs `d_in`, `dz_in`, `dh_in` | Stage `k=j-1` |
| Control/control-rate weights and bounds | Stage `k=j-1` |
| State bounds and their slack penalties | Stage `k=j` |
| Output/custom matrices, tracking reference, output/custom weights, bounds and penalties | Stage `k=j`; eligible stage-0 rows also reuse entry 1 |

For output/custom profiles, entry `N` supplies eligible terminal terms.
Their known-input forecasts still start at stage 0: for example, the
interior output at stage 1 uses `C(:,:,1)` and `d_in(:,2)` when those
entries are supplied. These are `N`-entry conventions, not independent
`N+1`-entry profiles.

Online updates preserve the horizon, dimensions, enabled terms, bound sides,
and endpoint row selections established during building. In particular,
changing coefficients does not reselect which output/custom rows contribute
at the current or terminal stage. Redefine and rebuild the controller for
structural changes, and regenerate any compiled solver.

Call updaters only for ingredients already configured, and always keep their
returned `mpc`. In selective updates, `[]` leaves that argument unchanged;
zero is numerical data, not a request to remove a term. For the single-weight
tracking, control-rate, and terminal updaters, supply a weight when calling
and skip the call when no update is needed. Use the public update functions
so dependent solver data is refreshed automatically.

## Updating the model

Update a model when its coefficients change, for example after computing a
new LPV schedule. Update a known-input forecast through `mpc_solve` even when
the matrices multiplying it remain unchanged.

### Dynamics

$$x_{k+1}=A_kx_k+B_ku_k+B_{d,k}d_k.$$

```matlab
mpc = update_mpc_dynamics(mpc, A, B, Bd);
```

`A`, `B`, and `Bd` have per-stage sizes `nx`-by-`nx`, `nx`-by-`nu`, and
`nx`-by-`nd`. Each can be constant or use third-dimension pages for stages
`0` through `N-1`. Use `[]` for any unchanged matrix.

The latest state enters as `x_prev`; the known-input forecast enters as
`d_in`, with its first column at stage 0. Updating `Bd` changes how that
forecast affects the state, while updating `d_in` changes the forecast itself.
CHRONOS uses the new matrices and refreshes the dynamics right-hand side
at the next solve. A dynamics update does not recompute a terminal weight
previously obtained from DLQR.

Reference: [`update_mpc_dynamics`](../api/online/update_mpc_dynamics.md).

(online-output-model)=
### Output model

$$y_k=C_kx_k+D_ku_k+D_{d,k}d_k.$$

```matlab
mpc = update_mpc_output_vector(mpc, C, D, Dd);
```

`C`, `D`, and `Dd` have per-stage sizes `ny`-by-`nx`, `ny`-by-`nu`, and
`ny`-by-`nd`. Use constant matrices or pages for stages `1` through `N`,
with eligible current-stage rows reusing page 1. `d_in` supplies the same
known input used by the dynamics, with columns starting at stage 0.

Both tracking and output constraints use the updated output model. Changes
to `C` or `D` request automatic refresh of active tracking derivatives and
the cost Hessian at the next solve. Changing `Dd` or `d_in` changes the
known offset in the output. To change the tracking target, supply `r_in`;
to change output limits, use the {ref}`output-limit updater <online-output-limits>`.

Reference: [`update_mpc_output_vector`](../api/online/update_mpc_output_vector.md).

(online-custom-cost-model)=
### Custom-cost model

$$z_k=C_{z,k}x_k+D_{z,k}u_k+D_{xu,z,k}xu_k+D_{dz,k}d_{z,k}.$$

```matlab
mpc = update_mpc_Custom_cost_vector(mpc, Cz, Dz, Dxuz, Ddz);
```

| Coefficient | Per-stage size | Multiplies |
| --- | --- | --- |
| `Cz` | `nz`-by-`nx` | State `x` |
| `Dz` | `nz`-by-`nu` | Current control `u` |
| `Dxuz` | `nz`-by-`nu` | Previous control `xu` |
| `Ddz` | `nz`-by-`ndz` | Dedicated known input `dz` |

Each coefficient accepts a constant matrix or pages following the custom
profile convention: stages `1` through `N`, with eligible stage-0 rows
reusing page 1. `xu_k=u_{k-1}`; at stage 0 it uses `u_prev`. Supply the
independent `dz_in` forecast to `mpc_solve`, starting at stage 0.

For example, when `z=u-dz` was configured, a new control-reference forecast
only requires changing `dz_in`. Changing how states or controls enter `z`
requires the model updater. CHRONOS refreshes the affected custom-cost
derivatives automatically; the existing `Qz` and `qz` still weight the new
signal. Change them separately under {ref}`custom-cost weights <online-custom-cost-weights>`.

Reference: [`update_mpc_Custom_cost_vector`](../api/online/update_mpc_Custom_cost_vector.md).

(online-custom-constraint-model)=
### Custom-constraint model

$$h_k=C_{h,k}x_k+D_{h,k}u_k+D_{xu,h,k}xu_k+D_{dh,k}d_{h,k}.$$

```matlab
mpc = update_mpc_custom_cnstr_vector(mpc, Ch, Dh, Dxuh, Ddh);
```

| Coefficient | Per-stage size | Multiplies |
| --- | --- | --- |
| `Ch` | `nh`-by-`nx` | State `x` |
| `Dh` | `nh`-by-`nu` | Current control `u` |
| `Dxuh` | `nh`-by-`nu` | Previous control `xu` |
| `Ddh` | `nh`-by-`ndh` | Dedicated known input `dh` |

Use constant matrices or pages for stages `1` through `N`; eligible
stage-0 rows reuse page 1. The previous-control signal uses `u_prev` at
stage 0. Supply `dh_in` separately to `mpc_solve` as an `ndh`-row forecast
starting at stage 0; it is independent of `d_in` and `dz_in`.

Choose the update according to what changed: use `Ch`, `Dh`, `Dxuh`, or
`Ddh` to change the combination being constrained; use `dh_in` to change
its known input; use the {ref}`custom-limit updater <online-custom-constraint-limits>`
to change its permitted range. The next solve evaluates the updated signal
and constraint equations automatically. Existing limits and slack penalties
continue to apply.

Reference: [`update_mpc_custom_cnstr_vector`](../api/online/update_mpc_custom_cnstr_vector.md).

## Updating the costs

Weights change the importance of an objective; references and known inputs
change its target or offset. Supply symmetric positive-semidefinite stage
quadratic weights and a symmetric positive-definite terminal `P`.

### Reference tracking

Tracking penalizes $\tfrac12(r_k-y_k)^TQ_{e,k}(r_k-y_k)$.
Supply a new `r_in` directly to `mpc_solve`. To change the tracking weight:

```matlab
mpc = update_mpc_Tracking_cost(mpc, Qe);
```

`Qe` is `ny`-by-`ny`, or `ny`-by-`ny`-by-`L`; `r_in` is `ny`-by-1 or
`ny`-by-`L`. Both follow the output profile convention: entry 1 applies
at stage 1 and eligible stage-0 rows, and entry `N` supplies the eligible
terminal term. Increasing later weight pages emphasizes tracking further
ahead. CHRONOS refreshes tracking derivatives and the cost Hessian after a
weight update; a reference change enters the tracking error at the next solve.

Reference: [`update_mpc_Tracking_cost`](../api/online/update_mpc_Tracking_cost.md).

### Control effort

Control effort penalizes $\tfrac12u_k^TR_{u,k}u_k+r_{u,k}^Tu_k$.

```matlab
mpc = update_mpc_Control_cost(mpc, Ru, ru);
```

`Ru` uses `nu`-by-`nu` matrix pages and `ru` uses `nu`-component columns,
starting at stage 0 and ending at stage `N-1`. A single matrix or column
applies throughout. Use `[]` to retain either weight. The quadratic and
linear terms must each have been enabled to update them. There is no
separate reference signal in this cost. CHRONOS refreshes the Hessian when
`Ru` changes; an `ru`-only update changes the gradient.

Reference: [`update_mpc_Control_cost`](../api/online/update_mpc_Control_cost.md).

### Control rate

The rate cost penalizes $\tfrac12\Delta u_k^TR_{du,k}\Delta u_k$, where
$\Delta u_k=u_k-u_{k-1}$ and $\Delta u_0=u_0-u_{\mathrm{prev}}$.

```matlab
mpc = update_mpc_Control_Rate_cost(mpc, Rdu);
```

`Rdu` is `nu`-by-`nu`, or a sequence of pages for stages `0` through `N-1`.
Supply the last applied input through `u_prev` on every solve; there is no
separate rate-reference input. These are increments per sample, without
division by the sampling period. The updater requests a cost Hessian refresh.

Reference: [`update_mpc_Control_Rate_cost`](../api/online/update_mpc_Control_Rate_cost.md).

(online-custom-cost-weights)=
### Custom-cost weights

The custom cost penalizes $\tfrac12z_k^TQ_{z,k}z_k+q_{z,k}^Tz_k$.

```matlab
mpc = update_mpc_Custom_cost(mpc, Qz, qz);
```

`Qz` uses `nz`-by-`nz` matrix pages; `qz` uses `nz`-component columns.
Both accept a constant value or a profile for stages `1` through `N`, with
eligible stage-0 terms reusing entry 1. Use `[]` to retain either weight.
CHRONOS refreshes quadratic derivatives and the Hessian for `Qz` changes,
and linear gradients for `qz` changes.

The signal being weighted is defined by the {ref}`custom-cost model <online-custom-cost-model>`.
Its `dz_in` forecast is passed to `mpc_solve`; it can encode a changing
target when the configured model expresses `z` as a deviation from that target.

Reference: [`update_mpc_Custom_cost`](../api/online/update_mpc_Custom_cost.md).

### Terminal cost

The terminal cost is $(x_{N,\mathrm{ref}}-x_N)^TP(x_{N,\mathrm{ref}}-x_N)$.

```matlab
mpc = update_mpc_ter_ingredients(mpc, P);
```

Supply one `nx`-by-`nx` matrix, not a horizon profile. It weights stage `N`
only and can change between solves. Supply its `nx`-by-1 target through
`xN_ref_in`, subject to the tracking-reference shortcut explained above.

CHRONOS updates the terminal derivatives and requests a Hessian refresh.
This updater accepts an already computed `P`; it does not rerun DLQR or
update the initialization rollout gain.

Reference: [`update_mpc_ter_ingredients`](../api/online/update_mpc_ter_ingredients.md).

### Soft-constraint penalties

State, output, and custom constraints allow violations with a linear penalty
on each bound's relaxation. To change that tradeoff, select the configured
constraint family:

```matlab
mpc = update_mpc_slack_cost(mpc, mpc.x_cnstr, qv_min, qv_max);
% Use mpc.y_cnstr for output or mpc.h_cnstr for custom constraints.
```

Each penalty accepts a scalar, a component column, or a profile with `nx`,
`ny`, or `nh` rows for the selected family. State profiles start at stage 1;
output/custom profiles follow their endpoint convention. Use `[]` to retain
one side. Larger positive penalties make violations more costly, while the
bounds remain soft. Control and control-rate limits have no such penalties.

The updater changes the penalty gradients used by the solver directly.
Changing tracking weights does not automatically retune the slack penalties
resolved during building; call this updater when those penalties should change.

Reference: [`update_mpc_slack_cost`](../api/online/update_mpc_slack_cost.md).

## Updating the constraint limits

Limit updates change the allowed range of an existing signal. They do not
change its model, add a bound side, or change whether it is hard or soft.
Use `[]` to leave one side unchanged. Each limit accepts a scalar for all
components, a constant component column, or a component-by-`L` profile;
short profiles hold their last column through the remaining horizon.

### State limits

```matlab
mpc = update_mpc_state_cnstr_limits(mpc, x_min, x_max);
```

These soft limits bound predicted `x` at stages `1` through `N`. Profiles
have `nx` rows; column 1 bounds the first predicted state, not `x_prev`.
For example, later columns can tighten a state envelope in anticipation of
future conditions. No separate constraint matrix is needed. Adjust violation
penalties with `mpc.x_cnstr` in the slack-cost updater.

Reference: [`update_mpc_state_cnstr_limits`](../api/online/update_mpc_state_cnstr_limits.md).

### Control limits

```matlab
mpc = update_mpc_control_cnstr_limits(mpc, u_min, u_max);
```

These hard limits bound `u` at stages `0` through `N-1`. Profiles have
`nu` rows; column 1 limits the control about to be applied. There is no
terminal control or additional constraint matrix. Keep the limits compatible
with any control-rate limits and the previously applied input.

Reference: [`update_mpc_control_cnstr_limits`](../api/online/update_mpc_control_cnstr_limits.md).

### Control-rate limits

```matlab
mpc = update_mpc_control_rate_cnstr_limits(mpc, du_min, du_max);
```

These hard limits bound $\Delta u_k=u_k-u_{k-1}$ at stages `0` through
`N-1`. Profiles have `nu` rows; column 1 limits `u0-u_prev`. Supply the
actual last applied input as `u_prev` to the solve. Limits are input
increments per sample, not derivatives per second; no additional matrix
or rate forecast is supplied.

Reference: [`update_mpc_control_rate_cnstr_limits`](../api/online/update_mpc_control_rate_cnstr_limits.md).

(online-output-limits)=
### Output limits

```matlab
mpc = update_mpc_output_cnstr_limits(mpc, y_min, y_max);
```

These soft limits bound the output `y=C*x+D*u+Dd*d`. Profiles have `ny`
rows and entries for stages `1` through `N`; eligible stage-0 rows reuse
column 1. CHRONOS applies the endpoint selections fixed during building.
For example, lower later upper-bound columns to anticipate a tighter output
limit. Update `C`, `D`, or `Dd` with the {ref}`output-model updater <online-output-model>`,
and supply the known input through `d_in`. Adjust violation penalties with
`mpc.y_cnstr` in the slack-cost updater.

Reference: [`update_mpc_output_cnstr_limits`](../api/online/update_mpc_output_cnstr_limits.md).

(online-custom-constraint-limits)=
### Custom-constraint limits

```matlab
mpc = update_mpc_custom_cnstr_limits(mpc, h_min, h_max);
```

These soft limits bound `h=Ch*x+Dh*u+Dxuh*xu+Ddh*dh`. Profiles have `nh`
rows and entries for stages `1` through `N`; eligible stage-0 rows reuse
column 1. Terminal limits apply only to the state-only rows selected during
building. Update the coefficients with the {ref}`custom-model updater <online-custom-constraint-model>`,
and supply the dedicated known-input forecast through `dh_in`.

For example, if `h` is a configured sum of states, changing its upper limit
requires only `update_mpc_custom_cnstr_limits(mpc, [], h_max)`. Changing
which combination of states is constrained requires updating `Ch` as well.
Adjust violation penalties with `mpc.h_cnstr` in the slack-cost updater.
The next solve uses the updated signal and limits automatically.

Reference: [`update_mpc_custom_cnstr_limits`](../api/online/update_mpc_custom_cnstr_limits.md).
