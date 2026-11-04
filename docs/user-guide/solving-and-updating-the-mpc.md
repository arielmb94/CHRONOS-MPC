# Solving and updating the MPC

[Configuring the MPC controller](configuring-the-mpc-controller.md) prepares
an `mpc` structure with your prediction model, costs, and constraints. This
chapter shows how to reuse that controller at each sampling instant: update
its data, solve with the current state and references, and apply the
control action. The [online API reference](../api/online/index.md) provides
the full argument descriptions.

The following sections explain each ingredient's online interface. Continue
to [Your first controller](your-first-controller.md) for a plant-based
walkthrough that brings configuration and online operation together.

```{contents}
:local:
:depth: 2
```

## CHRONOS online control workflow

At each sampling instant:

1. Obtain the current state estimate and the input applied at the previous sample.
2. Prepare the current references, known-input forecasts, and any new LPV schedule.
3. Update the configured models, weights, or limits that have changed using the update API functions.
4. Solve with `mpc_solve`, apply the first control action, and retain the returned `mpc`.

## Solving the MPC problem

At each sampling instant, solve the current problem with
[`mpc_solve`](../api/online/mpc_solve.md). Pass your `mpc` structure, which
holds the configured problem and the latest predictions and solver state.
Keep the updated structure returned by the solver for the next sample.

When calling `mpc_solve`, supply the following signals:

| Input | Size | Meaning | When to use |
| --- | --- | --- | --- |
| `x_prev` | $n_x\times 1$ | State at the start of the current sample; initializes the prediction | Every solve |
| `u_prev` | $n_u\times 1$ | Control action applied during the previous sample | Every solve; use the value actually applied by the plant or actuator |
| `r_in` | $n_y\times L$ | Reference for the configured output tracking cost | When tracking is configured and its reference needs updating; pass `[]` to retain the stored reference |
| `xN_ref_in` | $n_x\times 1$ | Target for the terminal predicted state | When terminal ingredients are configured; required each solve unless `xN_ref_is_y` reuses the terminal tracking reference |
| `d_in` | $n_d\times L$ | Known-input forecast shared by the dynamics and output models, starting at stage 0 | When a configured dynamics or output model uses this input and its forecast needs updating; pass `[]` to retain the stored forecast |
| `dz_in` | $n_{dz}\times L$ | Known-input forecast for the configured custom-cost signal, starting at stage 0 | When the custom-cost signal uses this input and its forecast needs updating; pass `[]` to retain the stored forecast |
| `dh_in` | $n_{dh}\times L$ | Known-input forecast for the configured custom-constraint signal, starting at stage 0 | When the custom-constraint signal uses this input and its forecast needs updating; pass `[]` to retain the stored forecast |

For `r_in`, `d_in`, `dz_in`, and `dh_in`, `L` is the number of supplied
columns, from 1 to `N`. Use `L=1` for a constant reference or forecast. If
`L<N`, CHRONOS repeats the last supplied column over the remaining stages;
columns beyond `N` are ignored. Each signal can have its own number of
supplied columns `L`.

Pass all argument positions, using `[]` for unused signals. For `r_in`,
`d_in`, `dz_in`, and `dh_in`, `[]` retains the stored values; it does not
reset them to zero. Supply nonempty data only for configured ingredients.
The three known inputs are forecasts supplied by your application, not
variables optimized by CHRONOS.

With terminal ingredients enabled, supply `xN_ref_in` on each solve. If you
configured `xN_ref_is_y=true`, `[]` instead reuses the terminal tracking
reference; this requires the tracking output to be the full state in the
same order. An explicit terminal state reference takes precedence.

`mpc_solve` returns the following data:

| Output | Meaning |
| --- | --- |
| `u0` | First control action in the solution; apply it to the plant. |
| `mpc` | Updated controller structure with the latest predictions and solver state; retain it for the next solve. |
| `iter` | Number of Newton iterations used. |

### Preserve the `mpc` structure

Keep the returned `mpc` after every solver or online update call. Each call
stores its changes in that structure, so pass the latest version to the
next update or solve.

## Updating the MPC problem online

The CHRONOS online API lets you update the numerical data of configured
models, cost terms, and constraints between solves. The update functions
also set internal update flags when dependent quantities need to be
recomputed during the setup phase of the next `mpc_solve` call. Use these functions
rather than modifying `mpc` fields directly, so CHRONOS can also compute
any related data needed for the next solve.

Online updates change numerical values only; the problem dimensions remain
fixed. During the build phase, CHRONOS allocates and precomputes internal
data for the dimensions specified when defining each ingredient. To change
those dimensions, redefine and rebuild the controller.

| What to update | When to use |
| --- | --- |
| {ref}`Models <online-models>` | Change prediction dynamics or the coefficients defining output, custom-cost, or custom-constraint signals, for example after evaluating a new LPV schedule. |
| {ref}`Cost terms <online-costs>` | Change objective weights or tracking and terminal targets. |
| {ref}`Constraint limits <online-constraint-limits>` | Change the permitted ranges of states, controls, control increments, outputs, or custom signals, or the penalties for violating soft bounds. |

(online-models)=
## Updating models

The model coefficients may change as your system moves to a new operating
point or as you evaluate a new LPV schedule. Use the model update functions
below to pass these new coefficients to CHRONOS before solving again.

You only need to supply coefficients that have changed: pass `[]` for the
others. Each updater keeps the dimensions, enabled terms, and stage
selections established when the controller was built.

| Model type | When to update it |
| --- | --- |
| {ref}`Dynamics <online-dynamics-model>` | Change the prediction-model matrices, for example after evaluating a new LPV schedule. |
| {ref}`Output model <online-output-model>` | Change the coefficients of the output used by tracking and output constraints. |
| {ref}`Custom-cost signal <online-custom-cost-model>` | Change the coefficients defining the custom signal being penalized. |
| {ref}`Custom-constraint signal <online-custom-constraint-model>` | Change the coefficients defining the custom signal being constrained. |

(online-dynamics-model)=
### Dynamics

Supply `A`, `B`, and `Bd` to
[`update_mpc_dynamics`](../api/online/update_mpc_dynamics.md).
Use discrete-time matrices at the sample time used to configure the
controller. For LPV models, evaluate and discretize the model before the
update, as described in
[From nonlinear to LPV models](from-nonlinear-problems-to-lpv-mpc.md).

**Model coefficients**

| Input | Size | Meaning |
| --- | --- | --- |
| `A` | $n_x\times n_x\times L$ | State matrix |
| `B` | $n_x\times n_u\times L$ | Control-input matrix |
| `Bd` | $n_x\times n_d\times L$ | Known-input matrix |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
matrix (`L=1`) for a constant coefficient. If `L<N`, the last matrix is
repeated over the remaining stages; any pages beyond `N` are ignored.

The first supplied matrix applies at stage 0. Each following page applies
to the next stage, through stage `N-1`.

If only the known-input forecast changes, pass the new values as `d_in`
directly to [`mpc_solve`](../api/online/mpc_solve.md). The `Bd` matrix
describes how the model uses that input, so it only needs updating when
those coefficients change.

(online-output-model)=
### Output model

Supply `C`, `D`, and `Dd` to
[`update_mpc_output_vector`](../api/online/update_mpc_output_vector.md).
Both tracking and output constraints use the updated coefficients.

**Model coefficients**

| Input | Size | Meaning |
| --- | --- | --- |
| `C` | $n_y\times n_x\times L$ | Coefficient on the state |
| `D` | $n_y\times n_u\times L$ | Coefficient on the current control input |
| `Dd` | $n_y\times n_d\times L$ | Coefficient on the shared known input |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
matrix (`L=1`) for a constant coefficient. If `L<N`, the last matrix is
repeated over the remaining stages; any pages beyond `N` are ignored.

Profiles start at stage 1. CHRONOS also updates the related data for the
stage-0 and terminal output rows selected during building. Stage-0 rows
use the control input through `D`, while terminal rows depend exclusively
on the states.

For a new known-input forecast or tracking reference, pass `d_in` or `r_in`
directly to [`mpc_solve`](../api/online/mpc_solve.md). These inputs can
change without updating the model coefficients. The `d_in` forecast is
shared with the dynamics and starts at stage $k=0$.

(online-custom-cost-model)=
### Custom-cost signal

Supply `Cz`, `Dz`, `Dxuz`, and `Ddz` to
[`update_mpc_Custom_cost_vector`](../api/online/update_mpc_Custom_cost_vector.md).

**Model coefficients**

| Input | Size | Meaning |
| --- | --- | --- |
| `Cz` | $n_z\times n_x\times L$ | Coefficient on the state |
| `Dz` | $n_z\times n_u\times L$ | Coefficient on the current control input |
| `Dxuz` | $n_z\times n_u\times L$ | Coefficient on the previous control input |
| `Ddz` | $n_z\times n_{dz}\times L$ | Coefficient on the dedicated known input |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
matrix (`L=1`) for a constant coefficient. If `L<N`, the last matrix is
repeated over the remaining stages; any pages beyond `N` are ignored.

Profiles start at stage 1. Any signal rows included at stage 0 reuse the
first page, while state-only rows included at stage `N` use the terminal
page of `Cz`.

For a new known-input forecast, pass `dz_in` directly to
[`mpc_solve`](../api/online/mpc_solve.md), with its first column at stage 0.
This forecast is independent of `d_in` and `dh_in`. For example, if your
model uses `dz_in` for a changing target, supply the new target there while
keeping the model coefficients unchanged. The previous-control term uses
the `u_prev` input at stage 0.

CHRONOS automatically computes any related cost data needed for the next
solve. The `Qz` and `qz` weights keep their current values; change them
separately with the {ref}`custom-cost weight updater <online-custom-cost-weights>`.

(online-custom-constraint-model)=
### Custom-constraint signal

Supply `Ch`, `Dh`, `Dxuh`, and `Ddh` to
[`update_mpc_custom_cnstr_vector`](../api/online/update_mpc_custom_cnstr_vector.md).

**Model coefficients**

| Input | Size | Meaning |
| --- | --- | --- |
| `Ch` | $n_h\times n_x\times L$ | Coefficient on the state |
| `Dh` | $n_h\times n_u\times L$ | Coefficient on the current control input |
| `Dxuh` | $n_h\times n_u\times L$ | Coefficient on the previous control input |
| `Ddh` | $n_h\times n_{dh}\times L$ | Coefficient on the dedicated known input |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
matrix (`L=1`) for a constant coefficient. If `L<N`, the last matrix is
repeated over the remaining stages; any pages beyond `N` are ignored.

Profiles start at stage 1. Any signal rows included at stage 0 reuse the
first page, while state-only rows included at stage `N` use the terminal
page of `Ch`.

For a new known-input forecast, pass `dh_in` directly to
[`mpc_solve`](../api/online/mpc_solve.md), with its first column at stage 0.
This forecast is independent of `d_in` and `dz_in`. If only its values
change, keep the model coefficients unchanged. The previous-control term
uses the `u_prev` input at stage 0.

The next solve uses the updated signal model with the existing bounds and
violation penalties. Change those separately with the
{ref}`custom-limit updater <online-custom-constraint-limits>` and
{ref}`slack-penalty updater <online-slack-penalties>`.

(online-costs)=
## Updating cost terms

As operating conditions or priorities change, you may want to adjust the
weights of the configured costs. Use the cost update functions below to
pass the new weights to CHRONOS before solving again.

Only cost terms enabled during definition can be updated. The signal
models and problem dimensions stay as configured, and CHRONOS computes any
related data needed for the next solve.

| Cost type | When to update it |
| --- | --- |
| {ref}`Tracking <online-tracking>` | Change the target or the importance of output errors. |
| {ref}`Control effort <online-control-effort>` | Change the penalty on input magnitude or its linear cost. |
| {ref}`Control rate <online-control-rate>` | Change how strongly the controller discourages input increments. |
| {ref}`Terminal state <online-terminal-cost>` | Change the final-state target or weight. |
| {ref}`Custom cost <online-custom-cost-weights>` | Change the weights on a configured custom signal. |

Use symmetric positive semidefinite matrices for stage quadratic weights
and a symmetric positive definite matrix for the terminal weight `P`.

To change the penalties for soft-constraint violations, see
{ref}`slack penalties <online-slack-penalties>` in the constraint update section.

(online-tracking)=
### Tracking

Supply the new tracking weight `Qe` to
[`update_mpc_Tracking_cost`](../api/online/update_mpc_Tracking_cost.md).

**Cost weights**

| Input | Size | Meaning |
| --- | --- | --- |
| `Qe` | $n_y\times n_y\times L$ | Quadratic weight on the output tracking error |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
matrix (`L=1`) for a constant weight. If `L<N`, the last matrix is repeated
over the remaining stages; any pages beyond `N` are ignored.

Profiles start at stage 1. Output rows included at stage 0 reuse the first
weight matrix, while state-only rows included at stage `N` use the terminal
weight matrix. For example, increasing the later weights gives more
importance to tracking errors further along the horizon.

If only the tracking reference changes, pass the new values as `r_in`
directly to [`mpc_solve`](../api/online/mpc_solve.md). You can keep `Qe`
unchanged when updating the reference.

(online-control-effort)=
### Control effort

Supply `Ru`, `ru`, or both to
[`update_mpc_Control_cost`](../api/online/update_mpc_Control_cost.md).
Pass `[]` for either weight you want to keep. The quadratic and linear
terms can be updated independently, provided each was enabled during
definition.

**Cost weights**

| Input | Size | Meaning |
| --- | --- | --- |
| `Ru` | $n_u\times n_u\times L$ | Quadratic weight on the control input |
| `ru` | $n_u\times L$ | Linear weight on the control input |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
matrix for `Ru` or a single column for `ru` (`L=1`) for a constant weight.
If `L<N`, the last matrix or column is repeated over the remaining stages;
any pages or columns beyond `N` are ignored.

The first matrix of `Ru` and first column of `ru` apply at stage 0, to the
control action about to be computed. Profiles continue through stage
`N-1`. To change only the linear cost, supply the new `ru` and pass `[]`
for `Ru`.

(online-control-rate)=
### Control rate

Supply the new weight `Rdu` to
[`update_mpc_Control_Rate_cost`](../api/online/update_mpc_Control_Rate_cost.md).

**Cost weights**

| Input | Size | Meaning |
| --- | --- | --- |
| `Rdu` | $n_u\times n_u\times L$ | Quadratic weight on control-input changes per sample |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
matrix (`L=1`) for a constant weight. If `L<N`, the last matrix is repeated
over the remaining stages; any pages beyond `N` are ignored.

The first matrix applies at stage 0, followed by one matrix per stage
through `N-1`. Choose the weights for input increments per sample; CHRONOS
does not divide those increments by the sampling period.

Supply the last applied control input as `u_prev` directly to
[`mpc_solve`](../api/online/mpc_solve.md) on every solve. A change in
`u_prev` does not require updating `Rdu`; it changes the starting point
from which the first input increment is measured.

(online-custom-cost-weights)=
### Custom cost

Supply `Qz`, `qz`, or both to
[`update_mpc_Custom_cost`](../api/online/update_mpc_Custom_cost.md).
Pass `[]` for either weight you want to keep. The quadratic and linear
terms can be updated independently, provided each was enabled during
definition.

**Cost weights**

| Input | Size | Meaning |
| --- | --- | --- |
| `Qz` | $n_z\times n_z\times L$ | Quadratic weight on the custom signal |
| `qz` | $n_z\times L$ | Linear weight on the custom signal |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
matrix for `Qz` or a single column for `qz` (`L=1`) for a constant weight.
If `L<N`, the last matrix or column is repeated over the remaining stages;
any pages or columns beyond `N` are ignored.

Profiles start at stage 1. Signal rows included at stage 0 reuse the first
weights, while state-only rows included at stage `N` use the terminal
weights.

These updates keep the custom signal model unchanged. To change its
coefficients, use the {ref}`custom-cost model updater <online-custom-cost-model>`.
If only its known-input forecast or a target encoded in that forecast
changes, pass the new `dz_in` directly to
[`mpc_solve`](../api/online/mpc_solve.md).

(online-terminal-cost)=
### Terminal state

Supply the new terminal weight `P` to
[`update_mpc_ter_ingredients`](../api/online/update_mpc_ter_ingredients.md).

**Cost weights**

| Input | Size | Meaning |
| --- | --- | --- |
| `P` | $n_x\times n_x$ | Quadratic weight on the terminal state error |

Supply one symmetric positive definite matrix for stage `N`. This weight
has no time-varying dimension within the horizon, but you can change it
between solves. If your weight was originally obtained with DLQR, compute
the new `P` before calling the updater; the updater does not rerun DLQR.

Supply the terminal target as an `nx`-by-1 `xN_ref_in` directly to
[`mpc_solve`](../api/online/mpc_solve.md) on every solve. If you configured
`xN_ref_is_y=true`, pass `[]` to reuse the terminal tracking reference
instead. This shortcut requires the tracking output to be the full state
in the same order. An explicit `xN_ref_in` takes precedence.

(online-constraint-limits)=
## Updating constraint limits

As operating conditions change, you may need to adjust the permitted range
of a configured signal. Use the limit update functions below to pass the
new bounds to CHRONOS before solving again. You can also update the
{ref}`slack penalties <online-slack-penalties>` for soft bounds.

Only lower and upper bounds enabled during definition can be updated.
Pass `[]` for either bound you want to keep; it does not remove that bound.
The dimensions, signal models, and hard or soft behavior stay as configured.

| Constraint type | When to update it | Enforcement |
| --- | --- | --- |
| {ref}`State <online-state-limits>` | Change the predicted state operating range. | Soft |
| {ref}`Control input <online-control-limits>` | Change actuator limits. | Hard |
| {ref}`Control rate <online-rate-limits>` | Change the permitted input increment per sample. | Hard |
| {ref}`Output <online-output-limits>` | Change limits on the configured output. | Soft |
| {ref}`Custom signal <online-custom-constraint-limits>` | Change limits on a configured custom quantity. | Soft |

(online-state-limits)=
### State

Supply `x_min`, `x_max`, or both to
[`update_mpc_state_cnstr_limits`](../api/online/update_mpc_state_cnstr_limits.md).

**Constraint limits**

| Input | Size | Meaning |
| --- | --- | --- |
| `x_min` | $n_x\times L$ | Lower bound on the predicted state |
| `x_max` | $n_x\times L$ | Upper bound on the predicted state |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
column (`L=1`) for limits that stay constant across the horizon, or a
scalar to use one value for all components and stages. If `L<N`, the last
column is repeated over the remaining stages; columns beyond `N` are
ignored.

The first column applies to the first predicted state at stage 1, rather
than the current state `x_prev`. Profiles continue through the terminal
state at stage `N`.

To update the violation penalties too, select `mpc.x_cnstr` in the
{ref}`slack-penalty updater <online-slack-penalties>`.

(online-control-limits)=
### Control input

Supply `u_min`, `u_max`, or both to
[`update_mpc_control_cnstr_limits`](../api/online/update_mpc_control_cnstr_limits.md).

**Constraint limits**

| Input | Size | Meaning |
| --- | --- | --- |
| `u_min` | $n_u\times L$ | Lower bound on the control input |
| `u_max` | $n_u\times L$ | Upper bound on the control input |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
column (`L=1`) for limits that stay constant across the horizon, or a
scalar to use one value for all components and stages. If `L<N`, the last
column is repeated over the remaining stages; columns beyond `N` are
ignored.

The first column applies at stage 0, to the control action about to be
computed. Profiles continue through stage `N-1`. Keep the new limits
compatible with any {ref}`control-rate limits <online-rate-limits>` and
the previously applied control input.

(online-rate-limits)=
### Control rate

Supply `du_min`, `du_max`, or both to
[`update_mpc_control_rate_cnstr_limits`](../api/online/update_mpc_control_rate_cnstr_limits.md).

**Constraint limits**

| Input | Size | Meaning |
| --- | --- | --- |
| `du_min` | $n_u\times L$ | Lower bound on the control-input change per sample |
| `du_max` | $n_u\times L$ | Upper bound on the control-input change per sample |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
column (`L=1`) for limits that stay constant across the horizon, or a
scalar to use one value for all components and stages. If `L<N`, the last
column is repeated over the remaining stages; columns beyond `N` are
ignored.

Profiles start at stage 0 and continue through stage `N-1`. The first
column limits the change from the last applied input to the next control
action. Supply that last applied input as `u_prev` directly to
[`mpc_solve`](../api/online/mpc_solve.md) on every solve.

Set these limits in units of input increment per sample; CHRONOS does not
divide the increment by the sampling period. Keep the new limits
compatible with the {ref}`control-input bounds <online-control-limits>`
and the value actually applied at the previous sample.

(online-output-limits)=
### Output

Supply `y_min`, `y_max`, or both to
[`update_mpc_output_cnstr_limits`](../api/online/update_mpc_output_cnstr_limits.md).

**Constraint limits**

| Input | Size | Meaning |
| --- | --- | --- |
| `y_min` | $n_y\times L$ | Lower bound on the output |
| `y_max` | $n_y\times L$ | Upper bound on the output |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
column (`L=1`) for limits that stay constant across the horizon, or a
scalar to use one value for all components and stages. If `L<N`, the last
column is repeated over the remaining stages; columns beyond `N` are
ignored.

Profiles start at stage 1. Output rows included at stage 0 reuse the first
column, while state-only rows included at stage `N` use the terminal
column. For example, you can lower later `y_max` columns to describe a
tighter limit expected further along the horizon.

To change `C`, `D`, or `Dd`, use the
{ref}`output-model updater <online-output-model>`. If only the known-input
forecast changes, pass the new `d_in` directly to
[`mpc_solve`](../api/online/mpc_solve.md); the existing limits continue to
apply. To update violation penalties, select `mpc.y_cnstr` in the
{ref}`slack-penalty updater <online-slack-penalties>`.

(online-custom-constraint-limits)=
### Custom signal

Supply `h_min`, `h_max`, or both to
[`update_mpc_custom_cnstr_limits`](../api/online/update_mpc_custom_cnstr_limits.md).

**Constraint limits**

| Input | Size | Meaning |
| --- | --- | --- |
| `h_min` | $n_h\times L$ | Lower bound on the custom signal |
| `h_max` | $n_h\times L$ | Upper bound on the custom signal |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
column (`L=1`) for limits that stay constant across the horizon, or a
scalar to use one value for all components and stages. If `L<N`, the last
column is repeated over the remaining stages; columns beyond `N` are
ignored.

Profiles start at stage 1. Signal rows included at stage 0 reuse the first
column, while state-only rows included at stage `N` use the terminal
column. To change only the upper bound, supply the new `h_max` and pass
`[]` for `h_min`.

These updates keep the custom signal model unchanged. To change its
coefficients, use the
{ref}`custom-model updater <online-custom-constraint-model>`. If only its
known-input forecast changes, pass the new `dh_in` directly to
[`mpc_solve`](../api/online/mpc_solve.md). To update violation penalties,
select `mpc.h_cnstr` in the
{ref}`slack-penalty updater <online-slack-penalties>`.

(online-slack-penalties)=
### Slack penalties

Supply a constraint selection `cnstr` and the new penalties `qv_min` and
`qv_max` to
[`update_mpc_slack_cost`](../api/online/update_mpc_slack_cost.md).
Pass `[]` for either penalty you want to keep. The corresponding lower
or upper bound must already be configured.

**Constraint selection**

| Supply as `cnstr` | Constraint family | Number of components `n` |
| --- | --- | --- |
| `mpc.x_cnstr` | Predicted states | $n_x$ |
| `mpc.y_cnstr` | Output | $n_y$ |
| `mpc.h_cnstr` | Custom signal | $n_h$ |

**Violation penalties**

| Input | Size | Meaning |
| --- | --- | --- |
| `qv_min` | $n\times L$ | Linear penalty for violating the lower bound |
| `qv_max` | $n\times L$ | Linear penalty for violating the upper bound |

Here, `L` is the number of supplied stages, from 1 to `N`. Use a single
column (`L=1`) for penalties that stay constant across the horizon, or a
scalar to use one value for all components and stages. If `L<N`, the last
column is repeated over the remaining stages; columns beyond `N` are
ignored.

State penalty profiles start at stage 1 and continue through stage `N`.
Output and custom-signal profiles also start at stage 1, with the first
column reused for any configured stage-0 rows and the terminal column used
for configured stage-`N` rows.

Increasing a positive penalty makes violations more costly while keeping
the bound soft. Control and control-rate limits have no slack penalties.
Changing tracking weights does not automatically change the default slack
penalties selected during building; use this updater when those penalties
should change too.

## Next step

You now know how to
[define your MPC problem](configuring-the-mpc-controller.md) and update its
data online. Continue with [Your first controller](your-first-controller.md)
to put these steps together in a complete closed-loop simulation, from
configuring the controller to updating, solving, and applying the control
action.
