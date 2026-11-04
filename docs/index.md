# CHRONOS

```{image} ../logo.png
:alt: CHRONOS logo
:width: 400px
:align: center
```

**Predictive control for changing operating conditions.**

CHRONOS is a Model Predictive Control (MPC) solver for nonlinear systems through the use of linear
parameter-varying (LPV) and linear time-varying (LTV) models. Define your
prediction model, control objectives, and operating limits, CHRONOS assembles
and solves the MPC problem to compute the next control action.

For nonlinear plants, express the dynamics and signals in LPV form and
update their coefficients as operating conditions change. Each solve uses
a fixed model schedule, with coefficients that can vary along the prediction
horizon. The online solver is designed for C/C++ generation with MATLAB Coder.

## Start here

**Learn the approach.** Follow the [User Guide](user-guide/index.md) from the
[CHRONOS overview](user-guide/chronos-overview.md) through nonlinear-to-LPV
modeling, controller configuration, and online operation. It ends with a
complete closed-loop simulation.

**Start with an application.** [Your first controller](user-guide/your-first-controller.md)
walks through a two-tank system: derive the prediction model, track a water-level
reference, configure operating limits, and simulate the nonlinear plant. The
chapter includes MATLAB path setup and a complete controller script.

## Shape the controller around your plant

- **Predict changing dynamics.** Use stage-varying model coefficients and
  known-input forecasts to describe the operating conditions ahead.
- **Choose what matters.** Combine output tracking, control-effort and
  control-rate penalties, terminal objectives, and costs on custom signals.
- **Express operating limits.** Set bounds on controls and control
  increments, on states, outputs, and custom signals.

The [modeling chapter](user-guide/from-nonlinear-problems-to-lpv-mpc.md)
explains how to prepare the model and signals. The
[configuration chapter](user-guide/configuring-the-mpc-controller.md)
shows how to combine them into a controller.

## From definition to the control loop

Configure and build the controller once. At each sampling instant, update
any model coefficients, weights, or limits that have changed, then solve
using the current state, previously applied input, references, and forecasts.
Apply the first control action and retain the returned controller structure
for the next sample.

[Solving and updating the MPC](user-guide/solving-and-updating-the-mpc.md)
connects these steps and explains which data can change online.

## Find a function

The [API reference](api/index.md) provides function signatures, argument
descriptions, and usage examples.

| Task | Reference |
| --- | --- |
| Define dynamics, outputs, objectives, and constraints; build the controller | [Initialization API](api/init/index.md) |
| Update configured problem data between samples | [Online API](api/online/index.md) |
| Compute the next control action | [`mpc_solve`](api/online/mpc_solve.md) |

## Source and further examples

Get CHRONOS from the [source repository](https://github.com/arielmb94/CHRONOS-MPC).
Explore the [application examples](https://github.com/arielmb94/CHRONOS-MPC/tree/docs/examples)
for additional plant models and controller designs. CHRONOS is distributed
under the [MIT license](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/LICENSE).

For questions or collaboration, contact
[Ariel Medero Borrell](https://www.linkedin.com/in/ariel-medero-borrell).

```{toctree}
:hidden:
:maxdepth: 2

Home <self>
User Guide <user-guide/index>
API reference <api/index>
```
