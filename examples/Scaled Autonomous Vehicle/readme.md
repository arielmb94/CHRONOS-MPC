# Scaled Autonomous Vehicle Example

## Folder contents

- *main_SAV_Simulator.m*: loads the simulator data, builds the controllers,
  launches the simulation, and plots the results.
- *sav_lat_controller_init.m*: defines and builds the CHRONOS lateral MPC.
- *update_BM.m*: evaluates the LPV bicycle model at the measured longitudinal
  speed.
- *SAV Simulator*: contains the remaining models and software blocks required
  to simulate the complete Scaled Autonomous Vehicle (SAV).

The online lateral-control code runs inside a MATLAB Function block in the
simulator. The initialized `mpc_initial` structure is copied into a persistent
controller state when the simulation starts.

## Example overview

This example controls the trajectory of a 1:15-scale autonomous vehicle. The
high-fidelity simulator reproduces both the identified vehicle dynamics and
the software architecture of the real platform, including path planning,
longitudinal control, motor control, and lateral control.

Only the lateral-control block is modified in this example. The path planner
provides a yaw-rate reference, and CHRONOS computes the steering command needed
to track it:

```text
Planned path ---> Pure Pursuit ---> Yaw-rate reference
                                            |
                                            v
                                  CHRONOS lateral MPC
                                            |
                                            v
                                    Steering command ---> SAV simulator
```

## Simulator fidelity and model mismatch

The simulator is considerably more detailed than the prediction model used by
the MPC:

| Aspect | SAV simulator | MPC prediction model |
| --- | --- | --- |
| Vehicle dynamics | Coupled nonlinear lateral and longitudinal dynamics with identified nonlinear tyre forces | Linear bicycle structure scheduled by $v_x$ |
| Steering actuator | Slow second-order servo with a 0.18 s input delay | Ignored |
| Propulsion | Rear-motor dynamics, delays, noisy low-resolution encoders, current sensing, observers, and motor controllers | Not included |
| Software | Position and velocity processing, Pure-Pursuit path planning, differential action, and longitudinal control | Receives the resulting yaw-rate reference and measured states |

Consequently, the MPC does not predict several important effects present in
the closed loop: nonlinear tyre behaviour, lateral-longitudinal coupling, or
the steering-actuator dynamics and delay. This intentional mismatch makes the
example a demanding test of CHRONOS. Successful trajectory and yaw-rate
tracking illustrate its practical robustness when controlling a realistic
nonlinear system with significant unmodelled behaviour.

## LPV lateral model

The controller uses the single-track bicycle model

```math
\begin{bmatrix}
\dot v_y\\
\ddot\psi
\end{bmatrix}
=
\begin{bmatrix}
-\frac{C_f+C_r}{m v_x} &
-v_x-\frac{C_fl_f-C_rl_r}{m v_x}\\
-\frac{C_fl_f-C_rl_r}{I_zv_x} &
-\frac{C_fl_f^2+C_rl_r^2}{I_zv_x}
\end{bmatrix}
\begin{bmatrix}
v_y\\
\dot\psi
\end{bmatrix}
+
\begin{bmatrix}
\frac{C_f}{m}\\
\frac{C_fl_f}{I_z}
\end{bmatrix}\delta.
```

The model signals are:

| Symbol | Meaning | Unit |
| --- | --- | --- |
| $v_x$ | Longitudinal speed and measured scheduling variable | m/s |
| $v_y$ | Lateral speed | m/s |
| $\psi$ | Yaw angle | rad |
| $\dot\psi$ | Yaw rate | rad/s |
| $\delta$ | Steering command | rad |

The front and rear cornering stiffnesses, $C_f$ and $C_r$, also vary with
$v_x$. Consequently, `update_BM` evaluates the LPV matrices at the current
longitudinal speed before every MPC solve.

## MPC definition

The MPC state, tracked output, and reference are

$$
x=\begin{bmatrix}v_y&\dot\psi\end{bmatrix}^T,
\qquad y=\dot\psi,
\qquad y^{ref}=\dot\psi^{ref}.
$$

Since only yaw rate is tracked, the output matrix is explicitly defined as
$C=\begin{bmatrix}0&1\end{bmatrix}$. The steering command is bounded by

$$
-0.7\leq\delta_k\leq0.7\ \text{rad}.
$$

With $N=50$ and $T_s=0.02$ s, the prediction horizon covers one second. The
configured objective is

$$
J=\frac{1}{2}\sum_{k=1}^{N}
(y_k^{ref}-y_k)^TQ_e(y_k^{ref}-y_k)
+\frac{1}{2}\sum_{k=0}^{N-1}\delta_k^TR_u\delta_k,
$$


## Online use

At each lateral-control sample, the MATLAB Function block:

1. Reads $v_x$, $v_y$, the yaw rate, its reference, and the previous steering
   command.
2. Evaluates and discretizes the bicycle model at the measured $v_x$.
3. Updates the CHRONOS prediction matrices.
4. Solves the MPC problem for the steering command.
5. Stores the updated MPC structure for the next sample.

The reported controller time includes both the LPV-model update and
`mpc_solve`. The simulator also records the Newton iteration count. The sample
time used by the launcher, initializer, and online block must remain equal.

Run *main_SAV_Simulator.m* to initialize the complete simulator, execute the
closed-loop test, and plot the trajectory, lateral-control signals, computation
time, and iterations.

## References

[1] Medero Borrell, A. (2023). *LPV lateral control of autonomous and automated
vehicles*. Universitat Politècnica de Catalunya.
https://doi.org/10.5821/dissertation-2117-412168

[2] Ariel M. Borrell, Vicenç Puig, and Olivier Sename. Fixed-structure
parameter-dependent state feedback controller: A scaled autonomous vehicle
path-tracking application. *Control Engineering Practice*, 2024.
https://doi.org/10.1016/j.conengprac.2024.105911
