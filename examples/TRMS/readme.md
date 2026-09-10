# Twin Rotor MIMO System Examples

The Twin Rotor MIMO System (TRMS) is a nonlinear plant inspired by helicopter
pitch and yaw dynamics. These examples use CHRONOS MPC-LPV controllers to solve
the same angle-tracking problem with three different control architectures.

## Plant model

The nonlinear plant used in every simulation is implemented in *TRMS.m*. Its
state and physical control vectors are

$$
x=\begin{bmatrix}
\omega_h&\Omega_h&\theta_h&\omega_v&\Omega_v&\theta_v
\end{bmatrix}^T,
\qquad
u=\begin{bmatrix}u_h&u_v\end{bmatrix}^T.
$$

| Symbol | Meaning | Unit |
| --- | --- | --- |
| $\omega_h$ | Tail-rotor angular speed | rad/s |
| $\Omega_h$ | Horizontal yaw rate | rad/s |
| $\theta_h$ | Horizontal yaw angle | rad |
| $\omega_v$ | Main-rotor angular speed | rad/s |
| $\Omega_v$ | Vertical pitch rate | rad/s |
| $\theta_v$ | Vertical pitch angle | rad |
| $u_h$ | Tail-rotor motor voltage | V |
| $u_v$ | Main-rotor motor voltage | V |

## Three MPC-LPV formulations

The user specifies angle references, but effective TRMS control also requires
appropriate rotor speeds. Those rotor-speed references are not directly known
from the angle commands, so each example uses a different strategy to obtain
them: nonlinear feedforward equations, virtual MPC actions, or an outer MPC
whose predicted actions become references for inner rotor controllers.

| Example | MPC decision inputs | Rotor-speed references | Main purpose |
| --- | --- | --- | --- |
| [Basic MPC](TRMS_basic/readme.md) | Motor voltages $u_h,u_v$ | Computed externally from the nonlinear equilibrium equations | Classical centralized MIMO MPC |
| [Virtual-reference MPC](TRMS_refMPC/readme.md) | Motor voltages and virtual references $[u_h,u_v,\omega_h^{ref},\omega_v^{ref}]^T$ | Optimized inside the MIMO MPC using a CHRONOS custom cost | Flexible modeling with virtual control actions |
| [Cascade MPC](TRMS_cascade_mpc/readme.md) | Outer MPC: rotor-speed references; inner MPCs: motor voltages | Complete outer-MPC sequences are passed to two SISO rotor MPCs | Predictive cascade with dedicated nonlinear actuator control |

### Basic MPC

A single six-state MIMO MPC computes the two motor voltages directly. Because
the full state is tracked, references are required for the angles, body rates,
and rotor speeds. The rotor-speed references are obtained from nonlinear TRMS
equilibrium equations in *compute_ref.m*.

### Virtual-reference MPC

A single six-state MIMO MPC still computes the physical motor voltages, but it
also introduces two virtual control actions representing the rotor-speed
references. A CHRONOS custom cost makes the rotor speeds track these virtual
inputs, so nonlinear feedforward equations are not required to generate the
rotor references.

### Cascade MPC

The outer four-state MIMO MPC models the TRMS body dynamics and uses the rotor
speeds as its control inputs. Two dedicated SISO MPCs then control the nonlinear
rotor dynamics and compute the motor voltages. Each inner MPC receives the full
rotor-reference sequence predicted by the outer MPC, making use of the complete
MPC horizon rather than only its first action.

## Folder contents

- *TRMS.m*: nonlinear TRMS plant used by all simulations.
- *TRMS_basic*: centralized MIMO MPC with externally computed rotor references.
- *TRMS_refMPC*: centralized MIMO MPC with optimized virtual rotor references.
- *TRMS_cascade_mpc*: outer MIMO MPC with two inner SISO rotor MPCs.

## Reference

[1] Rotondo, D., Nejjari, F., & Puig, V. (2013). Quasi-LPV modeling,
identification and control of a Twin Rotor MIMO System. *Control Engineering
Practice, 21*(6), 829-846.
