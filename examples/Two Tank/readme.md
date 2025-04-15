# Two Tank System Example

### Folder structure

In this folder you will find 3 files with the following objectives:

* *two_tank_init.m*: script to define the MPC problem using the CHRONOS init functions.
* *sim_two_tank_lpv.m*: script to simulate the Two Tank system in closed-loop using the CHRONOS mpc solver, at each iteration we use the CHRONOS update functions to adapt its internal Linear Parameter Varying model to the instantaneous water height level.
* *sim_two_tank_lti.m*: this script is identical to *sim_two_tank_lpv.m*, except that the CHRONOS update step is skipped. This lead to a stable MPC which however has a significant tracking offset due to the differences between the linearization point of the internal MPC model and the state of "real" non-linear system.

### Example introduction

The two tanks system is a classical example used in control lectures. It consists of two inter-connected cylindrical tanks, where the first tank receives a steady flow of water, which can be regulated, and the second tank has a downstream discharge of water at its base. The dynamics of the two tank system is the following:

$$ \dot h_1 = u/A_b-\sqrt{2gh_1}/A_b $$
$$ \dot h2 = \sqrt{2gh_1}/A_b-\sqrt{2gh_2}/A_b $$

where $h_1$ and $h_2$ are the water heights of each tank, $u$ is the controlled water massflow into tank 1 and $A_b$ is the tank area, equal for both tanks. The regulation objective is to control the water level on the second tank, e.g. our tracking target is:

$$ y = h_2$$

### From non-linear to linear time varying system description

The CHRONOS solver is designed to solve Nonlinear MPC by making use of Linear Parameter Varying (LPV) models. In order to transform the non-linear system into a LPV one, we do a linear embeddeding by dividing and multiplying the square root terms by the respective tank height:

$$ \sqrt{2gh_i} := \frac{\sqrt{2gh_i}}{h_i}h_i  $$

Substituting the linear embeddings on the non-linear dynamics equation, we arrive at the following state-space LPV description of the Two Tank system dynamics:

$$ \left [\begin{array}{c} \dot h_1\\\ \dot h2 \end{array} \right ] =
\left [\begin{array}{cc}  -\sqrt{2gh_1}/(h_1A_b) & 0\\\ \sqrt{2gh_1}/(h_1A_b & -\sqrt{2gh_2}/(h_2A_b)\end{array} \right ]
\left [\begin{array}{c} h_1\\\  h2 \end{array} \right ] + 
\left [\begin{array}{c} 1/A_b\\\ 0 \end{array} \right ]u$$

Note that if we expand the state-space LPV model, we recover the exact nonlinear dynamics of the two-tank system. This highlights a key advantage of the LPV representation: it captures the full nonlinear behavior of the system while casting it in a form compatible with convex optimization. As a result, CHRONOS can solve the nonlinear MPC problem exactly, using fast, reliable, and well-established convex optimization algorithms.

### MPC Definition

In order to control the height of the second tank we solve at each iteration the following MPC problem using the CHRONOS solver:

$$\min_{u,x}J = (x_{ref}-x_N)^TP(x_{ref}-x_N) + \sum_{i=1}^{N-1} (r-y_i)^TQ_{e}(r-y_i) + \sum_{i=0}^{N_{ctr-1}}\Delta u_i^TdR_u\Delta u_i $$

s.t.

$$ x^+=A(h_1,h_2)x+Bu$$
$$ h_1,h_2 \in [0.01,\\\ 1]$$
$$ u \in [0,\\\ 10] \\\ \Delta u \in [-0.1,\\\ 0.1] $$

