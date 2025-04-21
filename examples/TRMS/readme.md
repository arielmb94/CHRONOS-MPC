# Twin Rotor MIMO System Examples

### Folder structure

The Twin Rotor MIMO System (TRMS) is a challenging mechanical system to control due its high nonlinearity. This makes it a perfect demonstrator of how CHRONOS can solve complex Nonlinear MPC problems thanks to the use of Linear Parameter Varying models. We present three different variations of the MPC control problem for the TRMS system in order to showcase the versatility of CHRONOS. On this folder you will find the following files and folders:

* *TRMS.m*: function implementing the nonlinear dynamics of the TRMS System, used for simulations. The TRMS model equations and parameters are borrowed from [1].
* *TRMS_basic*: implementation of a standard centralized MPC approach to solve the TRMS control problem. The MPC takes into account all states of the TRMS system and the state references are computed from desired setpoints and equilibrium equations.
* *TRMS_refMPC*: the MPC implementation in this case uses an extended control input vector, where the extended inputs act as references for some of the TRMS states. Compared with the basic MPC example, this implementation does not require the computation of state referenes based on the nonlinear system equilibrium equations as the MPC computes them internally. This aim to showcase an interesting use of CHRONOS's performance vectors $z$ to achieve personalized optimization behaviour. 
* *TRMS_cascade_mpc*: for this implementation we break down the TRMS control problem into three different MPC. A MIMO MPC defines the references and two independent SISO MPC controllers compute the final control actions. Different to other cascade controller strategies, the inner SISO MPC controllers do not take a single static reference point, instead they make use of the full reference sequence computed by the outer MIMO MPC.

### References

[1] Rotondo, D., Nejjari, F., & Puig, V. (2013). Quasi-LPV modeling, identification and control of a Twin Rotor MIMO System. Control Engineering Practice, 21(6), 829-846.

