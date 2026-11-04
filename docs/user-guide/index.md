# User Guide

The User Guide introduces the MPC problem and ingredients supported by
CHRONOS, then explains how to translate nonlinear models into its LPV form
and configure, solve, and update the controller.

Read the chapters in order:

1. **CHRONOS overview:** identify the supported models, objectives, and constraints.
2. **From nonlinear to LPV models:** prepare the prediction and signal models.
3. **Configuring the MPC controller:** register the ingredients and build the controller.
4. **Solving and updating the MPC:** update problem data and compute a control action at each sample.
5. **Your first controller:** combine these steps in a two-tank simulation.

Start with [CHRONOS overview](chronos-overview.md).

```{toctree}
:maxdepth: 1

CHRONOS overview <chronos-overview>
from-nonlinear-problems-to-lpv-mpc
Configuring the MPC controller <configuring-the-mpc-controller>
Solving and updating the MPC <solving-and-updating-the-mpc>
Your first controller <your-first-controller>
```
