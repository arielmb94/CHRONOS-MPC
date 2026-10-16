# CHRONOS

**Model Predictive Control for linear parameter-varying systems in MATLAB.**

CHRONOS lets you define a prediction model, objectives, and constraints in
MATLAB, then update the problem data and solve for the next control action at
each sampling instant. LPV models let you represent changing operating
conditions within a structured MPC formulation.

```{note}
This is an initial documentation preview. Tutorials, the user guide, and the
API reference will be added as the documentation develops.
```

## Controller workflow

1. **Define** the prediction horizon, dynamics, outputs, costs, and constraints.
2. **Build** the MPC structure and fixed-size solver workspace once.
3. **Update** model matrices and problem data as operating conditions change.
4. **Solve** using the current state, previous control action, and references.
5. **Apply** the first control action and repeat at the next sample.

CHRONOS is developed in MATLAB, with an online solver designed for C/C++
generation using MATLAB Coder.

## Explore CHRONOS

- [Source code and installation overview](https://github.com/arielmb94/CHRONOS-MPC)
- [Introductory tutorial](https://github.com/arielmb94/CHRONOS-MPC/tree/docs/tutorials/01%20-%20First%20Step)
- [Application examples](https://github.com/arielmb94/CHRONOS-MPC/tree/docs/examples)
- [MIT license](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/LICENSE)

## Contact

For technical inquiries and collaboration, contact
[Ariel Medero Borrell](https://www.linkedin.com/in/ariel-medero-borrell).
