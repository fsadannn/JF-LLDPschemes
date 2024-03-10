# Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince for large systems of differential equations

## This repo hold the original source code used for the Local Linearization papers

This Matlab toolbox provides the <strong>Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince (LLDP)</strong> schemes LLDP1 and LLDP2 described in [1] for the integration of large systems of initial value problems.

[1] Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince for large systems of differential equations
    by F.S. Naranjo-Noda and J.C. Jimenez
## <strong>Jacobian-free HOLL</strong>

### [```LLDP1```](./llint/LLDP1.m) variable step-size Jacobian-free Locally Locally Linearized Runge-Kutta method of Dormand and Prince with variable order finite difference

### [```LLDP2```](./llint/LLDP2.m) variable step-size Jacobian-free Locally Locally Linearized Runge-Kutta method of Dormand and Prince with order 1 finite difference

<br/>

### <strong>Demos</strong>

[`run_examples_lldp_fj.m`](./demo/run_examples_lldp_fj.m) generates the Tables 3-5 of [1] illustrating the performance of the Jacobian-free Locally Locally Linearized Runge-Kutta method of Dormand and Prince schemes LLDP1 and LLDP2 in the integration of test examples.

