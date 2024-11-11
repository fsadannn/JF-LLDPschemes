# Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince for large systems of differential equations

## This repo holds the original source code used for the Local Linearization papers

This Matlab toolbox provides the two implementations of the <strong> Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince (LLDP)</strong> described in [1] for the integration of large systems of initial value problems.

Furthermore, the toolbox provides the implementation of the Concurrent Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince (LLDP) described in [2] for the integration of the same type of equations.

[1] [Naranjo-Noda, F. S., & Jimenez, J. C. (2024). Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince for large systems of differential equations. Journal of Computational and Applied Mathematics, 449, 115974.](https://doi.org/10.1016/j.cam.2024.115974)

[2] New class of Jacobian-free high order local Linearization methods for differential equations by F.S. Naranjo-Noda and J.C. Jimenez

## <strong>Jacobian-free HOLL</strong>

### [```LLDP1```](./llint/LLDP1.m) variable step-size Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince with variable order finite difference in the Arnoldi algorithm

### [```LLDP2```](./llint/LLDP2.m) variable step-size Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince with fixed order 1 finite difference in the Arnoldi algorithm

### [```LLDP```](./llint/LLDP.m) variable step-size Concurrent Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince

<br/>

### <strong>Demos</strong>

[`run_examples_lldp_fj.m`](./demo/run_examples_lldp_fj.m) generates the Tables 3-5 of [1] illustrating the performance of the Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince schemes LLDP1 and LLDP2 in the integration of test examples.

[`run_large_examples_lldp_fj.m`](./demo/run_large_examples_lldp_fj.m) generates the Table 7 of [1] illustrating the performance of the Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince schemes LLDP1 and LLDP2 in the integration of test examples with larger dimensions.

[`run_examples_lldp_fj.m`](./demo_lldp/run_examples_lldp_fj_i.m) generates the Tables 1-2 of [2] illustrating the performance of the Concurrent Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince scheme LLDP in the integration of test examples.
