# `pc19` – Adaptive Variable-Order Predictor-Corrector ODE Solver

This MATLAB package implements **`pc19`**, a variable-order, variable-step solver for systems of non-stiff ordinary differential equations (ODEs). It is based on a **predictor-corrector strategy** using the **Adams-Bashforth-Moulton method**, with automatic step-size and order control up to 9th order.

The solver is optimized for high accuracy and efficiency by adapting both the order of the method and the time step size based on local error estimation.

---

## File Overview

| File | Description |
|------|-------------|
| `pc19.m` | Main solver routine implementing the predictor-corrector algorithm |
| `OneStep.m` | Initializes integration using Heun's method and estimates initial step size |
| `estimateNextH.m` | Newton’s method-based predictor for the next optimal step size |
| `choosefunc.m` | Defines error prediction and derivative functions up to 9th order |
| `getSteps.m` | Computes previous time step differences for the Adams method |
| `getParams2.m` to `getParams9.m` | Adams-Bashforth-Moulton coefficients and error estimators for orders 2–9 |

---

## Usage

### Basic Syntax

```matlab
[u, t, counter] = pc19(f, t0, y0, tf, tol);
[u, t, counter] = pc19(f, t0, y0, tf, tol, maxord);
```
### Input Arguments

| Argument | Type     | Default | Description                                                                 |
|----------|----------|---------|-----------------------------------------------------------------------------|
| `f`      | function handle | **required** | Function handle representing the ODE system, i.e., `f(t, y)` |
| `t0`     | float    | **required** | Initial time                                                              |
| `y0`     | vector   | **required** | Initial condition(s) for the ODE system                                   |
| `tf`     | float    | **required** | Final time                                                                |
| `tol`    | float    | **required** | Global error tolerance                                                    |
| `maxord` | integer  | `9`     | Maximum order of the method to be used (allowed: 2–9)                      |

---

### Output Values

| Output     | Type     | Description                                                                 |
|------------|----------|-----------------------------------------------------------------------------|
| `u`        | matrix   | Solution values; each row corresponds to the solution vector at a time step |
| `t`        | vector   | Time steps corresponding to the solution                                    |
| `counter`  | integer  | Number of function evaluations performed during integration                 |

---

## Examples
Two test cases are included. To run both test cases, simply run the full script containing the example calls to `pc19` and the visualizations.
### Three-Body Orbital Problem
Simulates a celestial mechanics system with gravitational interactions. Compares `pc19` to MATLAB’s `ode45` and `ode113`.
### Highly Oscillatory "Wiggly" Function
Solves a nonlinear ODE with a known exact solution to test solver accuracy.

---

## Acknowledgements
Thank you to Professor Mark Gockenbach for the invaluable guidance in preparing this algorithm. If you find this code useful in your research or publication, please cite this repository as appropriate.
