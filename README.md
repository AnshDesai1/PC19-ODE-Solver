# `pc19` – Adaptive Variable-Order Predictor-Corrector ODE Solver

This MATLAB package implements **`pc19`**, a variable-order, variable-step solver for systems of non-stiff ordinary differential equations (ODEs). It is based on a **predictor-corrector strategy** using the **Adams-Bashforth-Moulton method**, with automatic step-size and order control up to 9th order.

The solver is optimized for high accuracy and efficiency by adapting both the order of the method and the time step size based on local error estimation.

---

## Features

- Variable-order Adams methods (orders 2 through 9)
- Adaptive time-stepping based on local truncation error (LTE)
- Initial step taken with **Heun's Method**
- Local extrapolation using PLTE (Predictor Local Truncation Error)
- Flexible for both scalar and system ODEs
- Compatible with MATLAB’s `f(t, y)` function handle format

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
