# sn2 — Numerical Methods in MATLAB

MATLAB scripts developed during a numerical methods course. Organized by topic for quick access during classes and exams.

## Structure

```
sn2/
├── root_finding/          # Methods for solving f(x) = 0
│   ├── bisection.m        # Bisection method
│   ├── newton_raphson.m   # Newton–Raphson method
│   └── secant.m           # Secant method
├── integration/           # Numerical quadrature
│   ├── trapezoidal.m      # Composite trapezoidal rule
│   └── simpsons.m         # Composite Simpson's 1/3 rule
├── differentiation/       # Finite-difference derivatives
│   └── finite_diff.m      # Forward, backward, and central differences
├── linear_systems/        # Ax = b solvers
│   ├── gauss_elimination.m  # Gaussian elimination with partial pivoting
│   └── lu_decomposition.m   # LU decomposition (Doolittle)
├── interpolation/         # Polynomial interpolation
│   ├── lagrange.m         # Lagrange interpolating polynomial
│   └── newton_divided_diff.m  # Newton's divided-difference interpolation
└── ode_solvers/           # Ordinary differential equations
    ├── euler.m            # Explicit Euler method
    └── runge_kutta4.m     # Classical 4th-order Runge–Kutta
```

## Quick Usage

Every function follows a consistent calling convention. Examples:

```matlab
% Root finding
root = bisection(@(x) x^3 - x - 2, 1, 2, 1e-6, 100);

% Numerical integration
I = trapezoidal(@(x) sin(x), 0, pi, 100);

% ODE solving
[t, y] = runge_kutta4(@(t,y) -y, 0, 1, 0.1, 1);
```

See the header comment in each `.m` file for the full function signature and a worked example.

## Topics Covered

| Topic | Methods |
|---|---|
| Root Finding | Bisection, Newton–Raphson, Secant |
| Integration | Trapezoidal rule, Simpson's 1/3 rule |
| Differentiation | Forward / Backward / Central differences |
| Linear Systems | Gaussian elimination, LU decomposition |
| Interpolation | Lagrange, Newton divided differences |
| ODEs | Euler, Runge–Kutta 4th order |
