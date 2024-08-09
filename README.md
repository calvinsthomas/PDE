# PDE
PDE - Boundary - .c

Latex Final Formulas:

1. Explicit Finite Difference Method:
u_i^{n+1} = u_i^n + \Delta t \left( \frac{u_{i+1}^n - 2u_i^n + u_{i-1}^n}{\Delta x^2} - u_i^n \cdot \frac{u_{i+1}^n - u_{i-1}^n}{2\Delta x} \right)

2. Crank-Nicolson Method:
\frac{u_i^{n+1} - u_i^n}{\Delta t} + \frac{u_i^n + u_i^{n+1}}{2} \cdot \frac{u_{i+1}^n - u_{i-1}^n + u_{i+1}^{n+1} - u_{i-1}^{n+1}}{4\Delta x} = \frac{1}{2} \left( \frac{u_{i+1}^n - 2u_i^n + u_{i-1}^n}{\Delta x^2} + \frac{u_{i+1}^{n+1} - 2u_i^{n+1} + u_{i-1}^{n+1}}{\Delta x^2} \right)

3. Forward-Time Centered-Space (FTCS) Scheme for Linear Diffusion:
u_i^{n+1} = u_i^n + \Delta t \cdot \frac{u_{i+1}^n - 2u_i^n + u_{i-1}^n}{\Delta x^2}

4. Forward Elimination:
Modified coefficients
b'_i = b_i - \frac{a_i c_{i-1}}{b'_{i-1}}, \quad d'_i = d_i - \frac{a_i d'_{i-1}}{b'_{i-1}}

5. Backward Substitution:
u_{n-1} = \frac{d'_{n-1}}{b'_{n-1}}, \quad u_i = \frac{d'_i - c_i u_{i+1}}{b'_i}, \quad \text{for } i = n-2, n-3, \ldots, 1
