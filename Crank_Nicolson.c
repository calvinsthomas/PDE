#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1000
#define T_MAX 0.1
#define DX (1.0 / (N - 1))
#define DT 0.00001  // Crank-Nicolson unconditionally stable, no CFL.
#define STEPS (int)(T_MAX / DT)

void solve_pde(double *u, int n, int steps, double dx, double dt) {
    double alpha = dt / (dx * dx);
    double *a = (double*)malloc((n-2) * sizeof(double));
    double *b = (double*)malloc((n-2) * sizeof(double));
    double *c = (double*)malloc((n-2) * sizeof(double));
    double *d = (double*)malloc((n-2) * sizeof(double));
    double *u_next = (double*)malloc(n * sizeof(double));

    // Initialize the tridiagonal matrix coefficients
    for (int i = 0; i < n-2; i++) {
        a[i] = -alpha;
        b[i] = 1 + 2 * alpha;
        c[i] = -alpha;
    }

    // Time-stepping loop
    for (int t = 0; t < steps; t++) {
        // Set up the right-hand side vector d
        for (int i = 1; i < n-1; i++) {
            double u_xx = (u[i+1] - 2*u[i] + u[i-1]) / (dx * dx);
            double u_x = (u[i+1] - u[i-1]) / (2 * dx);
            d[i-1] = u[i] + dt * (u_xx - u[i] * u_x) / 2.0;
        }

        // Apply boundary conditions in the RHS vector
        d[0] += alpha * 1;  // u(0, t) = 1
        d[n-3] += alpha * exp(-1 - t * dt);  // u(1, t) = e^(-1-t)

        // Solve the tridiagonal system (using Thomas algorithm)
        for (int i = 1; i < n-2; i++) {
            double m = a[i] / b[i-1];
            b[i] -= m * c[i-1];
            d[i] -= m * d[i-1];
        }

        u_next[n-2] = d[n-3] / b[n-3];
        for (int i = n-4; i >= 0; i--) {
            u_next[i+1] = (d[i] - c[i] * u_next[i+2]) / b[i];
        }

        // Apply boundary conditions
        u_next[0] = 1;  // u(0, t) = 1
        u_next[n-1] = exp(-1 - t * dt);  // u(1, t) = e^(-1-t)

        // Swap pointers to avoid copying arrays
        double *temp = u;
        u = u_next;
        u_next = temp;
    }

    // Copy the result back to the original array
    memcpy(u, u_next, n * sizeof(double));

    // Free allocated memory
    free(a);
    free(b);
    free(c);
    free(d);
    free(u_next);
}

int main() {
    double *u = (double*)malloc(N * sizeof(double));
    if (!u) {
        perror("Failed to allocate memory for u");
        exit(EXIT_FAILURE);
    }

    // Initialize the solution array with the initial condition u(x, 0) = e^x
    for (int i = 0; i < N; i++) {
        u[i] = exp(i * DX);
    }

    // Solve the PDE using the Crank-Nicolson method
    solve_pde(u, N, STEPS, DX, DT);

    // Output the final solution (save to file)
    for (int i = 0; i < N; i++) {
        printf("u(%f, T_MAX) = %f\n", i * DX, u[i]);
    }

    // Free allocated memory
    free(u);

    return 0;
}
