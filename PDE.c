#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>  // For parallelization

// Constants definition for problem parameters
#define N 1000  // Number of spatial points
#define T_MAX 1.0  // Maximum time
#define DX (1.0 / (N - 1))  // Spatial step size
#define DT 0.0001  // Time step size
#define STEPS (int)(T_MAX / DT)  // Number of time steps

/**
 * @brief Solves the PDE using an explicit finite difference method.
 * 
 * This function computes the solution to the PDE u_t + u u_x = u_xx
 * over the spatial domain [0,1] for a given time period. The method
 * uses a simple forward time-centered space scheme.
 *
 * @param u Pointer to the array holding the initial condition and
 *          which will be updated with the solution.
 * @param n Number of spatial points.
 * @param steps Number of time steps.
 * @param dx Spatial step size.
 * @param dt Time step size.
 */
void solve_pde(double *u, int n, int steps, double dx, double dt) {
    // Allocate memory for the next time step array
    double *u_next = (double*)malloc(n * sizeof(double));
    if (!u_next) {
        perror("Failed to allocate memory for u_next");
        exit(EXIT_FAILURE);
    }

    // Time-stepping loop
    for (int t = 0; t < steps; t++) {
        // Parallelize the spatial loop for better performance
        #pragma omp parallel for
        for (int i = 1; i < n-1; i++) {
            // Calculate second derivative u_xx and first derivative u_x
            double u_xx = (u[i-1] - 2*u[i] + u[i+1]) / (dx * dx);
            double u_x = (u[i+1] - u[i-1]) / (2 * dx);
            // Update the solution using the finite difference scheme
            u_next[i] = u[i] + dt * (u_xx - u[i] * u_x);
        }

        // Apply boundary conditions
        u_next[0] = 1;  // u(0, t) = 1
        u_next[n-1] = exp(-1 - t * dt);  // u(1, t) = e^(-1-t)

        // Swap pointers to avoid copying arrays, enhancing performance
        double *temp = u;
        u = u_next;
        u_next = temp;
    }

    // Copy the result back to the original array
    memcpy(u, u_next, n * sizeof(double));

    // Free allocated memory for the next time step array
    free(u_next);
}

/**
 * @brief The main function that initializes the problem and solves the PDE.
 * 
 * This function sets up the initial condition, calls the solve_pde
 * function, and then outputs the final solution.
 *
 * @return int Returns 0 on successful execution.
 */
int main() {
    // Allocate memory for the solution array
    double *u = (double*)malloc(N * sizeof(double));
    if (!u) {
        perror("Failed to allocate memory for u");
        exit(EXIT_FAILURE);
    }

    // Initialize the solution array with the initial condition u(x, 0) = e^x
    for (int i = 0; i < N; i++) {
        u[i] = exp(i * DX);
    }

    // Solve the PDE using the finite difference method
    solve_pde(u, N, STEPS, DX, DT);

    // Output the final solution (save to file)
    for (int i = 0; i < N; i++) {
        printf("u(%f, T_MAX) = %f\n", i * DX, u[i]);
    }

    // Free allocated memory for the solution array
    free(u);

    return 0;
}
