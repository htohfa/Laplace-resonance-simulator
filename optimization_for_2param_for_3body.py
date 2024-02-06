import rebound
import numpy as np
sim = rebound.Simulation()
sim.add(m=1)
sim.add(e=0.4, a=0.4, inc=0.2, f=0.43, Omega=0.82, omega=2.98)
sim.add(e=0.2, a=1.0, pomega=2.14)
sim.add(e=0.2, a=1.5, omega=1.14)
op = rebound.OrbitPlot(sim)


def initialize_simulation(semi_major_axes, initial_eccentricities):
    sim = rebound.Simulation()
    sim.add(m=2)
    for i in range(2):
        sim.add(m=3e-5, a=semi_major_axes[i], e=initial_eccentricities[i])

    sim.move_to_com()
    return sim

def objective_function(sim, initial_eccentricities):
    P_inner = 2 * np.pi * np.sqrt(sim.particles[1].a**3 / sim.G / (sim.particles[0].m + sim.particles[1].m))
    sim.integrate(2 * P_inner, exact_finish_time=0)
    de = [sim.particles[i+1].e - initial_eccentricities[i] for i in range(2)]
    return np.array(de)

def newtons_method(sim, initial_guess, tol=1e-5, max_iter=100):
    x = initial_guess
    for i in range(max_iter):
        sim_copy = sim.copy()
        set_eccentricities(sim_copy, x)
        F = objective_function(sim_copy, x)
        if np.linalg.norm(F) < tol:
            break

        J = np.zeros((len(x), len(x)))
        for j in range(len(x)):
            dx = np.zeros_like(x)
            dx[j] = x[j] * 1e-5 + 1e-8
            sim_copy = sim.copy()
            set_eccentricities(sim_copy, x + dx)
            F_dx = objective_function(sim_copy, x + dx)
            J[:, j] = (F_dx - F) / dx[j]

        delta_x = np.linalg.solve(J, -F)

        scale_factor = 0.1
        new_x = x + scale_factor * delta_x


        new_x = np.clip(new_x, 0, 1)

        if np.linalg.norm(new_x - x) < tol:
            break

        x = new_x

    return x



def set_eccentricities(sim, eccentricities):
    for i, e in enumerate(eccentricities):
        sim.particles[i+1].e = e

semi_major_axes = [1.0, (2**(2/3))]

initial_eccentricities = [0.1, 0.01]
sim = initialize_simulation(semi_major_axes, initial_eccentricities)

optimized_eccentricities = newtons_method(sim, initial_eccentricities)
optimized_eccentricities


sim = rebound.Simulation()

sim.add(m=2)

sim.add(m=6e-5, a=1.0, e=0.01)
sim.add(m=4e-5, a=(2**(2/3)), e=0.01)

rebound.OrbitPlot(sim)
