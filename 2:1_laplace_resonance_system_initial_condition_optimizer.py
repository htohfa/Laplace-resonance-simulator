import rebound
import numpy as np

def initialize_simulation(semi_major_axes, initial_eccentricities, mean_anomalies, longitude_of_pericenter):
    sim = rebound.Simulation()
    sim.add(m=2)  
    sim.add(m=3e-6, a=semi_major_axes[0], e=initial_eccentricities[0], M=mean_anomalies[0], pomega=longitude_of_pericenter[0])
    sim.add(m=5e-6, a=semi_major_axes[1], e=initial_eccentricities[1], M=mean_anomalies[1], pomega=longitude_of_pericenter[1])
    sim.move_to_com()
    return sim

def objective_function(sim, initial_parameters):
 
    P_inner = 2 * np.pi * np.sqrt(sim.particles[1].a**3 / sim.G / (sim.particles[0].m + sim.particles[1].m))
    P_outer = 2 * np.pi * np.sqrt(sim.particles[2].a**3 / sim.G / (sim.particles[0].m + sim.particles[2].m))
    current_period_ratio = P_outer / P_inner
    
    diffs = []
    for i in range(2):
        diffs.append(sim.particles[i+1].e - initial_parameters['eccentricities'][i])
        diffs.append(sim.particles[i+1].M - initial_parameters['mean_anomalies'][i])
        diffs.append(sim.particles[i+1].pomega - initial_parameters['longitude_of_pericenter'][i])
    
    diffs.append(current_period_ratio - initial_parameters['period_ratio'])
    
    return np.array(diffs)


def newtons_method(sim, initial_parameters, initial_semi_major_axes, tol=1e-12, max_iter=50):
    initial_f_outer = 1
    x = np.concatenate([
    initial_parameters['eccentricities'],
    initial_parameters['mean_anomalies'],
    initial_parameters['longitude_of_pericenter'],
    [initial_f_outer] 
])
    for i in range(max_iter):
        sim_copy = sim.copy()
        
        temp_parameters = {
            'eccentricities': x[:2],
            'mean_anomalies': x[2:4],
            'longitude_of_pericenter': x[4:6],
            'f_outer': x[6],
            'period_ratio': initial_parameters['period_ratio']
        }
        set_parameters(sim_copy, temp_parameters, initial_semi_major_axes) 
        F = objective_function(sim_copy, temp_parameters) 
        if np.linalg.norm(F) < tol:
            break

        J = np.zeros((len(x), len(x)))
        for j in range(len(x)):
            dx = np.zeros_like(x)
            dx[j] = x[j] * 1e-5 + 1e-8  # Perturbation
            sim_copy = sim.copy()
            
            temp_parameters_perturbed = {
                'eccentricities': x[:2] + dx[:2],
                'mean_anomalies': x[2:4] + dx[2:4],
                'longitude_of_pericenter': x[4:6] + dx[4:6],
                'f_outer': x[6] + dx[6], 
                'period_ratio': initial_parameters['period_ratio']
            }
            set_parameters(sim_copy, temp_parameters_perturbed, initial_semi_major_axes)  
            
            F_dx = objective_function(sim_copy, temp_parameters_perturbed)  
            J[:, j] = (F_dx - F) / dx[j]

        delta_x = np.linalg.solve(J, -F)

        scale_factor = 0.1
        new_x = x + scale_factor * delta_x
        new_x[:2] = np.clip(new_x[:2], 0, 1)

        if np.linalg.norm(new_x - x) < tol:
            break

        x = new_x

    optimized_parameters = {
        'eccentricities': x[:2],
        'mean_anomalies': x[2:4],
        'longitude_of_pericenter': x[4:6]
    }

    return optimized_parameters



def set_parameters(sim, parameters, initial_semi_major_axes):
    for i in range(2):
        sim.particles[i+1].e = parameters['eccentricities'][i]
        sim.particles[i+1].M = parameters['mean_anomalies'][i]
        sim.particles[i+1].pomega = parameters['longitude_of_pericenter'][i]
    f_outer = parameters['f_outer']
    sim.particles[2].a = initial_semi_major_axes[1] * f_outer

initial_eccentricities = [0.4, 0.5]
initial_mean_anomalies = [0.2, 0.5] 
initial_longitude_of_pericenter = [0.2, 1.0] 

initial_semi_major_axes = [2.0, (2**(2/3))]  
initial_parameters = {
    'eccentricities': initial_eccentricities,
    'mean_anomalies': initial_mean_anomalies,
    'longitude_of_pericenter': initial_longitude_of_pericenter,
    'period_ratio': 1.5  
}

sim = initialize_simulation(initial_semi_major_axes, initial_eccentricities, initial_mean_anomalies, initial_longitude_of_pericenter)
optimized_parameters = newtons_method(sim, initial_parameters, initial_semi_major_axes)
print(optimized_parameters)

sim = rebound.Simulation()

sim.add(m=2)

sim.add(m=3e-6, a=initial_semi_major_axes[0], M= optimized_parameters['mean_anomalies'][0], e=optimized_parameters['eccentricities'][0], pomega = optimized_parameters['longitude_of_pericenter'][0])
sim.add(m=4e-6, a= initial_semi_major_axes[1], M= optimized_parameters['mean_anomalies'][1], e=optimized_parameters['eccentricities'][1], pomega = optimized_parameters['longitude_of_pericenter'][1])


op = rebound.OrbitPlot(sim)
op.fig.savefig("orbit.pdf")
