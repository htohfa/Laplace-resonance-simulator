import rebound
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


def initialize_simulation(initial_e_1,initial_e_2,mean_anomaly_2, kappa, pomega_1,pomega_2):
  sim = rebound.Simulation()
  a_1= 1.389481
  a_2= 2.20566405
  mean_anomaly_1 = 0
  sim.add(m=2)
  sim.add(m=3e-6, a=a_1, e=initial_e_1, M=mean_anomaly_1, pomega=pomega_1)
  sim.add(m=5e-6, a=a_2, e=initial_e_2, M=mean_anomaly_2, pomega=pomega_2)
  sim.move_to_com()
  kappa= sim.particles[2].P/sim.particles[1].P

  return [initial_e_1,initial_e_2, mean_anomaly_2, kappa, pomega_1,pomega_2]


def end_simulation(initial_e_1,initial_e_2, mean_anomaly_2, kappa, pomega_1,pomega_2):
  sim = rebound.Simulation()
  a_1= 1.389481
  a_2= 2.20566405
  mean_anomaly_1 = 0
  sim.add(m=2)
  sim.add(m=3e-6, a=a_1, e=initial_e_1, M=mean_anomaly_1, pomega=pomega_1)
  sim.add(m=5e-6, a=a_2, e=initial_e_2, M=mean_anomaly_2, pomega=pomega_2)
  integrate_to_mean_anomaly_2pi(sim, 1)

  sim.move_to_com()

  par=[sim.particles[1].e,sim.particles[2].e, sim.particles[2].M, sim.particles[2].P/sim.particles[1].P, sim.particles[1].pomega,sim.particles[2].pomega]


  return par


def integrate_to_mean_anomaly_2pi(sim, body_index):
    target_mean_anomaly = 4 * np.pi
    initial_time_step = 1e-4
    time_step = initial_time_step

    while True:
        sim.integrate(sim.t + time_step)

        orbit = sim.particles[body_index].P
        current_mean_anomaly = sim.particles[1].M

        if np.abs(current_mean_anomaly  - 2*np.pi) < 1e-4:
          current_mean_anomaly = sim.particles[1].M +2*np.pi
          if np.abs((current_mean_anomaly  - target_mean_anomaly)) <1e-4:
            break


def diff(theta):
  initial_e_1,initial_e_2, mean_anomaly_2,kappa, pomega_1,pomega_2 =theta
  ini_param = initialize_simulation(initial_e_1,initial_e_2, mean_anomaly_2, kappa, pomega_1,pomega_2)
  end_param = end_simulation(initial_e_1,initial_e_2, mean_anomaly_2, kappa, pomega_1,pomega_2)
  difference =[]
  for i in range(0,4):
    difference.append(ini_param[i]-end_param[i])
  difference.append(np.abs(ini_param[4]-ini_param[5])-np.abs(end_param[4]-end_param[5]))

  return np.sqrt(difference[0]**2+difference[1]**2+difference[2]**2+difference[3]**2+difference[4]**2)

ini_guess=[ 0.05,0.7, 0 ,0.5,0.5, 1]
bounds = [ (0, 0.99),(0, 0.99), (None, None), (None, None), (None, None), (None, None)]
result = minimize(diff, ini_guess, bounds=bounds)


initial_e_1,initial_e_2, mean_anomaly_2, kappa, pomega_1,pomega_2 = (result.x)


sim = rebound.Simulation()
sim.add(m=2)
sim.add(m=3e-6, a= 1.389481, e=initial_e_1, M=0, pomega=pomega_1)
sim.add(m=5e-6, a=  2.20566405, e=initial_e_2, M=mean_anomaly_2, pomega=pomega_2)

sim.move_to_com()
rebound.OrbitPlot(sim)

initial_e_1,initial_e_2, mean_anomaly_2, kappa, pomega_1,pomega_2= (result.x)
sim = rebound.Simulation()
sim.add(m=2)
sim.add(m=3e-6, a=1.389481, e=initial_e_1, M=0, pomega=pomega_1)
sim.add(m=5e-6, a= 2.20566405, e=initial_e_2, M=mean_anomaly_2, pomega=pomega_2)


target_mean_anomaly = 4 * np.pi
initial_time_step = 1e-4
time_step = initial_time_step

while True:
  sim.integrate(sim.t + time_step)
  orbit = sim.particles[1].P
  current_mean_anomaly = sim.particles[1].M

  if np.abs(current_mean_anomaly  - 2*np.pi) < 1e-4:
      current_mean_anomaly = sim.particles[1].M +2*np.pi
      if np.abs((current_mean_anomaly  - target_mean_anomaly)) <1e-4:
        break

rebound.OrbitPlot(sim)
