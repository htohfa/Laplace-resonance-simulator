import rebound
import numpy as np
from scipy.optimize import minimize

def initialize_simulation(a_1,a_2, initial_e_1,initial_e_2, mean_anomaly_1,mean_anomaly_2, pomega_1,pomega_2):
  sim = rebound.Simulation()
  sim.add(m=2)
  sim.add(m=3e-6, a=a_1, e=initial_e_1, M=mean_anomaly_1, pomega=pomega_1)
  sim.add(m=5e-6, a=a_2, e=initial_e_2, M=mean_anomaly_2, pomega=pomega_2)
  sim.move_to_com()

  return [a_1,a_2, initial_e_1,initial_e_2, mean_anomaly_1,mean_anomaly_2, pomega_1,pomega_2]


def end_simulation(a_1,a_2, initial_e_1,initial_e_2, mean_anomaly_1,mean_anomaly_2, pomega_1,pomega_2):
  sim = rebound.Simulation()
  sim.add(m=2)
  sim.add(m=3e-6, a=a_1, e=initial_e_1, M=mean_anomaly_1, pomega=pomega_1)
  sim.add(m=5e-6, a=a_2, e=initial_e_2, M=mean_anomaly_2, pomega=pomega_2)
  P_inner = 2 * np.pi * np.sqrt(sim.particles[1].a**3 / sim.G / (sim.particles[0].m + sim.particles[1].m))

  sim.integrate(2 * P_inner)
  sim.move_to_com()

  par=[sim.particles[1].a, sim.particles[2].a, sim.particles[1].e,sim.particles[2].e, sim.particles[1].M, sim.particles[2].M, sim.particles[1].pomega,sim.particles[2].pomega]


  return par

def diff(theta):
  a_1,a_2, initial_e_1,initial_e_2, mean_anomaly_1,mean_anomaly_2, pomega_1,pomega_2 =theta
  ini_param = initialize_simulation(a_1,a_2, initial_e_1,initial_e_2, mean_anomaly_1,mean_anomaly_2, pomega_1,pomega_2)
  end_param = end_simulation(a_1,a_2, initial_e_1,initial_e_2, mean_anomaly_1,mean_anomaly_2, pomega_1,pomega_2)
  difference =[]
  for i in range(0,6):
    difference.append(ini_param[i]-end_param[i])
  difference.append(np.abs(ini_param[6]-ini_param[7])-np.abs(end_param[6]-end_param[7]))

  return np.sqrt(difference[0]**2+difference[1]**2+difference[2]**2+difference[3]**2+difference[4]**2+difference[5]**2)


ini_guess=[ (2**(2/3)),2, 0.05,0.7, 0 ,0.5,0.5, 1]
bounds = [(0, None),(0, None), (0, 0.99),(0, 0.99), (None, None), (None, None), (None, None), (None, None)]
result = minimize(diff, ini_guess, bounds=bounds)

a_1,a_2, initial_e_1,initial_e_2, mean_anomaly_1,mean_anomaly_2, pomega_1,pomega_2 =(result.x)
sim = rebound.Simulation()
sim.add(m=2)
sim.add(m=3e-6, a=a_1, e=initial_e_1, M=mean_anomaly_1, pomega=pomega_1)
sim.add(m=5e-6, a=a_2, e=initial_e_2, M=mean_anomaly_2, pomega=pomega_2)
sim.move_to_com()
op= rebound.OrbitPlot(sim)
op.fig.savefig("initial_sim.pdf")


sim = rebound.Simulation()
sim.add(m=2)
sim.add(m=3e-6, a=a_1, e=initial_e_1, M=mean_anomaly_1, pomega=pomega_1)
sim.add(m=5e-6, a=a_2, e=initial_e_2, M=mean_anomaly_2, pomega=pomega_2)
P_inner = 2 * np.pi * np.sqrt(sim.particles[1].a**3 / sim.G / (sim.particles[0].m + sim.particles[1].m))
sim.integrate(2 * P_inner)
op = sim.move_to_com()
rebound.OrbitPlot(sim)
