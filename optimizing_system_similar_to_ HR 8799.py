import rebound
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

def initialize_simulation(initial_e_1,initial_e_2, initial_e_3,initial_e_4,mean_anomaly_2, mean_anomaly_3,mean_anomaly_4, kappa_1,kappa_2,kappa_3, pomega_1,pomega_2, pomega_3, pomega_4):
  sim = rebound.Simulation()
  a_1= 1.389481
  a_2= 2.20566405
  a_3= 4.20566405
  a_4= 8.20566405
  mean_anomaly_1 = 0
  sim.add(m=2)
  sim.add(m=3e-6, a=a_1, e=initial_e_1, M=mean_anomaly_1, pomega=pomega_1)
  sim.add(m=5e-6, a=a_2, e=initial_e_2, M=mean_anomaly_2, pomega=pomega_2)
  sim.add(m=7e-6, a=a_3, e=initial_e_3, M=mean_anomaly_3, pomega=pomega_3)
  sim.add(m=9e-6, a=a_4, e=initial_e_4, M=mean_anomaly_4, pomega=pomega_4)
  sim.move_to_com()
  kappa_1= sim.particles[1].P/sim.particles[2].P
  kappa_2= sim.particles[2].P/sim.particles[3].P
  kappa_3= sim.particles[3].P/sim.particles[4].P
  return [initial_e_1,initial_e_2,initial_e_3,initial_e_4, mean_anomaly_2, mean_anomaly_3, mean_anomaly_4, kappa_1,kappa_2,kappa_3, pomega_1,pomega_2, pomega_3,pomega_4]



def end_simulation(initial_e_1,initial_e_2, initial_e_3,initial_e_4,mean_anomaly_2, mean_anomaly_3,mean_anomaly_4, kappa_1,kappa_2,kappa_3, pomega_1,pomega_2, pomega_3, pomega_4):
  sim = rebound.Simulation()
  a_1= 1.389481
  a_2= 2.20566405
  a_3= 4.20566405
  a_4= 8.20566405
  mean_anomaly_1 = 0
  sim.add(m=2)
  sim.add(m=3e-6, a=a_1, e=initial_e_1, M=mean_anomaly_1, pomega=pomega_1)
  sim.add(m=5e-6, a=a_2, e=initial_e_2, M=mean_anomaly_2, pomega=pomega_2)
  sim.add(m=7e-6, a=a_3, e=initial_e_3, M=mean_anomaly_3, pomega=pomega_3)
  sim.add(m=9e-6, a=a_4, e=initial_e_4, M=mean_anomaly_4, pomega=pomega_4)
  sim.move_to_com()
  kappa_1= sim.particles[1].P/sim.particles[2].P
  kappa_2= sim.particles[2].P/sim.particles[3].P
  kappa_3= sim.particles[3].P/sim.particles[4].P

  integrate_to_mean_anomaly_16pi(sim, 1)

  sim.move_to_com()

  par=[sim.particles[1].e,sim.particles[2].e,sim.particles[3].e,sim.particles[4].e, sim.particles[2].M,sim.particles[3].M,sim.particles[4].M, sim.particles[1].P/sim.particles[2].P, sim.particles[2].P/sim.particles[3].P,sim.particles[3].P/sim.particles[4].P, sim.particles[1].pomega,sim.particles[2].pomega, sim.particles[3].pomega,  sim.particles[4].pomega]


  return par


def integrate_to_mean_anomaly_16pi(sim, body_index):
    target_mean_anomaly = 16 * np.pi
    initial_time_step = 1e-4
    time_step = initial_time_step

    while True:
        sim.integrate(sim.t + time_step)

        orbit = sim.particles[body_index].P
        current_mean_anomaly = sim.particles[1].M

        if np.abs(current_mean_anomaly  - 2*np.pi) < 1e-4:
          current_mean_anomaly = sim.particles[1].M +2*np.pi
          if np.abs(current_mean_anomaly  - 4*np.pi) < 1e-4:
            current_mean_anomaly = sim.particles[1].M +4*np.pi
            if np.abs(current_mean_anomaly  - 6*np.pi) < 1e-4:
              current_mean_anomaly = sim.particles[1].M +6*np.pi
              if np.abs(current_mean_anomaly  - 8*np.pi) < 1e-4:
                current_mean_anomaly = sim.particles[1].M +8*np.pi
                if np.abs(current_mean_anomaly  - 10*np.pi) < 1e-4:
                  current_mean_anomaly = sim.particles[1].M +10*np.pi
                  if np.abs(current_mean_anomaly  - 12*np.pi) < 1e-4:
                    current_mean_anomaly = sim.particles[1].M +12*np.pi
                    if np.abs(current_mean_anomaly  - 14*np.pi) < 1e-4:
                      current_mean_anomaly = sim.particles[1].M +14*np.pi


                      if np.abs((current_mean_anomaly  - target_mean_anomaly)) <1e-4:
                        break


def diff(theta):
  initial_e_1,initial_e_2, initial_e_3,initial_e_4,mean_anomaly_2, mean_anomaly_3,mean_anomaly_4, kappa_1,kappa_2,kappa_3, pomega_1,pomega_2, pomega_3, pomega_4 =theta
  ini_param = initialize_simulation(initial_e_1,initial_e_2, initial_e_3,initial_e_4,mean_anomaly_2, mean_anomaly_3,mean_anomaly_4, kappa_1,kappa_2,kappa_3, pomega_1,pomega_2, pomega_3, pomega_4)
  end_param = end_simulation(initial_e_1,initial_e_2, initial_e_3,initial_e_4,mean_anomaly_2, mean_anomaly_3,mean_anomaly_4, kappa_1,kappa_2,kappa_3, pomega_1,pomega_2, pomega_3, pomega_4)
  difference =[]
 
  for i in range(0,8):
    difference.append(ini_param[i]-end_param[i])
  difference.append(np.abs(ini_param[10]-ini_param[11])-np.abs(end_param[10]-end_param[11]))
  difference.append(np.abs(ini_param[11]-ini_param[12])-np.abs(end_param[11]-end_param[12]))
  difference.append(np.abs(ini_param[12]-ini_param[13])-np.abs(end_param[12]-end_param[13]))
  
  return np.sqrt(difference[0]**2+difference[1]**2+difference[2]**2+difference[3]**2+difference[4]**2+difference[5]**2+difference[6]**2+difference[7]**2+difference[8]**2+difference[9]**2+difference[10]**2)



ini_guess=[ 0.05,0.7,0.5,0.5, 0 ,0.5,0.5, 1, 0 ,0.5,0.5, 1,1,1]
bounds = [ (0, 0.99),(0, 0.99),(0, 0.99),(0, 0.99), (None, None), (None, None), (None, None), (None, None), (None, None), (None, None), (None, None), (None, None), (None, None), (None, None)]
result = minimize(diff, ini_guess, bounds=bounds)


initial_e_1,initial_e_2, initial_e_3,initial_e_4,mean_anomaly_2, mean_anomaly_3,mean_anomaly_4, kappa_1,kappa_2,kappa_3, pomega_1,pomega_2, pomega_3, pomega_4 =(result.x)
sim = rebound.Simulation()
a_1= 1.389481
a_2= 2.20566405
a_3= 4.20566405
a_4= 8.20566405
mean_anomaly_1 = 0
sim.add(m=2)
sim.add(m=3e-6, a=a_1, e=initial_e_1, M=mean_anomaly_1, pomega=pomega_1)
sim.add(m=5e-6, a=a_2, e=initial_e_2, M=mean_anomaly_2, pomega=pomega_2)
sim.add(m=7e-6, a=a_3, e=initial_e_3, M=mean_anomaly_3, pomega=pomega_3)
sim.add(m=9e-6, a=a_4, e=initial_e_4, M=mean_anomaly_4, pomega=pomega_4)
sim.move_to_com()
rebound.OrbitPlot(sim)


initial_e_1,initial_e_2, initial_e_3,initial_e_4,mean_anomaly_2, mean_anomaly_3,mean_anomaly_4, kappa_1,kappa_2,kappa_3, pomega_1,pomega_2, pomega_3, pomega_4 =(result.x)
sim = rebound.Simulation()
a_1= 1.389481
a_2= 2.20566405
a_3= 4.20566405
a_4= 8.20566405
mean_anomaly_1 = 0
sim.add(m=2)
sim.add(m=3e-6, a=a_1, e=initial_e_1, M=mean_anomaly_1, pomega=pomega_1)
sim.add(m=5e-6, a=a_2, e=initial_e_2, M=mean_anomaly_2, pomega=pomega_2)
sim.add(m=7e-6, a=a_3, e=initial_e_3, M=mean_anomaly_3, pomega=pomega_3)
sim.add(m=9e-6, a=a_4, e=initial_e_4, M=mean_anomaly_4, pomega=pomega_4)


target_mean_anomaly = 16 * np.pi
initial_time_step = 1e-4
time_step = initial_time_step

while True:
  sim.integrate(sim.t + time_step)
  current_mean_anomaly = sim.particles[1].M

  if np.abs(current_mean_anomaly  - 2*np.pi) < 1e-4:
    current_mean_anomaly = sim.particles[1].M +2*np.pi
    if np.abs(current_mean_anomaly  - 4*np.pi) < 1e-4:
      current_mean_anomaly = sim.particles[1].M +4*np.pi
      if np.abs(current_mean_anomaly  - 6*np.pi) < 1e-4:
        current_mean_anomaly = sim.particles[1].M +6*np.pi
        if np.abs(current_mean_anomaly  - 8*np.pi) < 1e-4:
          current_mean_anomaly = sim.particles[1].M +8*np.pi
          if np.abs(current_mean_anomaly  - 10*np.pi) < 1e-4:
            current_mean_anomaly = sim.particles[1].M +10*np.pi
            if np.abs(current_mean_anomaly  - 12*np.pi) < 1e-4:
              current_mean_anomaly = sim.particles[1].M +12*np.pi
              if np.abs(current_mean_anomaly  - 14*np.pi) < 1e-4:
                current_mean_anomaly = sim.particles[1].M +14*np.pi

                if np.abs((current_mean_anomaly  - target_mean_anomaly)) <1e-4:
                  break
sim.move_to_com()
rebound.OrbitPlot(sim)
