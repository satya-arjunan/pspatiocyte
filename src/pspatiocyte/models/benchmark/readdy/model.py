import os
import numpy as np
import readdy
import sys
import math
import time
import csv

dt = 1e-3 # s
volume = 909.0
box_length = math.pow(volume,1./3.)
n_particles_e = 9090
n_particles_s = 90910
duration = 10. # s

system = readdy.ReactionDiffusionSystem(
    box_size=[box_length, box_length, box_length], # for volume 86.15 um^3
    unit_system={"length_unit": "micrometer", "time_unit": "second"}
)

system.add_species("E", diffusion_constant=10.)
system.add_species("S", diffusion_constant=10.)
system.add_species("ES", diffusion_constant=10.)
system.add_species("P", diffusion_constant=10.)

system.reactions.add("fwd: E +(0.03) S -> ES", rate=86.78638438)
system.reactions.add("back: ES -> E +(0.03) S", rate=1.)
system.reactions.add("prod: ES -> E +(0.03) P", rate=1.)

simulation = system.simulation(kernel="CPU")
simulation.kernel_configuration.n_threads = int(sys.argv[1])
simulation.output_file = "out.h5"
simulation.reaction_handler = "UncontrolledApproximation"

edge_length = system.box_size[0]
initial_positions_e = np.random.random(size=(n_particles_e, 3)) * edge_length - .5*edge_length
initial_positions_s = np.random.random(size=(n_particles_s, 3)) * edge_length - .5*edge_length
simulation.add_particles("E", initial_positions_e)
simulation.add_particles("S", initial_positions_s)

simulation.observe.number_of_particles(stride=10, types=["E", "S", "ES", "P"])

if os.path.exists(simulation.output_file):
  os.remove(simulation.output_file)

n_steps = int(duration / dt)

start = time.time()
simulation.run(n_steps=n_steps, timestep=dt)
end = time.time()

traj = readdy.Trajectory(simulation.output_file)
times, counts = traj.read_observable_number_of_particles()
times = times * dt

rows = zip(times, counts[:,0], counts[:,1], counts[:,2], counts[:,3])
with open("output.txt", "w+") as f:
  writer = csv.writer(f)
  writer.writerow(['time', 'E', 'S', 'ES', 'P'])
  for row in rows:
    writer.writerow(row)

print(end-start)
