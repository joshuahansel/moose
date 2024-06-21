#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from thm_utilities import readCSVFile

runs = [
  ('main_2d_out.csv', '2D', 'black', '-'),
  ('main_3d_out.csv', '3D', 'red', '--')
]

data = dict()
for run_filename, run_label, run_color, run_ls in runs:
  data[run_filename] = readCSVFile(run_filename)
  data[run_filename]['t_index'] = range(len(data[run_filename]['time']))
  data[run_filename]['fp_it_per_time'] = data[run_filename]['n_fp_iterations'] / data[run_filename]['dt']
  data[run_filename]['cum_fp_it'] = np.cumsum(data[run_filename]['n_fp_iterations'])

def initializePlot():
  plt.figure(figsize=(8, 6))
  plt.rc('text', usetex=True)
  plt.rc('font', family='sans-serif')
  ax = plt.subplot(1, 1, 1)
  ax.get_yaxis().get_major_formatter().set_useOffset(False)

def makeLinearPlot(x_name, y_name):
  for run_filename, run_label, run_color, run_ls in runs:
    plt.plot(data[run_filename][x_name], data[run_filename][y_name], linestyle=run_ls, color=run_color, marker='', label=run_label)

def makeLogPlot(x_name, y_name):
  for run_filename, run_label, run_color, run_ls in runs:
    plt.semilogy(data[run_filename][x_name], data[run_filename][y_name], linestyle=run_ls, color=run_color, marker='', label=run_label)

def finalizePlot(plot_filename):
  plt.legend(frameon=False, prop={'size':10})
  plt.tight_layout()
  plt.savefig(plot_filename, dpi=300)

# Time step size

initializePlot()
plt.xlabel("Time Step Index")
plt.ylabel("Time Step Size [s]")
makeLogPlot('t_index', 'dt')
finalizePlot('dt_by_index.png')

# Steady-state norm

initializePlot()
plt.xlabel("Time [s]")
plt.ylabel("Steady-State Norm")
makeLogPlot('time', 'ss_err')
finalizePlot('ss_err_by_time.png')

# Fixed point iterations

initializePlot()
plt.xlabel("Time [s]")
plt.ylabel("FP Iterations")
makeLinearPlot('time', 'n_fp_iterations')
finalizePlot('fp_it_by_time.png')

# Fixed point iterations per time

initializePlot()
plt.xlabel("Time Step Index")
plt.ylabel("FP Iterations per Time [1/s]")
makeLogPlot('t_index', 'fp_it_per_time')
finalizePlot('fp_it_per_time_by_index.png')

# Cumulative fixed point iterations

initializePlot()
plt.xlabel("Time [s]")
plt.ylabel("Cumulative FP Iterations")
makeLinearPlot('time', 'cum_fp_it')
finalizePlot('cum_fp_it_by_time.png')

# Wall time

initializePlot()
plt.xlabel("Time [s]")
plt.ylabel("Wall Time [s]")
makeLinearPlot('time', 'wall_time')
finalizePlot('wall_time.png')
