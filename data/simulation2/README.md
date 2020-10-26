# Simulation data from the PPRZ logs
The original files are in the directory [simulation1](../simulation1).

The mission was the original agricultural survey with Paparazzi autopilot and Opterra fixed-wing UAV (credit goes to **Amit** Ferenz Appel and **Hector** Garcia de Marina).

Files here have been cut to the actual survey. You can find how to interpret the data here, so next time you don't forget; as you always do.

## Files
### [path_sim.csv](path_sim.csv). 
Data were sampled at 1 sec, we started at 680.732 and finished at 999.732. 

Contains 3 columns: 
1. the x related to origin. The [max,min] is x = [-123,85].
2. the y, y = [-66,221].
3. the z, z = [10.18,23.37].


### [max_qos_tees.csv](max_qos_tees.csv). 
Data from simulation of the algorithm with epsilon = 1. The QoS is max all the time and the energy is derived from the throttle (approximated). 


To generate:
* call Matlab script, 
* select 'Energy', 
* order is default, 
* characteristic time is default, 
* sensor energy data are [pprz_throttle.csv](../simulation1/pprz_throttle.csv), 
* column is 3, 
* time step 500 ms, 
* computational model and mission specification are default, 
* select OP3, 
* epsilon is 1.

Contains 4 columns:
1. time (sec). The [max,min] is = [.5,320].
2. data from sensor = [29.368,40.381].
3. estimated data [29.65,40.381].
4. zero column (just for gnuplot)

### [est_vs_joules.csv](est_vs_joules.csv)
Data which correlates when the estimation has stopped in seconds and how the energy evolves from that point on, against the actual measured energy.

Contains 3 columns:
1. time (sec). The [max,min] is = [.5,320].
2. data from sensor = [19297,24635].
3. data from sensor, constant 22334.
