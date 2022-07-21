
# Energy-aware planning-scheduling for autonomous aerial robots 

This is the camera ready version of our planning-scheduling [paper](pdf/energy-planning_3.pdf) that we will be presenting at IROS'22. 


## Reference

To reference the paper, it is sufficient to use the following BibTeX entry:
```bibtex
@inproceedings{seewald2022energy,
  title={Energy-aware coverage planning and scheduling for autonomous aerial robots},
  author={Seewald, Adam and Garc{\'i}a de Marina, H{\'e}ctor and Midtiby, Henrik Skov and Schultz, Ulrik Pagh},
  booktitle={Proceedings of the IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS'22)},
  pages={8},
  year={2022},
  organization={IEEE},
  url={https://adamseewald.cc/short/energy2021},
  note={To appear}
}
```

## Simulator

The simulator that generates the data in the paper can be invoked [here](scripts/MAIN2.m). It is written in MATLAB (R); we have version 9.9.0.1467703 (R2020b).

## Data

* Data from the Opterra fixed-wing UAV flights are stored in [data/simulation1](data/simulation1)
* Data from different simulations are stored in [data/simulation3](data/simulation3)

### Plan I (static)

For the static plan I with no re-planning-scheduling; an agricultural survey (path) with hazard detection (computation), wind speed 5 m/s, and wind direction 0 degrees

* all data: [data/simulation3/raw5/new_physics/static2](data/simulation3/raw5/new_physics/static2)
    * position: [position_simulation3Ds.csv](data/simulation3/raw5/new_physics/static2/position_simulation3Ds.csv)
    * energy: [energy_simulation3Ds.csv](data/simulation3/raw5/new_physics/static2/energy_simulation3Ds.csv)
    * period: [perioddata_simulation3Ds.csv](data/simulation3/raw5/new_physics/static2/perioddata_simulation3Ds.csv)

To generate the data
* simulator: [scripts/SIM11.m](scripts/SIM11.m)
* Matlab data: [3Ds.mat](data/simulation3/raw5/new_physics/static2/3Ds.mat)

### Plan II (static)

For the static plan II, same as I but wind direction 90 degrees

* all data: [data/simulation3/raw1/new_physics/static2](data/simulation3/raw1/new_physics/static2)
    * position: [position_simulation3s.csv](data/simulation3/raw1/new_physics/static2/position_simulation3s.csv)
    * energy: [energy_simulation3s.csv](data/simulation3/raw1/new_physics/static2/energy_simulation3s.csv)
    * period: [perioddata_simulation3s.csv](data/simulation3/raw1/new_physics/static2/perioddata_simulation3s.csv)

To generate the data
* simulator: [scripts/SIM10.m](scripts/SIM10.m)
* Matlab data: [3s.mat](data/simulation3/raw1/new_physics/static2/3s.mat)


### Plan i (dynamic)

For the plan i with re-planning-scheduling of path and computation parameters simultaneously (the conditions are the same as for Plan I)

* all data: [data/simulation3/raw5/new_physics/dynamic_revised](data/simulation3/raw5/new_physics/dynamic_revised)
    * position: [position_simulation3D.csv](data/simulation3/raw5/new_physics/dynamic_revised/position_simulation3D.csv)
    * energy: [energy_simulation3D.csv](data/simulation3/raw5/new_physics/dynamic_revised/energy_simulation3D.csv)
    * period: [perioddata_simulation3D.csv](data/simulation3/raw5/new_physics/dynamic_revised/perioddata_simulation3D.csv)
    * parameters: [ctl_simulation3D.csv](data/simulation3/raw5/new_physics/dynamic_revised/ctl_simulation3D.csv)

To generate the data and see the algorithm
* simulator: [scripts/SIM8_revised.m](scripts/SIM8_revised.m)
* Matlab data: [3D.mat](data/simulation3/raw5/new_physics/dynamic_revised/3D.mat)

The plan starts at the highest configuration of parameters

In the plan there is an unexpected battery drop at 93 and 270 seconds of 6 and 9 % respectively

### Plan ii (dynamic)
For the plan i with re-planning-scheduling of path and computation paremeters simultaneously (the conditions are the same as for Plan II)

* all data: [data/simulation3/raw1/new_physics/dynamic_revised](data/simulation3/raw1/new_physics/dynamic_revised)
    * position: [position_simulation3.csv](data/simulation3/raw1/new_physics/dynamic_revised/position_simulation3.csv)
    * energy: [energy_simulation3.csv](data/simulation3/raw1/new_physics/dynamic_revised/energy_simulation3.csv)
    * period: [perioddata_simulation3.csv](data/simulation3/raw1/new_physics/dynamic_revised/perioddata_simulation3.csv)
    * parameters: [ctl_simulation3.csv](data/simulation3/raw1/new_physics/dynamic_revised/ctl_simulation3.csv)

To generate the data and see the algorithm
* simulator: [scripts/SIM9_revised.m](scripts/SIM9_revised.m)
* Matlab data: [3.mat](data/simulation3/raw1/new_physics/dynamic_revised/3.mat)

The plan starts at the lowest configuration of parameters and the battery behaves linearly.

## Plots

To generate the plots, we use gnuplot (version 5.2).

* Illustrative abstract (Figure 1), Figures 2--4 are all illustrations made using vector graphic software and imported as TiKZ pictures
* Energy and spectrum plots from Figure 1: [scripts/plots/abstract/graphs.plot](scripts/plots/abstract/graphs.plot) run `gnuplot graphs.plot`
* Figures 5--6 are then a merge of multiple plots. For trajectories: [scripts/plots/results/new_physics2/trajs.plot](scripts/plots/results/new_physics2/trajs.plot)
* Whereas for energies: [scripts/plots/results/new_physics2/ener.plot](scripts/plots/results/new_physics2/ener.plot)
* Finally for energy models: [scripts/plots/results/new_physics2/model_1.plot](scripts/plots/results/new_physics2/model_1.plot) and [scripts/plots/results/new_physics2/model_2.plot](scripts/plots/results/new_physics2/model_2.plot)


