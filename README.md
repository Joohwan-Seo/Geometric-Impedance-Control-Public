# Geometric-Impedance-Control
Tested in Matlab R2021a and R2021b

Written by Joohwan Seo, Ph.D. Student in Mechanical Engineering, UC Berkeley.

## Authors' comment
These are quite raw, unorganized files. I hope everyone can get some of the insights, and please forgive me for the nasty codes.

## Fixes & Update
### 2023/03/25
1. There was a serious error in the benchmark controller; now that it is fixed, it works far better than before.
2. For the readers to test controllers in various scenarios, more scenarios are added.
3. In `desired_trajectory2.m` and `trajectory_calculator.m`, we show a way to design a polynomial-based smooth trajectory generation technique also for the rotational part.

### 2023/04/07
0. Final version of the paper submitted.
1. Control equations are added on every main_geo_discrete file.
2. Geometric Impedance Control version 2 added
3. Additional updates on the trajectory tracking simulations.
4. Code variables now match with the variables in the paper.

### 2023/04/11
The arxiv paper is updated to the final version.

## Main files
`main_geo_discrete.m` runs the simulation for the proposed approach (intuitive geometric impedance) - Equation (30)\
`main_geo_discrete_v2.m` runs the simulation for the proposed approach (geometric impedance control-version 1, GIC-v1) - Equation (32)\
`main_geo_discrete_v3.m` runs the simulation for the proposed approach (geometric impedance control-version 2, GIC-v2) - Equation (47)\
`main_imp_discrete.m` runs the simulation for the benchmark approach (conventional impedance control, CIC)\
`main_imp_discrete_v2.m` runs the simulation for the benchmark approach, version 2 (Not rigorously tested)\
`plotter.m` for visualizing the result -- GIC-v1 vs CIC\
`plotter2.m` for visualizing the result -- GIC-v2 vs CIC\
### Note
Comparison results on `plotter2.m` are not presented in the paper since the two control laws are different, thus unable to do a fair comparison.

`plotter_geo_comp.m` is for comparison between intuitive geometric impedance control and GIC-v1. When the gains become different, they start to behave differently.\
When gains are scalar in translational and rotational ways, both controllers are identical.

## Objectives in Main files
Change `obj` variable to `tracking`, `tracking2`, `regulation`, and `regulation2`.\
`tracking`: sinusoidal trajectory tracking (presented in paper)\
`tracking2`: 3rd-order smooth polynomial trajectory tracking both in translational and rotational\
`regulation`: regulation task\
`regulation2`: Multi-point regulation, can be considered as step-input case.

`desired_trajectory.m` gives the trajectory utilized in `tracking` objective.\
`desired_trajectory2.m` gives the trajectory utilized in `tracking2` objective.\
`desired_trajectory_regulation.m` gives the desired setpoint utilized in `regulation` objective.\
`desired_trajectory_regulation2.m` gives the desired setpoints utilized in `regulation2` objective.

## "function_generator" folder
`RTB matlab` (Robotics Toobox Matlab) is needed to run the code. \
Build the UR5e robot model and associated dynamic parameter matrices as well as Jacobian matrices.

## "sub_direct" folder
Dynamic parameter matrices built from function_generator are saved here. Some miscellaneous functions are also defined here.

## "results" folder
The simulation result data generated by the main files are saved here.

## Presented at
IFAC World Congress 2023, Yokohama, Japan:
"Geometric Impedance Control on SE(3) for Robotic Manipulators"

Arxiv submitted version:
https://doi.org/10.48550/arXiv.2211.07945

Official Proceedings version:
https://doi.org/10.1016/j.ifacol.2023.10.1581

## Bibtex Citation
@article{seo2023geometric,
  title={Geometric impedance control on SE (3) for robotic manipulators},\
  author={Seo, Joohwan and Prakash, Nikhil Potu Surya and Rose, Alexander and Choi, Jongeun and Horowitz, Roberto},\
  journal={IFAC-PapersOnLine},\
  volume={56},\
  number={2},\
  pages={276--283},\
  year={2023},\
  publisher={Elsevier}\
}
