# Geometric-Impedance-Control
Tested in Matlab R2021a and R2021b

Written by Joohwan Seo, Ph.D. Student in Mechanical Engineering, UC Berkeley.

## Authors' comment
This is quite raw, unorganized files. I hope everyone can get some of the insights, and please forgive me for the nasty codes.


## Main files
`main_geo_discrete.m` runs the simulation for the proposed approach (intuitive geometric impedance)\
`main_imp_discrete_v2.m` runs the simulation for the proposed approach (geometric impedance)\
`main_imp_discrete.m` runs the simulation for the benchmark approach \
`main_imp_discrete_v2.m` runs the simulation for the benchmark approach, version 2 \
`plotter.m` for visualizing the result.

`plotter_geo_comp.m` for comparison between two geometric controllers. When the gains become different, they started to behave different.\
When gains are scalar in translational and rotational way, the both controllers are identical.

## Objectives in Main files
Change `obj` variable to `tracking`, `tracking2`, `regulation`, and `regulation2`.
`tracking`: sinusoidal trajectory tracking (presented in paper)\
`tracking2`: 3rd-order smooth polynomial trajectory tracking both in translational and rotational\
`regulation`: regulation task\
`regulation2`: Multi-point regulation, can be considered as step-input case.\

`desired_trajectory.m` gives the trajectory utilized in `tracking` objective.\
`desired_trajectory2.m` gives the trajectory utilized in `tracking2` objective.\
`desired_trajectory_regulation.m` gives the desired setpoint utilized in `regulation` objective.\
`desired_trajectory_regulation2.m` gives the desired setpoints utilized in `regulation2` objective.\

## "function_generator" folder
`RTB matlab` (Robotics Toobox Matlab) is needed to run the code. \
Build the UR5e robot model and associated dynamic parameter matrices as well as Jacobian matrices.

## "sub_direct" folder
Dynamic parameter matrices built from function_generator is saved here. Some miscelleneous functions are also defined here.

## Submitted to
IFAC world congress 2023:
"Geometric Impedance Control on SE(3) for Robotic Manipulators"

arxiv submitted version:
https://doi.org/10.48550/arXiv.2211.07945

## Bibtex Citation
@article{seo2022geometric, \
  title={Geometric Impedance Control on SE (3) for Robotic Manipulators},\
  author={Seo, Joohwan and Prakash, Nikhil Potu Surya and Rose, Alexander and Horowitz, Roberto},\
  journal={arXiv preprint arXiv:2211.07945},\
  year={2022}\
}
