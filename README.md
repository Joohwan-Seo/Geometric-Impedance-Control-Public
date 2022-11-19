# Geometric-Impedance-Control
Tested in Matlab R2021a and R2021b

Written by Joohwan Seo, Ph.D. Student in Mechanical Engineering, UC Berkeley.

## Authors' comment
This is quite raw, unorganized files. I hope everyone can get some of the insights, and please forgive me for the nasty codes.


## Main files
`main_geo_discrete.m` runs the simulation for the proposed approach \
`main_imp_discrete.m` runs the simulation for the benchmark approach \
`plotter.m` for visualizing the result.

Switch between the `tracking` and `regulation` to select the obejctive.

## Trajectory files
`desired_trajectory.m` gives the trajectory for the simulation in body-frame, utilized for the proposed control.\
`desired_trajectory_spatial.m` gives the trajectory for the simulation in spatial frame, utilized for the benchmark approach.\
`desired_trajectory_regulation.m` give just desired orientation and point, both in spatial frame. 

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
