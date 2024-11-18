# GCS
Geometric Configuration Similarity

This repository provides access to the Matlab code to reproduce the results presented in the paper:
"A Geometric Approach for the Comparison of Kinematic Synergy Postures" by Tessari, West, and Hogan (2024).

Requirements:
- Matlab2023b and next generations.
- Simulink, Simscape and the Multibody Toolboxes.

The repository contains the following files:
- "comparing_synergies.m" -> This is the main script. Run this to reproduce the papers' results.
- "hand_model_calibration.m" -> This is the calibration script. Run this everytime you make a change to the 3D model of the hand.
- "T_WH_piano_data.mat" -> This file contains the experimental data used to showcase the validity of the proposed approach and can be downloaded here: https://www.dropbox.com/scl/fi/z0ci1ac676r9lx74n1qi2/T_WH_piano_data.mat?rlkey=7brap1x8cf4dv9qzeyrgi7g6r&st=7bf7xenu&dl=0 *.
- "functions/GCS.m" -> This function computes the Geometric Configuration Similarity according to Equation 6 in the paper.
- "functions/geomSim.m" -> This function uses GCS to generate the figures of the paper and simulate the hand 3D model.
- "functions/joint_traj_gen.m" -> This function generates the joint trajectories for the hand 3D model
- "functions/getCosineSimilarity.m" -> This function computed the the Cosine Similarity.
- "functions/mean_d_max_15_final.mat" -> This file contains the calibrationd data required to correctly compute the GCS.
- "functions/HandModel/hand_simulator_v4.slx" -> This is the Simulink/Simscape model of the hand.
- "functions/HandModel/HumanHand/*CADfiles*.STEP" -> This folder contains all the CAD files of the the different hand parts required for the simulation.

*The data are protected by a password and available upon reasonable request.

For any questions or issues feel free to reach out to ftessari@mit.edu
