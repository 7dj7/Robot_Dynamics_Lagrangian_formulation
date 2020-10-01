# Robot_Dynamics_Lagrangian_formulation
## Calculates the dynamic equations of motion of a N DOF robot arm using the Lagrangian formulation

### Description of files
'compute_dh_matrix.m' - function to calculate the homogeneous transformation matrix according to the DH Parameter <br>
'inertia_tensor.m' - function to initialize the inertia tensor with respect to the body fixed frame <br>
'complete_dynamics.m' - main script file to calculate the dynamic equations of motion (torqe/force equations)

### Notes
As an example, the dynamics of a 7 DOF robot arm called the Swayer Arm has been calculated <br>
![sawyer_figure](https://user-images.githubusercontent.com/35372234/94787969-99a9ab80-03f0-11eb-802c-90024dafa024.png)

One needs to change the DOF value and the DH matrix in order obtain the desired equations of motion
