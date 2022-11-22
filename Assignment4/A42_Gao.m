%% Information
% Title: AESM1511-Assignment4.2
% Time: 10-11-2022
% Name: Boyu Gao + Miaoyang Yuan + Yongting Zhu
% Student Number: 5794676 + 5732476 + 5846528
% E-mail: B.Gao@student.tudelft.nl + yuanmy2318@gmail.com + Y.Zhu-44@student.tudelft.nl

%% Assignment 4.2
% Title: Scattering from a multilayered earth —— a numerical solution by finite differences
% 1) Draw a flowchart
% 2) Define the discrete parameters (electric permittivity, electric conductivity, magnetic permeability) 
%    models from layer 1 to layer 4 in the domain D2. 
%    All these 3 parameters modeles are time independent, only z dependent.
% 3) Precomupte the incident electric field for all discrete times.
% 4) Triple loop calculate the electric and magnetic field in D2 domain.
%    The triple for-loops should be written in a function —— A42_Task4
% 5) Achieve the same result in task 4 using a dual loop (drop the space-loop).
% 6) Halve the z-step and calculate the energy in the differences.
% 7) Visualization.

%% Global parameters
% 1) Constants
% c0              the free space velocity. [m/s]
% mu0             the free space magnetic permeability. [H/m]
% e0              the free space electric permittivity. [F/m]
% 2) Layers parameters (parameters in the table)
% e_r             relative values of the electric permittivity of all layers (6)
% e_a             actual values of the electric permittivity of all layers (6)
% mu_r            relative values of the magnetic permeability of all layers (6)
% mu_a            actual values of the magnetic permeability of all layers (6)
% rho             electric conductivity of all layers (6)
% d               thickness of all layers (6)  [m]
% c               velocities of all layers (6) [m/s]
% 3) Frequency parameters
% f_c             the centre frequencies of the Ricker wavelet [Hz]
% f_max           the max frequency = 4*f_c [Hz]
% 4) Grid parameetrs 
% t_w             the total time window of the computation [s]
% d2_layers       the total layers (4) of the D2 domain [m]
% 5) Figure parameters
% fontsize        font size = 16
% fontname        font name is 'Arial'
% cbar_scale      scale the colarbar axis to 25% of the max and min of the amplitude values 

%% Main and Function file lists:
% 1) This file is the main function. 
% 2) A42_Task4    Triple loop calculate the electric and magnetic field in D2 Domain
% 3) A42_Task5    Dual loop calculate the electric and magnetic field in D2 Domain 
%    Note: the results of A42_Task4 and A42_Task5 are the same!
% 4) A42_Task7    Plot function 


%% Clean sheet
clear variables  % remove items from workspace, freeing up system memory
close all        % close any open figures
clc              % and clear the command window


%% Global parameters
% 1) Constants
c0 = 299792458;         
mu0 = 4 * pi * 10^(-7); 
e0 = 1 / (mu0*c0^2);
% 2) Layers parameters (table)
e_r = [1, 6, 2, 16, 6, 9];  
e_a = e_r * e0;
mu_r = [1, 1, 1, 1, 1, 1];
mu_a = mu_r * mu0;
rho = [0, 10^(-3), 10^(-4), 10^(-2), 10^(-3), 10^(-3)]; 
d = [0.24, 0.10, 0.003, 0.05, 0.15];
c = c0 ./ sqrt(mu_r.*e_r); 
% 3) Frequency parameters
f_c = [500*10^6, 1*10^9, 2*10^9];
f_max = 4*f_c;
% 4) Grid parameters
t_w = 12 * 10^(-9);
d2_layers = 4;
% 5) Figure parameters
fontsize = 16;
fontname = 'Arial';
cbar_scale = 0.25; 


%% Task 2: Define the discrete parameters models
% 1) Define time-step and z-step
t_step_criterion = 1 / (2*max(f_max));             % the stability criterion of the max time-step
z_step_criterion = max(c) * t_step_criterion;      % the stability criterion of the max z-step: 0.0187
% 2) Initialize the actual time-step and z_step
% Note: to ensure there are at least 2 (z) grid points in each layer, and that
% these two grid points are not on the boundry, so we should multiply a
% value which is slightly less than 0.5 (we choose 0.49 here).
t_step = t_step_criterion;                         % Initialize the actual time-step 
z_step = 0.49 * min(d);                            % Initialize the actual z-stpe: 0.00147
% 3) Automatically check and adjust the t_step and z_step
disp('Task 2: Automatic check and adjust the z_step and t_steps.');
if z_step <= z_step_criterion
    t_step = z_step / max(c);
    disp('The current z_step is small enough to meet the stability criterion, and the code accordingly adjusts the t-step.');
else
    z_step = z_step_criterion;
    disp('The current z_step is too large, so the accordingly adjusts the z_step.');
end
% 4) Discretize the time and z of the D2 domain
% Calculate the depth of each layer (D2 domain: from layer 1 to layer 4) of the interface
d2_depths = cumsum(d) - d(1);               % the interfaces depth
d_d2 = d2_depths(end);                      % the thickness of d2 domain: layer1 to layer4 
z_points = linspace(0, d_d2, d_d2/z_step);  % z total points: 206             
t_points = linspace(0, t_w, t_w/t_step);    % t total points: 2447
total_z_points = length(z_points);
total_t_points = length(t_points);
% 5) Define the discrete parameters (e, rho, mu) models
% Note: all parameters (e, rho, mu) models here are time independent, only z dependent!
e_z = zeros(total_z_points, 1);             % Initialize the electric permittivity model
rho_z = zeros(total_z_points, 1);           % Initialize the electric conductivity model
mu_z = zeros(total_z_points, 1);            % Initialize the magnetic permeability model
% Loop: the outer is about layers, the inner is about models
for l_tmp = 2:length(d2_depths)             % 2 to 5 —— also corresponds to e_a, rho, mu_a
    for z_tmp = 1:total_z_points            % 1 to 206
       if z_points(z_tmp) <= d2_depths(l_tmp) && z_points(z_tmp) >= d2_depths(l_tmp-1)
            e_z(z_tmp) = e_a(l_tmp);
            rho_z(z_tmp) = rho(l_tmp);
            mu_z(z_tmp) = mu_a(l_tmp);
       end
    end
end
fprintf('The final z_step = %f (m), the final t_step = %e (s).\n', z_step, t_step);
fprintf('Total z points is %d, total t points is %d.\n', total_z_points, total_t_points);
disp('---');


%% Task 3: Precompute the incident electric field for all discrete times 
% 1) Define local parameters
Z0 = sqrt( mu0 / e0 );       % constant: 1
tau = sqrt(2) ./ f_c;        % constant: 3
tt = t_points - d(1)/c(1);   % the actual t used in the wavelet W function: 2447
% 2) Calculate the equation (38) using equation (39)
E_incident = zeros(1, total_t_points, length(f_c));   % Initialize the incident electric field: 1 x 2447 x 3
% Note1: we first calculated the derivative of W with respect to time, so we
% directly use the derivative in equation (38)
% E_incident will be used in Task 4, equation (37): 1 x 2447 x 3
for f_tmp = 1 : length(f_c)
    E_incident(1, :, f_tmp) = Z0 * (2*pi^2) .* (1./tau(f_tmp)) * ((tt./tau(f_tmp))-1) .* ...
                                exp(-2*pi^2*(tt./tau(f_tmp)-1).^2) .* (3-4*pi^2*(tt./tau(f_tmp)-1).^2);        
end


%% Task 4: Compute 3 time-domain solutions for the electric + magnetic field everywhere in D2
% 1) Initialize the electric + magnetic fields
Ey_field = zeros( total_z_points, total_t_points+2, length(f_c) );   % Initialize the electric field in D2: 206 x 2449 x 3
Hx_field = zeros( total_z_points, total_t_points+2, length(f_c) );   % Initialize the magnetic field in D2: 206 x 2449 x 3
% 2) Call the function A42_Task4 and get two return values
% Note: finally in task 7 we will plot with Ey_field and Hx_field
[Ey_field, Hx_field] = A42_Task4(Ey_field, Hx_field, ...
                                 mu_a, e_a, rho, c, f_c, z_step, t_step, ...
                                 e_z, rho_z, mu_z, ...
                                 total_t_points, total_z_points, ...
                                 E_incident);


%% Task 5: Write a function without an explicit loop over points in space
% Note: in task 4 we use a triple for-loop, and the innermost loop is about the z space points.
% So in task 5, we will just use a double for-loop without the sapce-loop to achieve the same result.
% 1) Initialize the electric + magnetic fields: the same size as in task 4
task5_Ey_field = zeros( total_z_points, total_t_points+2, length(f_c) );   % Initialize the electric field in D2: 206 x 2449 x 3
task5_Hx_field = zeros( total_z_points, total_t_points+2, length(f_c) );   % Initialize the magnetic field in D2: 206 x 2449 x 3
% 2) Call the function A42_Task5 and get two return values
[task5_Ey_field, task5_Hx_field] = A42_Task5(task5_Ey_field, task5_Hx_field, ...
                                             mu_a, e_a, rho, c, f_c, z_step, t_step, ...
                                             e_z, rho_z, mu_z, ...
                                             total_t_points, total_z_points, ...
                                             E_incident);


%% Task 6: Halve the z-step and calculate the energy in the differences
disp("Task 6: We don't have enough time to finish this question.");
disp('---');


%% Task 7: Visualize the electric and magnetic field
% The figure of task 4
ttitle_T4 = 'T4.2-4. The electric and magnetic field in D2 domain for f_c = 500 MHZ, 1 GHZ, 2 GHZ.';
A42_Task7(cbar_scale, f_c, t_points, z_points, ...
          fontsize, fontname, ttitle_T4, ...
          Ey_field, Hx_field);
% The figure of task 5
ttitle_T7 = 'T4.2-7. The electric and magnetic field in D2 domain for f_c = 500 MHZ, 1 GHZ, 2 GHZ.';
A42_Task7(cbar_scale, f_c, t_points, z_points, ...
          fontsize, fontname, ttitle_T7, ...
          task5_Ey_field, task5_Hx_field);

% Display
disp('Task 7: The higher the frequency, the higher the resolution.');
disp('Although task 4 and task 5 use different loop calculation methods, the final results are the same.');
disp('---');


