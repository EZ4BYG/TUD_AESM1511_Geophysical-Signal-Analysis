% Dual loop calculate the electric and magnetic field in D2 Domain
function [task5_Ey_field, task5_Hx_field] = A42_Task5(task5_Ey_field, task5_Hx_field, ...
                                                      mu_a, e_a, rho, c, f_c, z_step, t_step, ...
                                                      e_z, rho_z, mu_z, ...
                                                      total_t_points, total_z_points, ...
                                                      E_incident)
% 1) Define coefficients used in equation (26) + (27) and (34) + (37) —— totally Same as the function A42_Task4
% Coefficients in equation (37): z = 0 / upper boundry condition —— D1 domain
z0_bd_coef1 = mu_a(1) / mu_a(2);
z0_bd_coef2 = z_step / (c(1)*t_step);
z0_bd_coef3 = 2/3 * z_step / c(1);
z0_bd_coef4 = 1/3 * z0_bd_coef1;
z0_bd_coef5 = 1/3 * z0_bd_coef2;   
% Coefficients in equation (26) and (27): D2 domain
d2_bd_coef1 = t_step / z_step;
d2_bd_coef2 = mu_a(6) / mu_a(5);   
% Coefficients in equation (34): z = d / lower boundry conditon —— D3 domain
zd_bd_coef1 = z_step / (c(6)*t_step);
zd_bd_coef2 = rho(6) * sqrt(mu_a(6)/e_a(6)) * z_step / 3;  % Z3 = sqrt(mu3/e3)
zd_bd_coef3 = 1/3 * (mu_a(6) / mu_a(5));
zd_bd_coef4 = 1/3 * zd_bd_coef1;   
% 2) Dual Loop calculation: frequency > time. We drop the space loop!
% First loop: centre frequency
disp('Task 5: Dual Loop calculation of the electric and magnetic field in D2 domain.');
for f_tmp = 1:length(f_c)
    fprintf('Loop calculation begins: f_c = %f Hz\n', f_c(f_tmp));
    % Second loop: time
    for t_tmp = 3:total_t_points
        % z = 0 / upper boundry condition: equation (37) —— Same
        task5_Ey_field(1, t_tmp, f_tmp) = (z0_bd_coef4 * (4*task5_Ey_field(1+1,t_tmp,f_tmp) - task5_Ey_field(1+2,t_tmp,f_tmp)) +  ...
                                           z0_bd_coef5 * (4*task5_Ey_field(1,t_tmp-1,f_tmp) - task5_Ey_field(1,t_tmp-2,f_tmp)) +  ...
                                           z0_bd_coef3 * E_incident(1,t_tmp,f_tmp)) / (z0_bd_coef1 + z0_bd_coef2);    
       
        % Cuz we drop the space loop, we now need to calculate Hx and Ey at the same time
        % Note: many of the original coefficients become arrays !!!
        % when z = 1: equation (27)
        task5_Hx_field(1, t_tmp, f_tmp) = task5_Hx_field(1, t_tmp-1, f_tmp) + ...
                                          (d2_bd_coef1/mu_z(1)) * (task5_Ey_field(1+1, t_tmp, f_tmp) - task5_Ey_field(1, t_tmp, f_tmp));
        % when z > 1, equation (27)
        task5_Hx_field(2:end-1, t_tmp, f_tmp) = task5_Hx_field(2:end-1, t_tmp-1, f_tmp) + ...
                                                  (d2_bd_coef1./mu_z(2:end-1)) .* (task5_Ey_field(3:end, t_tmp, f_tmp) - task5_Ey_field(2:end-1, t_tmp, f_tmp));
        % when z > 1, equation (26)
        d2_bd_coef3 = e_z(2:end-1) - rho_z(2:end-1)*t_step/2;
        d2_bd_coef4 = e_z(2:end-1) + rho_z(2:end-1)*t_step/2;
        d2_bd_coef5 = 1 ./ d2_bd_coef4 * d2_bd_coef1;
        d2_bd_coef6 = d2_bd_coef3 ./ d2_bd_coef4;
        task5_Ey_field(2:end-1, t_tmp+1, f_tmp) = d2_bd_coef6 .* task5_Ey_field(2:end-1, t_tmp, f_tmp) + ...
                                                  d2_bd_coef5 .* ( task5_Hx_field(2:end-1, t_tmp, f_tmp) - task5_Hx_field(1:end-2, t_tmp, f_tmp) );

        % Back to the second loop: time
        % z = d / lower boundry condition: equation (34) —— Same 
        task5_Ey_field(total_z_points, t_tmp+1, f_tmp) = (zd_bd_coef3 * (4*task5_Ey_field(total_z_points-1, t_tmp+1, f_tmp) - task5_Ey_field(total_z_points-2, t_tmp+1, f_tmp)) + ...
                                                    zd_bd_coef4 * (4*task5_Ey_field(total_z_points, t_tmp, f_tmp) - task5_Ey_field(total_z_points, t_tmp-1, f_tmp))) / ...
                                                    (d2_bd_coef2 + zd_bd_coef1 - zd_bd_coef2);
    end
    fprintf('Loop calculation ends.\n');
end
% 3) Select the points in D2 domain —— Same
task5_Ey_field = task5_Ey_field(:, 3:end, :);
task5_Hx_field = task5_Hx_field(:, 3:end, :);
disp('---');
