% Plot function
function [] = A42_Task7(cbar_scale, f_c, t_points, z_points, ...
                        fontsize, fontname, ttitle, ...
                        Ey_field, Hx_field)
% Note: every time this plot function is called, a new figure will be created.
% This plot function doesn't have return values. 
figure();         % create a new figure 
sgtitle(ttitle);  % the total title
for f_tmp = 1:length(f_c)
    % The electric field: 1-3
    subplot(2, 3, f_tmp)
    imagesc(t_points, z_points, Ey_field(:, :, f_tmp));
    xlabel('Time [ns]', 'FontSize', fontsize);
    ylabel('Vertical distance [m]', 'FontSize', fontsize);
    title( ['E_y Field [V/m]: f_c = ' num2str(f_c(f_tmp)/1000000) ' MHz'], 'FontSize', fontsize, 'FontName', fontname );
    colormap('jet');        
    colorbar;
    caxis( [cbar_scale * min(min(Ey_field(:, :, f_tmp))), cbar_scale * max(max(Ey_field(:, :, f_tmp)))] ); 
    % The magnetic field: 4-6
    subplot(2, 3, f_tmp+3)
    imagesc(t_points, z_points, Hx_field(:, :, f_tmp));
    xlabel('Time [ns]', 'FontSize', fontsize);
    ylabel('Vertical distance [m]', 'FontSize', fontsize);
    title( ['H_x Field [A/m]: f_c = ' num2str(f_c(f_tmp)/1000000) ' MHz'], 'FontSize', fontsize, 'FontName', fontname );
    colormap('jet');        
    colorbar;
    caxis( [cbar_scale * min(min(Hx_field(:, :, f_tmp))), cbar_scale * max(max(Hx_field(:, :, f_tmp)))] ); 
end