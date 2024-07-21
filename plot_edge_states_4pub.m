

% vec_test = eigVecs_cell{1}(:,337);
% vec_test'*vec_test

% Identify edge-localized states
edge_criterion = 0.83;  % Define a criterion for edge localization
edge_states = false(num_eigenvalues, num_mu);

width_edge = 2;
for k = 1:num_mu
    if abs(mu_values(k)/t)>3.4
        width_edge = 3;
        edge_criterion = 0.85;
    else
        width_edge = 2;
        edge_criterion = 0.83;
    end

    for n = 1:num_eigenvalues
        pdf = pdf_p_n_h_real_space(:, :, n, k);
        % Calculate edge and center densities
        edge_density = sum(pdf([1:width_edge, end-width_edge+1:end], :), 'all') + sum(pdf(:, [1:width_edge, end-width_edge+1:end]), 'all');
        center_density = sum(pdf(width_edge+1:end-width_edge, width_edge+1:end-width_edge), 'all');
        % Compare densities to determine if state is edge-localized
        edge_states(n, k) = edge_density / (edge_density + center_density) > edge_criterion;
    end
end

% Plot 3D graph of edge states probability density for mu = 2.4
mu_target = 3.2;
[~, mu_idx] = min(abs(mu_values - mu_target));
num_edge_states = sum(edge_states(:,mu_idx));
fprintf('Total number of edge states for mu= %f: %d\n',mu_values(mu_idx), num_edge_states)

% Grid definition for plotting
[X, Y] = ndgrid(1:Nx, 1:Ny);
X = X';
Y = Y';

% f2 = figure;
% Activate figure 2
figure(2);
% Clear figure 2 only
clf;

% Set the new position and size for figure 2
% [left, bottom, width, height]
newPosition = [50, 50, 1000, 900];
set(gcf, 'Position', newPosition);

% % Retrieve the handle of the current axes
% ax = gca;
% 
% % Get the current viewpoint
% [az, el] = view(ax);
az =   64.5120;
el =   49.9054;

% Hold on to the current figure
hold on;

for n = 1:num_eigenvalues
    if edge_states(n, mu_idx)
        pdf = pdf_p_n_h_real_space(:, :, n, mu_idx);
        surf(X, Y, pdf,'FaceAlpha',0.5);% 'EdgeColor', 'none');
    end
end
xlabel('X');
ylabel('Y');
zlabel('Probability Density');
title([num2str(num_edge_states),' Edge States Probability Density at \mu = ', num2str(mu_target)]);
hold off;
view(az, el);


% Activate figure 3
figure(3);
% Clear figure 3 only
clf;

% Set the new position and size for figure 2
% [left, bottom, width, height]
newPosition = [50, 50, 1700, 900];
set(gcf, 'Position', newPosition);

% % Retrieve the handle of the current axes
% ax = gca;
% 
% % Get the current viewpoint
% [az, el] = view(ax);
az =   64.5120;
el =   49.9054;

% Create the first subplot (1 row, 2 columns, first subplot)
subplot(1, 2, 1);
hold on; % Hold on to plot multiple surfaces

for n = 1:num_eigenvalues
    if edge_states(n, mu_idx)

        eigVec = eigVecs_cell{mu_idx} ;
        psi_p = eigVec(1:Nx*Ny, n);  % Only take the electron part
        psi_p = reshape(psi_p,Nx,Ny)';
        % psi_h = eigVec(Nx*Ny+1:end, n);  % Only take the electron part
        % psi_h = reshape(psi_h,Nx,Ny)';
        surf(X, Y, abs(psi_p).^2,'FaceAlpha',0.5);% 'EdgeColor', 'none');
        break
    end
end
xlabel('X');
ylabel('Y');
zlabel('Probability Density');
title([num2str(num_edge_states),' Edge States Particle Density at \mu = ', num2str(mu_target)]);
hold off;
view(az, el);

% Create the second subplot (1 row, 2 columns, second subplot)
subplot(1, 2, 2);
hold on; % Hold on to plot multiple surfaces

for n = 1:num_eigenvalues
    if edge_states(n, mu_idx)

        eigVec = eigVecs_cell{mu_idx} ;
        % psi_p = eigVec(1:Nx*Ny, n);  % Only take the electron part
        % psi_p = reshape(psi_p,Nx,Ny)';
        psi_h = eigVec(Nx*Ny+1:end, n);  % Only take the electron part
        psi_h = reshape(psi_h,Nx,Ny)';
        surf(X, Y, abs(psi_h).^2,'FaceAlpha',0.5);% 'EdgeColor', 'none');
        break
    end
end
xlabel('X');
ylabel('Y');
zlabel('Probability Density');
title([num2str(num_edge_states),' Edge States Hole Density at \mu = ', num2str(mu_target)]);
hold off;
view(az, el);


