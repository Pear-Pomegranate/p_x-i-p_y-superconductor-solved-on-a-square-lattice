% main.m

clc;

% Step 1: Initialize Parameters
Nx = 20;                % Lattice size in x-direction
Ny = 20;                % Lattice size in y-direction
t = 1;                  % Hopping parameter
Delta_0 = 2;            % p-wave pairing potential
mu_start = 0;          % Start value of chemical potential
mu_end = 5;             % End value of chemical potential
mu_step = 0.2;          % Step size for chemical potential

mu_values = mu_start:mu_step:mu_end;  % Array of mu values
num_mu = length(mu_values);           % Number of mu values
num_eigenvalues = 2*Nx * Ny;            % Number of eigenvalues to compute

% Initialize array to store eigenvalues
eigE_matrix = zeros(num_eigenvalues, num_mu);
% Initialize cell array to store eigenvectors
eigVecs_cell = cell(1, num_mu);

%pdf_p_n_h_real_space
% Initialize array to store probability densities
pdf_p_n_h_real_space = zeros(Nx, Ny, num_eigenvalues, num_mu);
Dipole_x = zeros(1,1,num_eigenvalues,num_mu);

% [X,Y] = ndgrid(1:Nx,1:Ny);
% X = X';
% Y = Y';

% Loop over mu values
for k = 1:num_mu
    mu = mu_values(k);
    
    % Update parameters
    parameters = [Nx, Ny, t, mu, Delta_0];
    
    % fprintf('Generating H_BdG:\n')
    % Generate the BdG Hamiltonian
    H_BdG = Generate_H_BdG(parameters);
    
    % tic;
    % fprintf('Finding Eigen States for H_BdG:')
    % Find eigenvalues and eigenvectors
    [eigVec, eigE] = eigs(H_BdG, num_eigenvalues, 'smallestabs');
    % elapsed_time = toc;
    % disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);
    
    % Convert the diagonal matrix eigE to a vector
    eigE = diag(eigE);
    
    % Sort eigenvalues and eigenvectors
    [eigE, sortIdx] = sort(eigE);
    eigVec = eigVec(:, sortIdx);
    
    % Store eigenvalues and eigenvectors
    eigE_matrix(:, k) = eigE;
    eigVecs_cell{k} = eigVec;


    for n = 1:num_eigenvalues
        psi_p = eigVec(1:Nx*Ny, n);  % Only take the electron part
        psi_p = reshape(psi_p,Nx,Ny)';
        psi_h = eigVec(Nx*Ny+1:end, n);  % Only take the electron part
        psi_h = reshape(psi_h,Nx,Ny)';
        pdf_p_n_h_real_space(:, :, n, k) = abs(psi_p).^2+abs(psi_h).^2;
    end
    
    % Display progress every 10%
    if mod(k, round(num_mu / 10)) == 0
        fprintf('Progress: %.0f%% k=%d out of %d\n', k / num_mu * 100,k, length(mu_values));
    end
    
end


% Plot eigenvalues as a function of mu
figure;
hold on;
for i = 1:num_eigenvalues
    plot(mu_values, eigE_matrix(i, :));
end
xlabel('Chemical Potential \mu');
ylabel('Eigen Energies');
title('Eigen Energies as a Function of Chemical Potential \mu');
hold off;

% Identify edge-localized states
edge_criterion = 0.8;  % Define a criterion for edge localization
edge_states = false(num_eigenvalues, num_mu);

width_edge = 3;
for k = 1:num_mu
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
mu_target = 2.4;
[~, mu_idx] = min(abs(mu_values - mu_target));
fprintf('Total number of edge states for mu= %f: %d\n',mu_values(mu_idx), sum(edge_states(:,mu_idx)))

% Grid definition for plotting
[X, Y] = ndgrid(1:Nx, 1:Ny);
X = X';
Y = Y';

figure;
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
title(['Edge States Probability Density at \mu = ', num2str(mu_target)]);
hold off;

% Optional: Save the eigenvalues and eigenvectors for later use
save('eigE_matrix.mat', 'eigE_matrix');
save('eigVecs_cell.mat', 'eigVecs_cell');
