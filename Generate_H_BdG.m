function H_BdG = Generate_H_BdG(parameters)

% Step 1: Initialize Parameters
Nx = parameters(1);                % Lattice size in x-direction
Ny = parameters(2);                % Lattice size in y-direction
t = parameters(3);                 % Hopping parameter
mu = parameters(4);                % Chemical potential
Delta_0 = parameters(5);           % p-wave pairing potential


% Start timing
% tic;

% Total number of sites
num_sites = Nx * Ny;

i_x = (1:Ny-1)'+ Ny*(0:Nx-1);
i_x = i_x(:)';
idx_ix = [ i_x, i_x+1];
idx_jx = [ i_x+1, i_x];

ones_i_x = ones(1,length(i_x));
val_tx = repmat(-t, 1, 2*length(i_x));
val_px = [ones_i_x*(-Delta_0/2), ones_i_x*(Delta_0/2)];

i_y = 1:(Nx-1)*Ny;
idx_iy = [i_y, i_y+Ny];
idx_jy = [i_y+Ny, i_y];

ones_i_y = ones(1,length(i_y));
val_py = [(1i*Delta_0/2)*ones_i_y, (-1i*Delta_0/2)*ones_i_y];
val_ty = [(-t)*ones_i_y, (-t)*ones_i_y];

idx_ii = 1:Nx*Ny;
val_mu = -mu*ones(1, length(idx_ii));

% sparse matrix for H_0
r_idx_H0 = [idx_ix, idx_iy, idx_ii];
c_idx_H0 = [idx_jx,idx_jy, idx_ii];
val_H0   = [val_tx, val_ty, val_mu];
H0 = sparse(r_idx_H0,c_idx_H0,val_H0,num_sites,num_sites);

%sparse matrix for px+ipy pairing
r_idx_Delta_p = [idx_ix, idx_iy];
c_idx_Delta_p = [idx_jx, idx_jy];
val_Delta_p   = [val_px, val_py];
Delta_p = sparse(r_idx_Delta_p,c_idx_Delta_p,val_Delta_p,num_sites,num_sites);

% Construct the full BdG Hamiltonian
H_BdG = [H0, Delta_p; Delta_p', -H0];

% End timing and display the elapsed time
% elapsed_time = toc;
% disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);

end

