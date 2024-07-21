
# px+ipy Superconductor Solved on a Square Lattice

Tight Binding Model for 2D px+ipy superconductor solved on a square lattice.

The Code includes a highly efficient function to generate the tight binding Hamiltonian. Diagnoalization, identification of edge modes, and plotting the probability density for each states under the Nambu spinor basis. (The first half of the components of the eigen vectors are for electron on corresponding lattice, and the second half is for hole-like states on the corresponding sites.)

The model Hamiltonian is of the form 

$$
\hat{H} = 
\begin{pmatrix}
\hat{H}_0 & \hat{\Delta} \\
\hat{\Delta}^{\dagger} & -\hat{H}_0
\end{pmatrix}
$$

Here 


$$H_0 = -\mu - t\sum_{\vec{n}}(c_{\vec{n}+\vec{\delta}}^{\dagger}c_{\vec{n}} + h.c.)$$

where $\vec{\delta}$  denotes the nearest neighbor displacement.

$$
\hat{\Delta} = \Delta_0  \sum_{\vec{n}} (c_{\vec{n}} c_{\vec{n}+\vec{\delta}x} + c_{\vec{n}+\vec{\delta}x}^{\dagger} c_{\vec{n}}^{\dagger}) + i \Delta_0 \sum_{\vec{n}} (c_{\vec{n}} c_{\vec{n}+\vec{\delta}y} - c_{\vec{n}+\vec{\delta}y}^{\dagger} c_{\vec{n}}^{\dagger}) 
$$

where the indices

$$
\vec{n} = (n_x, n_y) \quad \text{denote a position on the lattice.}
$$

For a lattice size of $N_x$ by $N_y$, the $H_{BdG}$ is $2N_xN_y$ by $2N_xN_y$.

The model Hamiltonian is related to the BdG Hamiltonian via  $H = \frac{1}{2} C^{\dagger} H_{BdG} C$, where the $C$ is the Nambu spinor $C = (c_1, \ldots, c_N, c_1^{\dagger}, \ldots, c_N^{\dagger})^T$, and $N = N_x N_y$

Therefore the first half of the eigenvector for $H_{BdG}$ stands for particle (holing) state, the second half corresponds to the hole-like state.

Let $\Phi = (\vec{u}, \vec{v})^T$ where $\vec{u} = (u_1, \ldots, u_N)$, $\vec{v} = (v_1, \ldots, v_N)$, be the eigenvector of $H_{BdG}$.

The probability distribution for each sites is $p_i =  (|u_i|^2 + |v_i|^2)$

To identify an edge state one requirement is that the probability distribution on the edges are significantly greater than the rest, i.e.

$$
p_{\text{edge}} = \sum_{i \in \text{edge}} p_i > P_0,
$$

where $P_0 \in (0,1)$ is some threshold probability to distinguish the edged states from the rest.







## Authors

- [@Pear-Pomegranate](https://github.com/Pear-Pomegranate)


## Demo

- The energy bands as a function of the chemical potential
![E_vs_mu](https://github.com/user-attachments/assets/5ff50005-9a1b-4f64-a794-9178e6ab47f9)
- The probability density distribution of the edge states on the lattices
![edge_states_density_plot](https://github.com/user-attachments/assets/410824fe-8ad9-4571-aa24-54588071da11)
- The probability density distribution of the edge states on the lattices for electron and hole states respectively.
![edge_states_p_h_density_plot](https://github.com/user-attachments/assets/d7a4b3e4-76a4-4083-8124-d4841d70116f)



## License

[Educational Community License v2.0](https://https://choosealicense.com/licenses/ecl-2.0/)

