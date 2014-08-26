%% For running on SciClone
Lx=9; % The number of lattice sites in the x direction

N_up=5; % The number of spin-up electrons
N_dn=5; % The number of spin-down electrons

U=4.0; % The on-site repulsion strength in the Hubbard Hamiltonian
t=1; % The hopping amplitude between nearest-neighbor sites in the x direction

tau = 0.5;

maxNumCompThreads(8)

[ exactUp, exactDn ] = unequalTimeGF( t, U, tau, Lx, N_up, N_dn );
