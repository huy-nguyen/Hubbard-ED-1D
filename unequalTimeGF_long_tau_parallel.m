function output_files = unequalTimeGF_long_tau_parallel( t, U, tau_start, tau_end, tau_step, noOfSites, noOfUp, noOfDn, NUM_OF_EIGEN_VALUES, sector, method, commit_number, need_profiling )
% calculate the unequal time GF by using series expansion

OPTS.issym = 1;
OPTS.isreal = 1;

expanded_space_size_up = nchoosek(noOfSites,noOfUp+1)*nchoosek(noOfSites,noOfDn);
expanded_space_size_dn = nchoosek(noOfSites,noOfUp)*nchoosek(noOfSites,noOfDn + 1);
format compact;

profile_directory_name = strcat('profile_',num2str(noOfSites, '%02d'),...
                                    '_sites_',num2str(noOfUp, '%02d'),...
                                    'u',num2str(noOfDn, '%02d'),...
                                    'd_U_',num2str(U, '%4.2f'),...
                                    '_t_',num2str(t),...
                                    '_eigen_', num2str(NUM_OF_EIGEN_VALUES, '%04d'),...
                                    ' ',datestr(now,'_yymmdd_HHMMSS'));
if strcmp( need_profiling, 'Yes' )
    fprintf('Profile directory: %s.\n', profile_directory_name )
end

                                
output_files = {};
list_of_taus = tau_start:tau_step:tau_end;
for i_filename = 1:length(list_of_taus)
    tau = list_of_taus(i_filename);
    output_files{i_filename} = strcat('ED_',num2str(noOfSites, '%02d'),...
                                        '_sites_',num2str(noOfUp, '%02d'),...
                                        'u',num2str(noOfDn, '%02d'),...
                                        'd_U_',num2str(U, '%4.2f'),...
                                        '_tau_',num2str(tau, '%4.2f'),...
                                        '_t_',num2str(t),...
                                        '_eigen_', num2str(NUM_OF_EIGEN_VALUES, '%04d'),...
                                        ' ',datestr(now,'_yymmdd_HHMMSS'),'.mat');
end
    

fprintf('Begin calculations at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
fprintf('%03d data files:\n', length(list_of_taus))
for i_dat_files = 1:length(list_of_taus)
    disp(output_files{i_dat_files});
end

fprintf('\n')

aux_file_name = strcat('aux_',num2str(noOfSites, '%02d'),...
                                    '_sites_',num2str(noOfUp, '%02d'),...
                                    'u',num2str(noOfDn, '%02d'),...
                                    'd_U_',num2str(U, '%4.2f'),...
                                    '_tau_',num2str(tau, '%4.2f'),...
                                    '_t_',num2str(t),...
                                    '_eigen_', num2str(NUM_OF_EIGEN_VALUES, '%04d'),...
                                    ' ',datestr(now,'_yymmdd_HHMMSS'),'.mat');
                                

tic;
if strcmp( need_profiling, 'Yes' )
    profile -memory on;
end

if (noOfUp < noOfSites) && (noOfDn < noOfSites)
    fprintf('Begin diagonalizing firstHamiltonian at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
    [groundState,groundStateEnergy]=eigs( hubbardHamiltonian_parallel( t, U, noOfSites, noOfUp, noOfDn ),...
                                               1,'sa'); %ASSUMING THAT THE HAMILTONIAN IS REAL SYMMETRIC
                                           
    save( aux_file_name, 'groundState', 'groundStateEnergy', '-mat', '-v7.3');
    aux_file = matfile(aux_file_name);
    
    for i_f = 1:length(output_files)
        save(output_files{i_f},'groundState','groundStateEnergy', '-v7.3');            
    end     
    clearvars i_f;
    clearvars groundState;
    fprintf('Done with diagonalization at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
    
%% SPIN UP:
    if strcmp( sector, 'up' ) || strcmp( sector, 'both' )
        fprintf('Begin diagonalizing spin-up of secondHamiltonian at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        if NUM_OF_EIGEN_VALUES >= expanded_space_size_up
            NUM_OF_EIGEN_VALUES_UP = expanded_space_size_up - 1;
            fprintf('NUM_EIGEN_VALUES exceeds dimension of spin-up matrix. Now set to %d\n', NUM_OF_EIGEN_VALUES_UP)
        else
            NUM_OF_EIGEN_VALUES_UP = NUM_OF_EIGEN_VALUES;
        end

        [eigenVectors_up, eigenValues_up] = eigs( hubbardHamiltonian_parallel( t, U, noOfSites, noOfUp+1, noOfDn ), ...
                                                NUM_OF_EIGEN_VALUES_UP, 'sa', OPTS);
        eigenValues_up = diag(eigenValues_up);
        fprintf('Done with diagonalization at time %s.\n', datestr(now,'yymmdd_HHMMSS'))

        save( aux_file_name, '-append', 'eigenVectors_up', 'eigenValues_up', '-mat', '-v7.3');
        for i_f = 1:length(output_files)
            save(output_files{i_f},'-append','eigenValues_up','eigenVectors_up', '-v7.3');            
        end     
        clearvars i_f eigenVectors_up;

        fprintf('Begin spin up calculations at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        for t_tau = 1:length(output_files) % loop over taus

            tau = list_of_taus(t_tau);
            fprintf('Working on tau = %4.2f     at time %s\n', tau, datestr(now,'yymmdd_HHMMSS'))
            spinUpGreenFunction=zeros(noOfSites);

            fprintf('        Begin calculating off-diagonal elements at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
            for i_site=1:noOfSites        
                fprintf('        Working on i = %3d     at time %s\n', i_site, datestr(now,'yymmdd_HHMMSS'))
                left_wave_function = (aux_file.groundState') * ...
                                                            creationOperator_parallel( noOfSites, noOfUp, noOfDn , i_site, 'up' )';
                save( aux_file_name, '-append', 'left_wave_function', '-mat', '-v7.3');
                clearvars left_wave_function;
                for j_site=(i_site+1):noOfSites 
                    fprintf('        Working on j =     %3d at time %s\n', j_site, datestr(now,'yymmdd_HHMMSS'))
                    right_wave_function = creationOperator_parallel( noOfSites, noOfUp, noOfDn , j_site, 'up' ) * ...
                                                        aux_file.groundState;
                    save( aux_file_name,'-append',  'right_wave_function', '-mat', '-v7.3');
                    clearvars right_wave_function;

                    k_sum = 0;
                    for k_eigenValues = 1:NUM_OF_EIGEN_VALUES_UP %sum over k
                        expo_factor = exp( tau*( groundStateEnergy - eigenValues_up(k_eigenValues)));                
                        i_total = aux_file.left_wave_function * aux_file.eigenVectors_up(:,k_eigenValues);
                        j_total = dot( aux_file.right_wave_function, conj(aux_file.eigenVectors_up(:,k_eigenValues)) );
                        k_sum = k_sum + expo_factor * i_total * j_total;

                        clearvars expo_factor i_total j_total;
                    end
                    spinUpGreenFunction(i_site,j_site) = k_sum;
                    clearvars creationMatrix right_wave_function k_sum;
                    right_wave_function = [];
                    save( aux_file_name, '-append', 'right_wave_function', '-mat', '-v7.3');
                    creationMatrix = []; % clear creationMatrix from mat file
                    save( aux_file_name, '-append', 'creationMatrix', '-mat', '-v7.3');

                end
                clearvars destructionMatrix left_wave_function;
                destructionMatrix = [];
                left_wave_function = [];
                save( aux_file_name, '-append', 'destructionMatrix', 'left_wave_function', '-mat', '-v7.3');
            end

            fprintf('        Begin calculating on-diagonal elements at time  %s.\n', datestr(now,'yymmdd_HHMMSS'))
            i_site = 1;
                left_wave_function = (aux_file.groundState') * ...
                                                        creationOperator_parallel( noOfSites, noOfUp, noOfDn , i_site, 'up' )';
                save( aux_file_name,'-append',  'left_wave_function', '-mat', '-v7.3');
                clearvars left_wave_function;   
                j_site = 1;  
                    right_wave_function = creationOperator_parallel( noOfSites, noOfUp, noOfDn , j_site, 'up' ) * ...
                        aux_file.groundState;
                    save( aux_file_name, '-append', 'right_wave_function', '-mat', '-v7.3');
                    clearvars right_wave_function;
                    k_sum = 0;
                    for k_eigenValues = 1:NUM_OF_EIGEN_VALUES_UP %sum over k
                        expo_factor = exp( tau*( groundStateEnergy - eigenValues_up(k_eigenValues)));
                        i_total = aux_file.left_wave_function * aux_file.eigenVectors_up(:,k_eigenValues);
                        j_total = dot( aux_file.right_wave_function, conj(aux_file.eigenVectors_up(:,k_eigenValues)) );
                        k_sum = k_sum + expo_factor * i_total * j_total;        
                        clearvars expo_factor i_total j_total;
                    end
                    diagonal_elem_up = k_sum;

            spinUpGreenFunction = spinUpGreenFunction + spinUpGreenFunction';
            for i_diag = 1:noOfSites
                spinUpGreenFunction(i_diag, i_diag) = diagonal_elem_up;
            end


            clearvars creationMatrix right_wave_function k_sum;
            right_wave_function = [];
            save( aux_file_name, '-append', 'right_wave_function', '-mat', '-v7.3');
            creationMatrix = []; % clear creationMatrix from mat file
            save( aux_file_name, '-append', 'creationMatrix', '-mat', '-v7.3');

            clearvars diagonal_elem_up;
            save(output_files{t_tau}, '-append', 'spinUpGreenFunction', '-v7.3');
            clearvars spinUpGreenFunction;
        end
        clearvars t_tau;
        clearvars sizeSpacePlusOne eigenVectors_up eigenValues_up ;
        for i_f = 1:length(output_files)
            save(output_files{i_f},'-append', 'NUM_OF_EIGEN_VALUES_UP', '-v7.3');
        end
    end
%% SPIN DOWN:    
    if strcmp( sector, 'dn' ) || strcmp( sector, 'both' )
        fprintf('Begin diagonalizing spin-dn of secondHamiltonian at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        if NUM_OF_EIGEN_VALUES >= expanded_space_size_dn
            NUM_OF_EIGEN_VALUES_DN = expanded_space_size_dn - 1;
            fprintf('NUM_EIGEN_VALUES exceeds dimension of spin-down matrix. Now set to %d\n', NUM_OF_EIGEN_VALUES_DN)
        else
            NUM_OF_EIGEN_VALUES_DN = NUM_OF_EIGEN_VALUES;
        end

        [eigenVectors_dn, eigenValues_dn] = eigs( hubbardHamiltonian_parallel( t, U, noOfSites, noOfUp, noOfDn + 1 ), ...
                                                NUM_OF_EIGEN_VALUES_DN, 'sa', OPTS);
        eigenValues_dn = diag(eigenValues_dn);
        fprintf('Done with diagonalization at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        save( aux_file_name, '-append', 'eigenVectors_dn', 'eigenVectors_dn', '-mat', '-v7.3');
        for i_f = 1:length(output_files)
            save(output_files{i_f},'-append','eigenVectors_dn','eigenValues_dn', '-v7.3');            
        end     
        clearvars i_f eigenVectors_dn;

        fprintf('Begin spin down calculations at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        for t_tau = 1:length(output_files) % loop over taus

            tau = list_of_taus(t_tau);
            fprintf('Working on tau = %4.2f     at time %s\n', tau, datestr(now,'yymmdd_HHMMSS'))
            spinDnGreenFunction=zeros(noOfSites);    

            fprintf('        Begin calculating off-diagonal elements at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
            for i_site=1:noOfSites 
                fprintf('        Working on i = %3d     at time %s\n', i_site, datestr(now,'yymmdd_HHMMSS'))
                left_wave_function = (aux_file.groundState') * ...
                                                      creationOperator_parallel( noOfSites, noOfUp, noOfDn , i_site, 'dn' )';
                save( aux_file_name, '-append', 'left_wave_function', '-mat', '-v7.3');
                clearvars left_wave_function;

                for j_site=(i_site+1):noOfSites 
                    fprintf('        Working on j =     %3d at time %s\n', j_site, datestr(now,'yymmdd_HHMMSS'))
                    right_wave_function = creationOperator_parallel( noOfSites, noOfUp, noOfDn , j_site, 'dn' ) * ...
                                                        aux_file.groundState;
                    save( aux_file_name,'-append',  'right_wave_function', '-mat', '-v7.3');
                    clearvars right_wave_function;
                    k_sum = 0;
                    for k_eigenValues = 1:NUM_OF_EIGEN_VALUES_DN %sum over k
                        expo_factor = exp( tau*( groundStateEnergy - eigenValues_dn(k_eigenValues)));                
                        i_total = aux_file.left_wave_function * aux_file.eigenVectors_dn(:,k_eigenValues);
                        j_total = dot( aux_file.right_wave_function, conj(aux_file.eigenVectors_dn(:,k_eigenValues)) );
                        k_sum = k_sum + expo_factor * i_total * j_total;

                        clearvars expo_factor i_total j_total;
                    end
                    spinDnGreenFunction(i_site,j_site) = k_sum;
                    clearvars creationMatrix right_wave_function;
                    right_wave_function = [];
                    save( aux_file_name,'-append',  'right_wave_function', '-mat', '-v7.3');
                    creationMatrix = []; % clear creationMatrix from mat file
                    save( aux_file_name,'-append',  'creationMatrix', '-mat', '-v7.3');
                end
                clearvars destructionMatrix left_wave_function;
                destructionMatrix = [];
                left_wave_function = [];
                save( aux_file_name, '-append', 'destructionMatrix', 'left_wave_function', '-mat', '-v7.3');
            end

            fprintf('        Begin calculating on-diagonal elements at time  %s.\n', datestr(now,'yymmdd_HHMMSS'))
            i_site = 1;
                left_wave_function = (aux_file.groundState') * ...
                                                    creationOperator_parallel( noOfSites, noOfUp, noOfDn , i_site, 'dn' )';
                save( aux_file_name,'-append',  'left_wave_function', '-mat', '-v7.3');
                clearvars left_wave_function;      
                j_site = 1;  
                    right_wave_function = creationOperator_parallel( noOfSites, noOfUp, noOfDn , j_site, 'dn' ) * ...
                        aux_file.groundState;
                    save( aux_file_name, '-append', 'right_wave_function', '-mat', '-v7.3');
                    clearvars right_wave_function;
                    k_sum = 0;
                    for k_eigenValues = 1:NUM_OF_EIGEN_VALUES_DN %sum over k
                        expo_factor = exp( tau*( groundStateEnergy - eigenValues_dn(k_eigenValues)));
                        i_total = aux_file.left_wave_function * aux_file.eigenVectors_dn(:,k_eigenValues);
                        j_total = dot( aux_file.right_wave_function, conj(aux_file.eigenVectors_dn(:,k_eigenValues)) );
                        k_sum = k_sum + expo_factor * i_total * j_total;        
                        clearvars expo_factor i_total j_total;
                    end
                    diagonal_elem_dn = k_sum;

            spinDnGreenFunction = spinDnGreenFunction + spinDnGreenFunction';
            for i_diag = 1:noOfSites
                spinDnGreenFunction(i_diag, i_diag) = diagonal_elem_dn;
            end

            clearvars creationMatrix right_wave_function k_sum;
            right_wave_function = [];
            save( aux_file_name, '-append', 'right_wave_function', '-mat', '-v7.3');
            creationMatrix = []; % clear creationMatrix from mat file
            save( aux_file_name, '-append', 'creationMatrix', '-mat', '-v7.3');

            clearvars diagonal_elem_dn;
            save(output_files{t_tau}, '-append', 'spinDnGreenFunction', '-v7.3');
            clearvars spinDnGreenFunction;
        end
        clearvars t_tau;
        clearvars sizeSpacePlusOne eigenVectors_dn eigenValues_dn;
        for i_f = 1:length(output_files)
            save(output_files{i_f},'-append', 'NUM_OF_EIGEN_VALUES_DN', '-v7.3');
        end
    end
else
    error('Error: cannot apply creation operator when number of electrons = number of sites');
end

time=toc

if strcmp( need_profiling, 'Yes' )
    profile off;
end

for i_f = 1:length(output_files)
    tau = list_of_taus(i_f);
    save(output_files{i_f},'-append','noOfSites','noOfUp','noOfDn','U','tau','t','time', 'method', 'commit_number', '-v7.3');            
end     
fprintf('Finish all calculations at time %s.\n', datestr(now,'yymmdd_HHMMSS'))

if strcmp( need_profiling, 'Yes' )
    profsave(profile('info'), profile_directory_name);
end
end
