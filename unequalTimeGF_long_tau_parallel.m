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
                               

tic;
if strcmp( need_profiling, 'Yes' )
    profile -memory on;
end

if (noOfUp < noOfSites) && (noOfDn < noOfSites)
    fprintf('Begin diagonalizing firstHamiltonian at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
    firstHamiltonian = hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn );
    [groundState,groundStateEnergy]=eigs( firstHamiltonian,...
                                               1,'sa'); %ASSUMING THAT THE HAMILTONIAN IS REAL SYMMETRIC
                                              
    for i_f = 1:length(output_files)
        save(output_files{i_f},'groundState','groundStateEnergy', 'firstHamiltonian', '-v7.3');            
    end     
    clearvars i_f;
    clearvars firstHamiltonian;
    fprintf('Done with diagonalization at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
    
%% SPIN UP:    
    if strcmp( sector, 'up' ) || strcmp( sector, 'both' )
        
        spinUpGreenFunctionCellArray = cell(1, length(output_files) );
        for i = 1:length(output_files)
            spinUpGreenFunctionCellArray{i} = zeros(noOfSites);
        end
        fprintf('Begin diagonalizing spin-up of secondHamiltonian at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        if NUM_OF_EIGEN_VALUES >= expanded_space_size_up
            NUM_OF_EIGEN_VALUES_UP = expanded_space_size_up - 1;
            fprintf('NUM_EIGEN_VALUES exceeds dimension of spin-up matrix. Now set to %d\n', NUM_OF_EIGEN_VALUES_UP)
        else
            NUM_OF_EIGEN_VALUES_UP = NUM_OF_EIGEN_VALUES;
        end
        secondHamiltonianUp = hubbardHamiltonian( t, U, noOfSites, noOfUp+1, noOfDn );
        [eigenVectors_up, eigenValues_up] = eigs( secondHamiltonianUp, ...
                                                NUM_OF_EIGEN_VALUES_UP, 'sa', OPTS);
        eigenValues_up = diag(eigenValues_up);
        fprintf('Done with diagonalization at time %s.\n', datestr(now,'yymmdd_HHMMSS'))

        for i_f = 1:length(output_files)
            save(output_files{i_f},'-append','eigenValues_up','eigenVectors_up', 'secondHamiltonianUp', '-v7.3');            
        end     
        clearvars i_f secondHamiltonianUp;

        fprintf('Begin spin up calculations at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        
        fprintf('        Begin calculating off-diagonal elements at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        for i_site=1:noOfSites        
            fprintf('        Working on i = %3d     at time %s\n', i_site, datestr(now,'yymmdd_HHMMSS'))
            destructionMatrixUp = creationOperator( noOfSites, noOfUp, noOfDn , i_site, 'up' )';
            left_wave_function = (groundState') * ...
                                                        destructionMatrixUp;
            clearvars destructionMatrixUp;
            for j_site=(i_site+1):noOfSites 
                fprintf('        Working on j =     %3d at time %s\n', j_site, datestr(now,'yymmdd_HHMMSS'))
                creationMatrixUp = creationOperator( noOfSites, noOfUp, noOfDn , j_site, 'up' );
                right_wave_function =  creationMatrixUp * ...
                                                    groundState;
                clearvars creationMatrixUp;  
                i_sum = sum(bsxfun(@times, left_wave_function', eigenVectors_up));
                j_sum = sum(bsxfun(@times, right_wave_function, eigenVectors_up));
                i_sum_times_j_sum = i_sum.*j_sum;
                aaa = tau_start:tau_step:tau_end;
                bbb = groundStateEnergy - eigenValues_up;
                expo_factor = exp( bbb * aaa)';
                ress = sum(bsxfun(@times, i_sum_times_j_sum,expo_factor), 2);
                for i_f = 1:length(output_files)
                    spinUpGreenFunctionCellArray{i_f}(i_site, j_site) = ress(i_f);
                end
                clearvars right_wave_function
            end
            clearvars left_wave_function;
        end

        fprintf('        Begin calculating on-diagonal elements at time  %s.\n', datestr(now,'yymmdd_HHMMSS'))
        i_site = 1;
            destructionMatrixUp = creationOperator( noOfSites, noOfUp, noOfDn , i_site, 'up' )';
            left_wave_function = (groundState') * ...
                                                        destructionMatrixUp;
            clearvars destructionMatrixUp;        
            j_site = 1;
                creationMatrixUp = creationOperator( noOfSites, noOfUp, noOfDn , j_site, 'up' );
                right_wave_function =  creationMatrixUp * ...
                                                    groundState;
                clearvars creationMatrixUp;                
                
                i_sum = sum(bsxfun(@times, left_wave_function', eigenVectors_up));
                j_sum = sum(bsxfun(@times, right_wave_function, eigenVectors_up));
                i_sum_times_j_sum = i_sum.*j_sum;
                aaa = tau_start:tau_step:tau_end;
                bbb = groundStateEnergy - eigenValues_up;
                expo_factor = exp( bbb * aaa)';
                ress = sum(bsxfun(@times, i_sum_times_j_sum,expo_factor), 2);
                for i_f = 1:length(output_files)
                    spinUpGreenFunctionCellArray{i_f}(i_site, j_site) = ress(i_f);
                end
        clearvars destructionMatrixUp creationMatrixUp;
        
        for t_tau = 1:length(output_files)
            spinUpGreenFunction = spinUpGreenFunctionCellArray{t_tau};    
            element_at_1_1 = spinUpGreenFunction(1, 1);
            spinUpGreenFunction = spinUpGreenFunction + spinUpGreenFunction';
            for i_diag = 1:noOfSites
                spinUpGreenFunction(i_diag, i_diag) = element_at_1_1;
            end
            save(output_files{t_tau},'-append', 'spinUpGreenFunction', 'NUM_OF_EIGEN_VALUES_UP', '-v7.3');
        end                
    end
%% SPIN DOWN:      
    if strcmp( sector, 'dn' ) || strcmp( sector, 'both' )
        
        spinDnGreenFunctionCellArray = cell(1, length(output_files) );
        for i = 1:length(output_files)
            spinDnGreenFunctionCellArray{i} = zeros(noOfSites);
        end
        fprintf('Begin diagonalizing spin-down of secondHamiltonian at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        if NUM_OF_EIGEN_VALUES >= expanded_space_size_dn
            NUM_OF_EIGEN_VALUES_DN = expanded_space_size_dn - 1;
            fprintf('NUM_EIGEN_VALUES exceeds dimension of spin-down matrix. Now set to %d\n', NUM_OF_EIGEN_VALUES_DN)
        else
            NUM_OF_EIGEN_VALUES_DN = NUM_OF_EIGEN_VALUES;
        end
        secondHamiltonianDn = hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn + 1);
        [eigenVectors_dn, eigenValues_dn] = eigs( secondHamiltonianDn, ...
                                                NUM_OF_EIGEN_VALUES_DN, 'sa', OPTS);
        eigenValues_dn = diag(eigenValues_dn);
        fprintf('Done with diagonalization at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        
        for i_f = 1:length(output_files)
            save(output_files{i_f},'-append','eigenValues_dn','eigenVectors_dn', 'secondHamiltonianDn', '-v7.3');            
        end     
        clearvars i_f secondHamiltonianDn;

        fprintf('Begin spin down calculations at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        
        fprintf('        Begin calculating off-diagonal elements at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        for i_site=1:noOfSites        
            fprintf('        Working on i = %3d     at time %s\n', i_site, datestr(now,'yymmdd_HHMMSS'))
            destructionMatrixDn = creationOperator( noOfSites, noOfUp, noOfDn , i_site, 'dn' )';
            left_wave_function = (groundState') * ...
                                                        destructionMatrixDn;
            clearvars destructionMatrixDn;
            for j_site=(i_site+1):noOfSites 
                fprintf('        Working on j =     %3d at time %s\n', j_site, datestr(now,'yymmdd_HHMMSS'))
                creationMatrixDn = creationOperator( noOfSites, noOfUp, noOfDn , j_site, 'dn' );
                right_wave_function =  creationMatrixDn * ...
                                                    groundState;
                clearvars creationMatrixDn;  
                i_sum = sum(bsxfun(@times, left_wave_function', eigenVectors_dn));
                j_sum = sum(bsxfun(@times, right_wave_function, eigenVectors_dn));
                i_sum_times_j_sum = i_sum.*j_sum;
                aaa = tau_start:tau_step:tau_end;
                bbb = groundStateEnergy - eigenValues_dn;
                expo_factor = exp( bbb * aaa)';
                ress = sum(bsxfun(@times, i_sum_times_j_sum,expo_factor), 2);
                for i_f = 1:length(output_files)
                    spinDnGreenFunctionCellArray{i_f}(i_site, j_site) = ress(i_f);
                end
                clearvars right_wave_function
            end
            clearvars left_wave_function;
        end

        fprintf('        Begin calculating on-diagonal elements at time  %s.\n', datestr(now,'yymmdd_HHMMSS'))
        i_site = 1;
            destructionMatrixDn = creationOperator( noOfSites, noOfUp, noOfDn , i_site, 'dn' )';
            left_wave_function = (groundState') * ...
                                                        destructionMatrixDn;
            clearvars destructionMatrixDn;        
            j_site = 1;
                creationMatrixDn = creationOperator( noOfSites, noOfUp, noOfDn , j_site, 'dn' );
                right_wave_function =  creationMatrixDn * ...
                                                    groundState;
                clearvars creationMatrixDn;                
                
                i_sum = sum(bsxfun(@times, left_wave_function', eigenVectors_dn));
                j_sum = sum(bsxfun(@times, right_wave_function, eigenVectors_dn));
                i_sum_times_j_sum = i_sum.*j_sum;
                aaa = tau_start:tau_step:tau_end;
                bbb = groundStateEnergy - eigenValues_dn;
                expo_factor = exp( bbb * aaa)';
                ress = sum(bsxfun(@times, i_sum_times_j_sum,expo_factor), 2);
                for i_f = 1:length(output_files)
                    spinDnGreenFunctionCellArray{i_f}(i_site, j_site) = ress(i_f);
                end
        clearvars left_wave_function right_wave_function;
        
        for t_tau = 1:length(output_files)
            spinDnGreenFunction = spinDnGreenFunctionCellArray{t_tau};    
            element_at_1_1 = spinDnGreenFunction(1, 1);
            spinDnGreenFunction = spinDnGreenFunction + spinDnGreenFunction';
            for i_diag = 1:noOfSites
                spinDnGreenFunction(i_diag, i_diag) = element_at_1_1;
            end
            save(output_files{t_tau},'-append', 'spinDnGreenFunction', 'NUM_OF_EIGEN_VALUES_DN', '-v7.3');
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
