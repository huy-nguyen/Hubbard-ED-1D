function output_files = unequalTimeGF_long_tau( t, U, tau_start, tau_end, tau_step, noOfSites, noOfUp, noOfDn, NUM_OF_EIGEN_VALUES )
% calculate the unequal time GF by using series expansion

OPTS.issym = 1;
OPTS.isreal = 1;

expanded_space_size_up = nchoosek(noOfSites,noOfUp+1)*nchoosek(noOfSites,noOfDn);
expanded_space_size_dn = nchoosek(noOfSites,noOfUp)*nchoosek(noOfSites,noOfDn + 1);
format compact;

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
                
                
tic;

if (noOfUp < noOfSites) && (noOfDn < noOfSites)
    fprintf('Begin calculations at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
    spinUpGreenFunction=zeros(noOfSites);
    spinDnGreenFunction=zeros(noOfSites);
    
    % with this line to save memory:
    [groundState,groundStateEnergy]=eigs( hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn ),...
                                               1,'sa'); %ASSUMING THAT THE HAMILTONIAN IS REAL SYMMETRIC
                                           
    for i_f = 1:length(output_files)
        save(output_files{i_f},'groundState','groundStateEnergy');            
    end     
    clearvars i_f;
    
    % NEED TO PRINT OUT DATA FILES
       
%% SPIN UP:
    fprintf('Begin spin-up calculations at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
    if NUM_OF_EIGEN_VALUES >= expanded_space_size_up
        NUM_OF_EIGEN_VALUES_UP = expanded_space_size_up - 1;
        fprintf('NUM_EIGEN_VALUES exceeds dimension of spin-up matrix. Now set to %d\n', NUM_OF_EIGEN_VALUES_UP)
    else
        NUM_OF_EIGEN_VALUES_UP = NUM_OF_EIGEN_VALUES;
    end
    
    [eigenVectors_up, eigenValues_up] = eigs( hubbardHamiltonian( t, U, noOfSites, noOfUp+1, noOfDn ), ...
                                            NUM_OF_EIGEN_VALUES_UP, 'sa', OPTS);
    eigenValues_up = diag(eigenValues_up);
    
    for i_f = 1:length(output_files)
        save(output_files{i_f},'-append','eigenValues_up','eigenVectors_up');            
    end     
    clearvars i_f;
    
    for t_tau = 1:length(output_files) % loop over taus
        
        tau = list_of_taus(t_tau);
        spinUpGreenFunction=zeros(noOfSites);
        
        fprintf('Begin calculating off-diagonal elements at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        for i_site=1:noOfSites        
            fprintf('Working on i = %3d     at time %s\n', i_site, datestr(now,'yymmdd_HHMMSS'))
            destructionMatrix=creationOperator( noOfSites, noOfUp, noOfDn , i_site, 'up' )'; 
            left_wave_function =  (groundState') * destructionMatrix; 
            for j_site=(i_site+1):noOfSites % maybe consider j_site=i_site:noOfSites ?
                fprintf('Working on j =     %3d at time %s\n', j_site, datestr(now,'yymmdd_HHMMSS'))
                right_wave_function = creationOperator( noOfSites, noOfUp, noOfDn , j_site, 'up' ) * ...
                                                    groundState;
                k_sum = 0;
                for k_eigenValues = 1:NUM_OF_EIGEN_VALUES_UP %sum over k
                    expo_factor = exp( tau*( groundStateEnergy - eigenValues_up(k_eigenValues)));                
                    i_total = left_wave_function * eigenVectors_up(:,k_eigenValues);
                    j_total = dot( right_wave_function, conj(eigenVectors_up(:,k_eigenValues)) );
                    k_sum = k_sum + expo_factor * i_total * j_total;

                    clearvars expo_factor i_total j_total;
                end
                spinUpGreenFunction(i_site,j_site) = k_sum;
                clearvars creationMatrix right_wave_function k_sum;
            end
            clearvars destructionMatrix left_wave_function;
        end

        fprintf('Begin calculating on-diagonal elements at time  %s.\n', datestr(now,'yymmdd_HHMMSS'))
        i_site = 1;
            destructionMatrix=creationOperator( noOfSites, noOfUp, noOfDn , i_site, 'up' )';
            left_wave_function =  (groundState') * destructionMatrix;    
            j_site = 1;  
                right_wave_function = creationOperator( noOfSites, noOfUp, noOfDn , j_site, 'up' ) * ...
                    groundState;
                k_sum = 0;
                for k_eigenValues = 1:NUM_OF_EIGEN_VALUES_UP %sum over k
                    expo_factor = exp( tau*( groundStateEnergy - eigenValues_up(k_eigenValues)));
                    i_total = left_wave_function * eigenVectors_up(:,k_eigenValues);
                    j_total = dot( right_wave_function, conj(eigenVectors_up(:,k_eigenValues)) );
                    k_sum = k_sum + expo_factor * i_total * j_total;        
                    clearvars expo_factor i_total j_total;
                end
                diagonal_elem_up = k_sum;

        spinUpGreenFunction = spinUpGreenFunction + spinUpGreenFunction';
        for i_diag = 1:noOfSites
            spinUpGreenFunction(i_diag, i_diag) = diagonal_elem_up;
        end

        clearvars diagonal_elem_up;
        save(output_files{t_tau}, '-append', 'spinUpGreenFunction');
        clearvars spinUpGreenFunction;
    end
    clearvars t_tau;
    clearvars sizeSpacePlusOne eigenVectors_up eigenValues_up ;
    
%% SPIN DOWN:        
    fprintf('Begin spin-down calculations at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
    if NUM_OF_EIGEN_VALUES >= expanded_space_size_dn
        NUM_OF_EIGEN_VALUES_DN = expanded_space_size_dn - 1;
        fprintf('NUM_EIGEN_VALUES exceeds dimension of spin-down matrix. Now set to %d\n', NUM_OF_EIGEN_VALUES_DN)
    else
        NUM_OF_EIGEN_VALUES_DN = NUM_OF_EIGEN_VALUES;
    end

    [eigenVectors_dn, eigenValues_dn] = eigs( hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn + 1 ), ...
                                            NUM_OF_EIGEN_VALUES_DN, 'sa', OPTS);
    eigenValues_dn = diag(eigenValues_dn);
    
    for i_f = 1:length(output_files)
        save(output_files{i_f},'-append','eigenVectors_dn','eigenValues_dn');            
    end     
    clearvars i_f;

    for t_tau = 1:length(output_files) % loop over taus
        
        tau = list_of_taus(t_tau);
        spinDnGreenFunction=zeros(noOfSites);    
    
        fprintf('Begin calculating off-diagonal elements at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
        for i_site=1:noOfSites 
            fprintf('Working on i = %3d     at time %s\n', i_site, datestr(now,'yymmdd_HHMMSS'))
            destructionMatrix=creationOperator( noOfSites, noOfUp, noOfDn , i_site, 'dn' )'; 
            left_wave_function =  (groundState') * destructionMatrix; 
            for j_site=(i_site+1):noOfSites % maybe consider j_site=i_site:noOfSites ?
                fprintf('Working on j =     %3d at time %s\n', j_site, datestr(now,'yymmdd_HHMMSS'))
                right_wave_function = creationOperator( noOfSites, noOfUp, noOfDn , j_site, 'dn' ) * ...
                                                    groundState;
                k_sum = 0;
                for k_eigenValues = 1:NUM_OF_EIGEN_VALUES_DN %sum over k
                    expo_factor = exp( tau*( groundStateEnergy - eigenValues_dn(k_eigenValues)));                
                    i_total = left_wave_function * eigenVectors_dn(:,k_eigenValues);
                    j_total = dot( right_wave_function, conj(eigenVectors_dn(:,k_eigenValues)) );
                    k_sum = k_sum + expo_factor * i_total * j_total;

                    clearvars expo_factor i_total j_total;
                end
                spinDnGreenFunction(i_site,j_site) = k_sum;
                clearvars creationMatrix right_wave_function;
            end
            clearvars destructionMatrix left_wave_function;
        end

        fprintf('Begin calculating on-diagonal elements at time  %s.\n', datestr(now,'yymmdd_HHMMSS'))
        i_site = 1;
            destructionMatrix=creationOperator( noOfSites, noOfUp, noOfDn , i_site, 'dn' )';
            left_wave_function =  (groundState') * destructionMatrix;    
            j_site = 1;  
                right_wave_function = creationOperator( noOfSites, noOfUp, noOfDn , j_site, 'dn' ) * ...
                    groundState;
                k_sum = 0;
                for k_eigenValues = 1:NUM_OF_EIGEN_VALUES_DN %sum over k
                    expo_factor = exp( tau*( groundStateEnergy - eigenValues_dn(k_eigenValues)));
                    i_total = left_wave_function * eigenVectors_dn(:,k_eigenValues);
                    j_total = dot( right_wave_function, conj(eigenVectors_dn(:,k_eigenValues)) );
                    k_sum = k_sum + expo_factor * i_total * j_total;        
                    clearvars expo_factor i_total j_total;
                end
                diagonal_elem_dn = k_sum;

        spinDnGreenFunction = spinDnGreenFunction + spinDnGreenFunction';
        for i_diag = 1:noOfSites
            spinDnGreenFunction(i_diag, i_diag) = diagonal_elem_dn;
        end

        clearvars diagonal_elem_dn;
        save(output_files{t_tau}, '-append', 'spinDnGreenFunction');
        clearvars spinDnGreenFunction;
    end
    clearvars t_tau;
    clearvars sizeSpacePlusOne eigenVectors_dn eigenValues_dn;


else
    error('Error: cannot apply creation operator when number of electrons = number of sites');
end

time=toc
    
for i_f = 1:length(output_files)
    save(output_files{i_f},'-append','noOfSites','noOfUp','noOfDn','U','tau','t','time', 'NUM_OF_EIGEN_VALUES_UP', 'NUM_OF_EIGEN_VALUES_DN');            
end     
fprintf('Finish all calculations at time %s.\n', datestr(now,'yymmdd_HHMMSS'))


end
