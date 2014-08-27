function [ spinUpGreenFunction, spinDnGreenFunction ] = unequalTimeGF_long_tau( t, U, tau, noOfSites, noOfUp, noOfDn, NUM_OF_EIGEN_VALUES )
% calculate the unequal time GF by using series expansion

OPTS.issym = 1;
OPTS.isreal = 1;

format compact;
savedFileName=strcat('ED_',int2str(noOfSites),'_sites_',int2str(noOfUp),'u',int2str(noOfDn),'d_U_',num2str(U, '%4.2f'),'_tau_',num2str(tau, '%4.2f'),'_t_',num2str(t),' ',datestr(now,'_yymmdd_HHMMSS'),'.mat')

tic;

if (noOfUp < noOfSites) && (noOfDn < noOfSites)
    spinUpGreenFunction=zeros(noOfSites);
    spinDnGreenFunction=zeros(noOfSites);
    
    % with this line to save memory:
    [groundState,groundStateEnergy]=eigs( hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn ),...
                                               1,'sa'); %ASSUMING THAT THE HAMILTONIAN IS REAL SYMMETRIC
    save(savedFileName,'groundState','groundStateEnergy'); %save variables...    
    disp('Saved groundState groundStateEnergy'); % for debugging
       
%% SPIN UP:
    sizeSpacePlusOne=nchoosek(noOfSites,noOfUp+1)*nchoosek(noOfSites,noOfDn);  % might not be needed. 
    [eigenVectors, eigenValues] = eigs( hubbardHamiltonian( t, U, noOfSites, noOfUp+1, noOfDn ), ...
                                            NUM_OF_EIGEN_VALUES, 'sa', OPTS);
    eigenValues = diag(eigenValues);
    for i_site=1:noOfSites        
        destructionMatrix=creationOperator( noOfSites, noOfUp, noOfDn , i_site, 'up' )'; 
        left_wave_function =  (groundState') * destructionMatrix; 
        for j_site=1:noOfSites % maybe consider j_site=i_site:noOfSites ?
            right_wave_function = creationOperator( noOfSites, noOfUp, noOfDn , j_site, 'up' ) * ...
                                                groundState;
            k_sum = 0;
            for k_eigenValues = 1:NUM_OF_EIGEN_VALUES %sum over k
                expo_factor = exp( tau*( groundStateEnergy - eigenValues(k_eigenValues)));                
                i_total = left_wave_function * eigenVectors(:,k_eigenValues);
                j_total = dot( right_wave_function, conj(eigenVectors(:,k_eigenValues)) );
                k_sum = k_sum + expo_factor * i_total * j_total;
                
                clearvars expo_factor i_total j_total;
            end
            spinUpGreenFunction(i_site,j_site) = k_sum;
            clearvars creationMatrix right_wave_function;
        end
        clearvars destructionMatrix left_wave_function;
    end
    
    
    clearvars sizeSpacePlusOne eigenVectors eigenValues;
%% SPIN DOWN:    
    sizeSpacePlusOne=nchoosek(noOfSites,noOfUp)*nchoosek(noOfSites,noOfDn + 1);  % might not be needed. 
    [eigenVectors, eigenValues] = eigs( hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn + 1 ), ...
                                            NUM_OF_EIGEN_VALUES, 'sa', OPTS);
    eigenValues = diag(eigenValues);
    for i_site=1:noOfSites        
        destructionMatrix=creationOperator( noOfSites, noOfUp, noOfDn , i_site, 'dn' )'; 
        left_wave_function =  (groundState') * destructionMatrix; 
        for j_site=1:noOfSites % maybe consider j_site=i_site:noOfSites ?
            right_wave_function = creationOperator( noOfSites, noOfUp, noOfDn , j_site, 'dn' ) * ...
                                                groundState;
            k_sum = 0;
            for k_eigenValues = 1:NUM_OF_EIGEN_VALUES %sum over k
                expo_factor = exp( tau*( groundStateEnergy - eigenValues(k_eigenValues)));                
                i_total = left_wave_function * eigenVectors(:,k_eigenValues);
                j_total = dot( right_wave_function, conj(eigenVectors(:,k_eigenValues)) );
                k_sum = k_sum + expo_factor * i_total * j_total;
                
                clearvars expo_factor i_total j_total;
            end
            spinDnGreenFunction(i_site,j_site) = k_sum;
            clearvars creationMatrix right_wave_function;
        end
        clearvars destructionMatrix left_wave_function;
    end
    
    
    clearvars sizeSpacePlusOne eigenVectors eigenValues;

else
    error('Error: cannot apply creation operator when number of electrons = number of sites');
end

time=toc
save(savedFileName,'-append','noOfSites','noOfUp','noOfDn','U','tau','t','time');
disp('saved noOfSites, noOfUp, noOfDn, U, tau, t, time');

save(savedFileName, '-append','spinUpGreenFunction', 'spinDnGreenFunction');
disp('saved spinUpGreenFunction, spinDnGreenFunction');



end
