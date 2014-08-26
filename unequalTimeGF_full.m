function [ spinUpGreenFunction, spinDnGreenFunction ] = unequalTimeGF_full( t, U, tau, noOfSites, noOfUp, noOfDn )
% diagonalizing the full matrix of secondHamiltonian instead of exponentiating a sparse matrix.

format compact;
savedFileName=strcat('ED_',int2str(noOfSites),'_sites_',int2str(noOfUp),'u',int2str(noOfDn),'d_U_',num2str(U, '%4.2f'),'_tau_',num2str(tau, '%4.2f'),'_t_',num2str(t),' ',datestr(now,'_yymmdd_HHMMSS'),'.mat')


TOTAL_UP_STATES=nchoosek(noOfSites,noOfUp);
TOTAL_DN_STATES=nchoosek(noOfSites,noOfDn);
TOTAL_ALL_STATES=TOTAL_UP_STATES*TOTAL_DN_STATES;

tic;

if (noOfUp < noOfSites) && (noOfDn < noOfSites)

    [totalHamiltonian] = hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn );
    [eigenVectors,eigenValues]=eigs(totalHamiltonian,1,'sa'); %ASSUMING THAT THE HAMILTONIAN IS REAL SYMMETRIC
%     eigenValues=diag(eigenValues);
    
    groundState=eigenVectors(:,1);
    groundStateEnergy = eigenValues(1);
    spinUpGreenFunction=zeros(noOfSites);
    spinDnGreenFunction=zeros(noOfSites);
    
    
    save(savedFileName,'eigenVectors','eigenValues','totalHamiltonian'); %save variables...    
    disp('Saved eigenVectors, eigenValues, totalHamiltonian'); % for debugging

    clearvars totalHamiltonian eigenVectors eigenValues; %...before clearing them to save memory
    
    sizeOriginalSpace=nchoosek(noOfSites,noOfUp)*nchoosek(noOfSites,noOfDn);
       
    % spin up:
    sizeSpacePlusOne=nchoosek(noOfSites,noOfUp+1)*nchoosek(noOfSites,noOfDn);    
    
    
    % the Hamiltonian in expanded space:
    secondHamiltonian=full(hubbardHamiltonian( t, U, noOfSites, noOfUp+1, noOfDn ));
    [eigenVectors,eigenValues]=eig(secondHamiltonian);
    clearvars secondHamiltonian;

    eigenValues=diag(eigenValues);
    middle_matrix = diag(exp((groundStateEnergy - eigenValues)*tau));
    
    for i=1:noOfSites
        
        % destruction operator:
        destructionMatrix=creationOperator( noOfSites, noOfUp, noOfDn , i, 'up' );
        destructionMatrix=destructionMatrix'; % need to take Hermitian conjugate to turn creation into destruction operator
            
        for j=1:noOfSites
            % creation operator
            creationMatrix=creationOperator( noOfSites, noOfUp, noOfDn , j, 'up' );
            
            
            % put them all together            
            right_temp = creationMatrix*groundState;
            right_wave_function = eigenVectors' * right_temp;
            clearvars creationMatrix right_temp;
            
            left_temp =  (groundState') * destructionMatrix;            
            left_wave_function = left_temp * eigenVectors;
            clearvars left_temp;
            
            spinUpGreenFunction(i,j) = left_wave_function * middle_matrix * right_wave_function;
            
        end
        
        clearvars destructionMatrix;
    end
    
    clearvars secondHamiltonian expmSecondHamiltonian;
    
    % spin dn:
    sizeSpacePlusOne=nchoosek(noOfSites,noOfUp)*nchoosek(noOfSites,noOfDn+1);    
    
    % the Hamiltonian in expanded space:
    secondHamiltonian=hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn+1 );
    expmSecondHamiltonian=expm( -tau*eye(sizeSpacePlusOne)*secondHamiltonian );
    
    for i=1:noOfSites
        
        % destruction operator:
        destructionMatrix=creationOperator( noOfSites, noOfUp, noOfDn , i, 'dn' );
        destructionMatrix=destructionMatrix'; % need to take Hermitian conjugate to turn creation into destruction operator
            
        for j=1:noOfSites
            % creation operator
            creationMatrix=creationOperator( noOfSites, noOfUp, noOfDn , j, 'dn' );            
          
            % put them all together            
            right_temp = creationMatrix*groundState;
            right_wave_function = eigenVectors' * right_temp;
            clearvars creationMatrix right_temp;
            
            left_temp =  (groundState') * destructionMatrix;            
            left_wave_function = left_temp * eigenVectors;
            clearvars left_temp;
            
            spinDnGreenFunction(i,j) = left_wave_function * middle_matrix * right_wave_function;
        end
        
        clearvars destructionMatrix;
    end    
    
    clearvars secondHamiltonian expmSecondHamiltonian;

else
    error('Error: cannot apply creation operator when number of electrons = number of sites');
end

time=toc
save(savedFileName,'-append','noOfSites','noOfUp','noOfDn','U','tau','t','time');
disp('saved noOfSites, noOfUp, noOfDn, U, tau, t, time');

save(savedFileName, '-append','spinUpGreenFunction', 'spinDnGreenFunction');
disp('saved spinUpGreenFunction, spinDnGreenFunction');



end
