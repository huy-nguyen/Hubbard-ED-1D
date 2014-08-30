function [ spinUpGreenFunction, spinDnGreenFunction ] = unequalTimeGF( t, U, tau, noOfSites, noOfUp, noOfDn )
% calculate the equal time GF by constructing separate matrices for the c_i and c_j^\dagger operators

format compact;
savedFileName=strcat('ED_',int2str(noOfSites),'_sites_',int2str(noOfUp),'u',int2str(noOfDn),'d_U_',num2str(U, '%4.2f'),'_tau_',num2str(tau, '%4.2f'),'_t_',num2str(t),' ',datestr(now,'_yymmdd_HHMMSS'),'.mat')

tic;

if (noOfUp < noOfSites) && (noOfDn < noOfSites)
    spinUpGreenFunction=zeros(noOfSites);
    spinDnGreenFunction=zeros(noOfSites);

    % Replace the following lines:
%     [totalHamiltonian] = hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn );
%     [eigenVectors,eigenValues]=eigs(totalHamiltonian,1,'sa'); %ASSUMING THAT THE HAMILTONIAN IS REAL SYMMETRIC
%     groundState=eigenVectors(:,1);
%     groundStateEnergy = eigenValues(1);
    
    % with this line to save memory:
    [groundState,groundStateEnergy]=eigs( hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn ),...
                                               1,'sa'); %ASSUMING THAT THE HAMILTONIAN IS REAL SYMMETRIC
    
    % Replace the following lines:
%     save(savedFileName,'eigenVectors','eigenValues','totalHamiltonian'); %save variables...    
%     disp('Saved eigenVectors, eigenValues, totalHamiltonian'); % for debugging
%     clearvars totalHamiltonian eigenVectors eigenValues; %...before clearing them to save memory
    
    % with these lines:
    save(savedFileName,'groundState','groundStateEnergy', '-v7.3'); %save variables...    
    disp('Saved groundState groundStateEnergy'); % for debugging
       
%% SPIN UP:
    sizeSpacePlusOne=nchoosek(noOfSites,noOfUp+1)*nchoosek(noOfSites,noOfDn);    
    
    % the Hamiltonian in expanded space:
    % Replace the following lines:
%     secondHamiltonian=hubbardHamiltonian( t, U, noOfSites, noOfUp+1, noOfDn );
%     expmSecondHamiltonian=expm( -tau*speye(sizeSpacePlusOne)*secondHamiltonian );
    
    % with this line:
    expmSecondHamiltonian = expm( -tau * speye(sizeSpacePlusOne) * ...
                                   hubbardHamiltonian( t, U, noOfSites, noOfUp+1, noOfDn ) );
                                   
    
    for i=1:noOfSites        
        % destruction operator:
        destructionMatrix=creationOperator( noOfSites, noOfUp, noOfDn , i, 'up' )'; % need to take Hermitian conjugate to turn creation into destruction operator
            
        for j=1:noOfSites
            % replace the following lines:
%             creationMatrix=creationOperator( noOfSites, noOfUp, noOfDn , j, 'up' );            
%             right_wave_function = creationMatrix*groundState;
%             clearvars creationMatrix;
            % with this line:
            right_wave_function = creationOperator( noOfSites, noOfUp, noOfDn , j, 'up' ) * ...
                                                groundState;
            left_wave_function =  (groundState') * destructionMatrix; 
            % replace the following lines:
%             temp = left_wave_function * expmSecondHamiltonian * right_wave_function;
%             spinUpGreenFunction(i,j) = exp(tau*groundStateEnergy) * temp;
            % with this line:
            spinUpGreenFunction(i,j) = (left_wave_function * expmSecondHamiltonian * right_wave_function) * exp(tau*groundStateEnergy);
            
            clearvars left_wave_function right_wave_function;
        end
        
        clearvars destructionMatrix;
    end
    
    clearvars secondHamiltonian expmSecondHamiltonian;
    
%% SPIN DOWN:
    sizeSpacePlusOne=nchoosek(noOfSites,noOfUp)*nchoosek(noOfSites,noOfDn+1);    
    
    % the Hamiltonian in expanded space:
    % replace the following lines:
%     secondHamiltonian=hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn+1 );
%     expmSecondHamiltonian=expm( -tau*speye(sizeSpacePlusOne)*secondHamiltonian );
    % with this line:
    expmSecondHamiltonian = expm( -tau * speye(sizeSpacePlusOne) * ...
                                   hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn+1 ));
    
    for i=1:noOfSites        
        % destruction operator:
        destructionMatrix=creationOperator( noOfSites, noOfUp, noOfDn , i, 'dn' )'; % need to take Hermitian conjugate to turn creation into destruction operator
            
        for j=1:noOfSites
            % replace the following lines:
%             creationMatrix=creationOperator( noOfSites, noOfUp, noOfDn , j, 'dn' );     
%             right_wave_function = creationMatrix*groundState;
%             clearvars creationMatrix;
            % with this line:
            right_wave_function = creationOperator( noOfSites, noOfUp, noOfDn , j, 'dn' ) * ...
                                        groundState;
            
            left_wave_function =  (groundState') * destructionMatrix;
            % replace the following lines:
%             temp = left_wave_function * expmSecondHamiltonian * right_wave_function;
%             spinDnGreenFunction(i,j) = exp(tau*groundStateEnergy) * temp;
            % with this line:
            spinDnGreenFunction(i,j) = (left_wave_function * expmSecondHamiltonian * right_wave_function) * exp(tau*groundStateEnergy); 
            
            clearvars left_wave_function right_wave_function;
        end
        
        clearvars destructionMatrix;
    end    
    
    clearvars secondHamiltonian expmSecondHamiltonian;

else
    error('Error: cannot apply creation operator when number of electrons = number of sites');
end

time=toc
save(savedFileName,'-append','noOfSites','noOfUp','noOfDn','U','tau','t','time','-v7.3');
disp('saved noOfSites, noOfUp, noOfDn, U, tau, t, time');

save(savedFileName, '-append','spinUpGreenFunction', 'spinDnGreenFunction', '-v7.3');
disp('saved spinUpGreenFunction, spinDnGreenFunction');



end
