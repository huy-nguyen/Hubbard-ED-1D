function [ spinUpGreenFunction, spinDnGreenFunction ] = unequalTimeGF( t, U, tau, noOfSites, noOfUp, noOfDn )
% calculate the equal time GF by constructing separate matrices for the c_i and c_j^\dagger operators

format compact;
savedFileName=strcat('ED_',int2str(noOfSites),'_sites_',int2str(noOfUp),'u',int2str(noOfDn),'d_U_',num2str(U, '%4.2f'),'_tau_',num2str(tau, '%4.2f'),'_t_',num2str(t),' ',datestr(now,'_yymmdd_HHMMSS'),'.mat');
sprintf('Saved data file: %s', savedFileName)


TOTAL_UP_STATES=nchoosek(noOfSites,noOfUp);
TOTAL_DN_STATES=nchoosek(noOfSites,noOfDn);
TOTAL_ALL_STATES=TOTAL_UP_STATES*TOTAL_DN_STATES;

% sprintf('Total num of states = %d',TOTAL_ALL_STATES)

tic;

if (noOfUp < noOfSites) && (noOfDn < noOfSites)

    [totalHamiltonian] = hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn );
    [eigenVectors,eigenValues]=eigs(totalHamiltonian,1,'sa'); %ASSUMING THAT THE HAMILTONIAN IS REAL SYMMETRIC
%     eigenValues=diag(eigenValues);
    
    groundState=eigenVectors(:,1);
    spinUpGreenFunction=zeros(noOfSites);
    spinDnGreenFunction=zeros(noOfSites);
    
%     disp('Done with diagonalization')
    
    save(savedFileName,'eigenVectors','eigenValues','totalHamiltonian'); %save variables...    
    disp('Saved eigenVectors, eigenValues, totalHamiltonian'); % for debugging

    clearvars totalHamiltonian eigenVectors eigenValues; %...before clearing them to save memory
    
    
            
    
    sizeOriginalSpace=nchoosek(noOfSites,noOfUp)*nchoosek(noOfSites,noOfDn);
    
    % the Hamiltonian in the original space
    firstHamiltonian=hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn );
    expmFirstHamiltonian= expm( tau*eye(sizeOriginalSpace)*firstHamiltonian );
    
    % spin up:
    sizeSpacePlusOne=nchoosek(noOfSites,noOfUp+1)*nchoosek(noOfSites,noOfDn);    
    
    
    % the Hamiltonian in expanded space:
    secondHamiltonian=hubbardHamiltonian( t, U, noOfSites, noOfUp+1, noOfDn );
    expmSecondHamiltonian=expm( -tau*eye(sizeSpacePlusOne)*secondHamiltonian );
    
    for i=1:noOfSites
        
        % destruction operator:
        destructionMatrix=creationOperator( noOfSites, noOfUp, noOfDn , i, 'up' );
        destructionMatrix=destructionMatrix'; % need to take Hermitian conjugate to turn creation into destruction operator
            
        for j=1:noOfSites
            % creation operator
            creationMatrix=creationOperator( noOfSites, noOfUp, noOfDn , j, 'up' );
            
            
            % put them all together
            spinUpGreenFunction(i,j)=(groundState')*expmFirstHamiltonian*destructionMatrix*expmSecondHamiltonian*creationMatrix*groundState;
            
            clearvars creationMatrix
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
            spinDnGreenFunction(i,j)=(groundState')*expmFirstHamiltonian*destructionMatrix*expmSecondHamiltonian*creationMatrix*groundState;
            
            clearvars creationMatrix
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
