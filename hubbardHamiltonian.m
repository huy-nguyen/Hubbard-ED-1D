function [ totalHamiltonian, kineticHamiltonian,  potentialHamiltonian] = hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn )
% Form the 1-D Hubbard Hamiltonian

[ combinedBasis, totalNoOfPossiblestates,totalNoOfUpStates, totalNoOfDnStates, upStates, dnStates ] = generateBasis( noOfSites, noOfUp, noOfDn );

% potentialHamiltonian=zeros(totalNoOfPossiblestates);
% kineticHamiltonian=zeros(totalNoOfPossiblestates);
% totalHamiltonian=zeros(totalNoOfPossiblestates);

noOfPar=noOfUp+noOfDn;


potentialHamiltonian=spalloc(totalNoOfPossiblestates,totalNoOfPossiblestates,noOfPar);
kineticHamiltonian=spalloc(totalNoOfPossiblestates,totalNoOfPossiblestates,2*noOfPar);
totalHamiltonian=spalloc(totalNoOfPossiblestates,totalNoOfPossiblestates,3*noOfPar);

% Form the potential Hamiltonian:
for j=1:totalNoOfPossiblestates % need to look at this again
%    for j=1:totalNoOfPossiblestates
       upSector= combinedBasis(j,4:noOfSites+3);
       dnSector=combinedBasis(j,noOfSites+4:end);
       doubleOccupancy=bitand(upSector,dnSector); % find all doubly-occupied sites
       potentialEnergy=sum(doubleOccupancy)*U; % sum up the number of doubly occupied sites and multiply by U
       potentialHamiltonian(j,j)=potentialEnergy;
%    end
end

% the number of electrons to be hopped over if the electrons hop around the lattice boundary (can easily see that this must be the case):
noOfUpInterior=noOfUp-1;
noOfDnInterior=noOfDn-1;

% site numbering: 1 2 3 4 5 ...
% Form the kinetic Hamiltonian:
for m=1:totalNoOfPossiblestates % go through each state in the basis:
    % save the unshifted spin up and spin down sectors:
    upSector= combinedBasis(m,4:noOfSites+3);
    dnSector=combinedBasis(m,noOfSites+4:end);
    % find the occupied lattice sites:    
    upNonZero=find(upSector);
    dnNonZero=find(dnSector);
    
    % shift for spin up:
    for n=upNonZero % for each occupied site
       % left shift:
       leftShiftResult=upSector;
       leftShiftedIndex=mod( n-2,noOfSites)+1;% figure out which site is the one to its left (periodic boundary condition)
       
       if upSector(leftShiftedIndex)~= 1
           % perform the shift:
           leftShiftResult(n)=0;           
           leftShiftResult(leftShiftedIndex)=1;
             
           % figure out where in the basis this shifted state is
           upIndexOfLeftShiftedResult=  find(upStates(:,1)== bi2de(leftShiftResult,'left-msb') );
           dnIndexOfLeftShiftedResult=mod( m-1,totalNoOfDnStates)+1;
           basisIndexOfLeftShiftedResult=(upIndexOfLeftShiftedResult-1)*totalNoOfDnStates+dnIndexOfLeftShiftedResult;
           % update that state:
           if leftShiftedIndex < n % if the electron does not hop around the boundary
               kineticHamiltonian(basisIndexOfLeftShiftedResult,m)= kineticHamiltonian(basisIndexOfLeftShiftedResult,m) - t;
           else % if the electron does hop around the boundary               
               if mod(noOfUpInterior,2)== 0 % if the number of electrons to be hopped over is even
                   kineticHamiltonian(basisIndexOfLeftShiftedResult,m)= kineticHamiltonian(basisIndexOfLeftShiftedResult,m) - t;
               else
                   kineticHamiltonian(basisIndexOfLeftShiftedResult,m)= kineticHamiltonian(basisIndexOfLeftShiftedResult,m) + t;
               end
           end
       end
       
       % right shift:
       rightShiftResult=upSector;
       rightShiftedIndex=mod( n,noOfSites)+1;     
       
       if upSector(rightShiftedIndex)~= 1
           % perform the shift:
           rightShiftResult(n)=0;           
           rightShiftResult(rightShiftedIndex)=1;
             
           % figure out where in the basis this shifted state is
           upIndexOfRightShiftedResult=  find(upStates(:,1)== bi2de(rightShiftResult,'left-msb') );
           dnIndexOfRightShiftedResult=mod( m-1,totalNoOfDnStates)+1;
           basisIndexOfRightShiftedResult=(upIndexOfRightShiftedResult-1)*totalNoOfDnStates+dnIndexOfRightShiftedResult;
           % update that state:
           if rightShiftedIndex > n
               kineticHamiltonian(basisIndexOfRightShiftedResult,m)= kineticHamiltonian(basisIndexOfRightShiftedResult,m) - t;
           else
               if mod(noOfUpInterior,2)== 0 
                   kineticHamiltonian(basisIndexOfRightShiftedResult,m)= kineticHamiltonian(basisIndexOfRightShiftedResult,m) - t;
               else
                   kineticHamiltonian(basisIndexOfRightShiftedResult,m)= kineticHamiltonian(basisIndexOfRightShiftedResult,m) + t;
               end
           end
       end
    end    
    % shift for spin down:
    for p=dnNonZero
       % left shift:
       leftShiftResult=dnSector;
       leftShiftedIndex=mod( p-2,noOfSites)+1;     
       
       if dnSector(leftShiftedIndex)~= 1
           % perform the shift:
           leftShiftResult(p)=0;           
           leftShiftResult(leftShiftedIndex)=1;
             
           % figure out where in the basis this shifted state is
           dnIndexOfLeftShiftedResult=  find(dnStates(:,1)== bi2de(leftShiftResult,'left-msb') );
           upIndexOfLeftShiftedResult= floor(( m - 1 )/totalNoOfDnStates)+1;
           basisIndexOfLeftShiftedResult=(upIndexOfLeftShiftedResult-1)*totalNoOfDnStates+dnIndexOfLeftShiftedResult;
           
           if leftShiftedIndex < p % if the electron does not hop around the boundary
               kineticHamiltonian(basisIndexOfLeftShiftedResult,m)= kineticHamiltonian(basisIndexOfLeftShiftedResult,m) - t;
           else % if the electron does hop around the boundary               
               if mod(noOfDnInterior,2)== 0 % if that number is even
                   kineticHamiltonian(basisIndexOfLeftShiftedResult,m)= kineticHamiltonian(basisIndexOfLeftShiftedResult,m) - t;
               else
                   kineticHamiltonian(basisIndexOfLeftShiftedResult,m)= kineticHamiltonian(basisIndexOfLeftShiftedResult,m) + t;
               end
           end             
       end
       % right shift:
       rightShiftResult=dnSector;
       rightShiftedIndex=mod( p,noOfSites)+1;     
       
       if dnSector(rightShiftedIndex)~= 1
           % perform the shift:
           rightShiftResult(p)=0;           
           rightShiftResult(rightShiftedIndex)=1;
           
           % figure out where in the basis this shifted state is
           dnIndexOfRightShiftedResult=  find(dnStates(:,1)== bi2de(rightShiftResult,'left-msb') );
           upIndexOfLeftShiftedResult= floor(( m - 1 )/totalNoOfDnStates)+1;
           basisIndexOfRightShiftedResult=(upIndexOfLeftShiftedResult-1)*totalNoOfDnStates+dnIndexOfRightShiftedResult;
           
           if rightShiftedIndex > p % if the electron does not hop around the boundary
               kineticHamiltonian(basisIndexOfRightShiftedResult,m)= kineticHamiltonian(basisIndexOfRightShiftedResult,m) - t;
           else % if the electron does hop around the boundary               
               if mod(noOfDnInterior,2)== 0 % if that number is even
                   kineticHamiltonian(basisIndexOfRightShiftedResult,m)= kineticHamiltonian(basisIndexOfRightShiftedResult,m) - t;
               else
                   kineticHamiltonian(basisIndexOfRightShiftedResult,m)= kineticHamiltonian(basisIndexOfRightShiftedResult,m) + t;
               end
           end
           
       end
    end
    
end

totalHamiltonian=kineticHamiltonian+potentialHamiltonian;


end

