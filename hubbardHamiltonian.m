function [ totalHamiltonian, kineticHamiltonian,  potentialHamiltonian] = hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn )


[ combinedBasis, totalNoOfPossiblestates,totalNoOfUpStates, totalNoOfDnStates, upStates, dnStates ] = generateBasis( noOfSites, noOfUp, noOfDn );

max_kinetic_num_non_zero_elems = 2*(noOfSites - noOfUp - 2) + 2*(noOfSites - noOfDn - 2);
KINETIC_COUNTER = 0;
kinetic_rows = zeros(max_kinetic_num_non_zero_elems, 1);
kinetic_cols = zeros(max_kinetic_num_non_zero_elems, 1);
kinetic_elems = zeros(max_kinetic_num_non_zero_elems, 1);


potential_elems = zeros(totalNoOfPossiblestates, 1);

extracted_up_states = combinedBasis(:,2); % for parfor
extracted_dn_states = combinedBasis(:,3);
for j=1:totalNoOfPossiblestates 
       upSectorDec= extracted_up_states(j);
       dnSectorDec=extracted_dn_states(j);
       upSector = de2bi_modified(upSectorDec, noOfSites);
       dnSector= de2bi_modified(dnSectorDec, noOfSites);       
       doubleOccupancy=bitand(upSector,dnSector); % find all doubly-occupied sites
       potentialEnergy=sum(doubleOccupancy)*U; % sum up the number of doubly occupied sites and multiply by U
       potential_elems(j) = potentialEnergy;
end
potentialHamiltonian = spdiags(potential_elems, 0, totalNoOfPossiblestates, totalNoOfPossiblestates);
clearvars potential_elems extracted_up_states extracted_dn_states;

% the number of electrons to be hopped over if the electrons hop around the lattice boundary (can easily see that this must be the case):
noOfUpInterior=noOfUp-1;
noOfDnInterior=noOfDn-1;



for m=1:totalNoOfPossiblestates % go through each state in the basis:
    % save the unshifted spin up and spin down sectors:
    upSectorDec = combinedBasis(m, 2);
    dnSectorDec = combinedBasis(m, 3);
    
    upSector= de2bi_modified(upSectorDec, noOfSites);
    dnSector= de2bi_modified(dnSectorDec, noOfSites);      
    
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
           upIndexOfLeftShiftedResult=  binaraysearchasc(upStates, bi2de_modified(leftShiftResult) );
           dnIndexOfLeftShiftedResult=mod( m-1,totalNoOfDnStates)+1;
           basisIndexOfLeftShiftedResult=(upIndexOfLeftShiftedResult-1)*totalNoOfDnStates+dnIndexOfLeftShiftedResult;
           % update that state:
           if leftShiftedIndex < n % if the electron does not hop around the boundary
               KINETIC_COUNTER = KINETIC_COUNTER + 1;
               kinetic_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
               kinetic_cols(KINETIC_COUNTER) = m;
               kinetic_elems(KINETIC_COUNTER) = -t;
               
           else % if the electron does hop around the boundary               
               if mod(noOfUpInterior,2)== 0 % if the number of electrons to be hopped over is even
                   KINETIC_COUNTER = KINETIC_COUNTER + 1;
                   kinetic_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
                   kinetic_cols(KINETIC_COUNTER) = m;
                   kinetic_elems(KINETIC_COUNTER) = -t;
               else
                   KINETIC_COUNTER = KINETIC_COUNTER + 1;
                   kinetic_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
                   kinetic_cols(KINETIC_COUNTER) = m;
                   kinetic_elems(KINETIC_COUNTER) = +t;
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
           upIndexOfRightShiftedResult=  binaraysearchasc(upStates, bi2de_modified(rightShiftResult) );
           dnIndexOfRightShiftedResult=mod( m-1,totalNoOfDnStates)+1;
           basisIndexOfRightShiftedResult=(upIndexOfRightShiftedResult-1)*totalNoOfDnStates+dnIndexOfRightShiftedResult;
           % update that state:
           if rightShiftedIndex > n
               KINETIC_COUNTER = KINETIC_COUNTER + 1;
               kinetic_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
               kinetic_cols(KINETIC_COUNTER) = m;
               kinetic_elems(KINETIC_COUNTER) = -t;
           else
               if mod(noOfUpInterior,2)== 0 
                   KINETIC_COUNTER = KINETIC_COUNTER + 1;
                   kinetic_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
                   kinetic_cols(KINETIC_COUNTER) = m;
                   kinetic_elems(KINETIC_COUNTER) = -t;
               else
                   KINETIC_COUNTER = KINETIC_COUNTER + 1;
                   kinetic_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
                   kinetic_cols(KINETIC_COUNTER) = m;
                   kinetic_elems(KINETIC_COUNTER) = +t;
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
           dnIndexOfLeftShiftedResult=  binaraysearchasc(dnStates, bi2de_modified(leftShiftResult) );
           upIndexOfLeftShiftedResult= floor(( m - 1 )/totalNoOfDnStates)+1;
           basisIndexOfLeftShiftedResult=(upIndexOfLeftShiftedResult-1)*totalNoOfDnStates+dnIndexOfLeftShiftedResult;
           
           if leftShiftedIndex < p % if the electron does not hop around the boundary
               KINETIC_COUNTER = KINETIC_COUNTER + 1;
               kinetic_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
               kinetic_cols(KINETIC_COUNTER) = m;
               kinetic_elems(KINETIC_COUNTER) = -t;
           else % if the electron does hop around the boundary               
               if mod(noOfDnInterior,2)== 0 % if that number is even
                   KINETIC_COUNTER = KINETIC_COUNTER + 1;
                   kinetic_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
                   kinetic_cols(KINETIC_COUNTER) = m;
                   kinetic_elems(KINETIC_COUNTER) = -t;
               else
                   KINETIC_COUNTER = KINETIC_COUNTER + 1;
                   kinetic_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
                   kinetic_cols(KINETIC_COUNTER) = m;
                   kinetic_elems(KINETIC_COUNTER) = +t;
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
           dnIndexOfRightShiftedResult=  binaraysearchasc(dnStates, bi2de_modified(rightShiftResult) );
           upIndexOfLeftShiftedResult= floor(( m - 1 )/totalNoOfDnStates)+1;
           basisIndexOfRightShiftedResult=(upIndexOfLeftShiftedResult-1)*totalNoOfDnStates+dnIndexOfRightShiftedResult;
           
           if rightShiftedIndex > p % if the electron does not hop around the boundary
               KINETIC_COUNTER = KINETIC_COUNTER + 1;
               kinetic_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
               kinetic_cols(KINETIC_COUNTER) = m;
               kinetic_elems(KINETIC_COUNTER) = -t;
           else % if the electron does hop around the boundary               
               if mod(noOfDnInterior,2)== 0 % if that number is even
                   KINETIC_COUNTER = KINETIC_COUNTER + 1;
                   kinetic_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
                   kinetic_cols(KINETIC_COUNTER) = m;
                   kinetic_elems(KINETIC_COUNTER) = -t;
               else
                   KINETIC_COUNTER = KINETIC_COUNTER + 1;
                   kinetic_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
                   kinetic_cols(KINETIC_COUNTER) = m;
                   kinetic_elems(KINETIC_COUNTER) = +t;
                   
               end
           end
           
       end
    end
    
end

kinetic_rows = kinetic_rows( kinetic_rows ~= 0);
kinetic_cols = kinetic_cols( 1:length(kinetic_rows));
kinetic_elems = kinetic_elems( 1:length(kinetic_rows));

kineticHamiltonian = sparse( kinetic_rows, kinetic_cols, kinetic_elems, totalNoOfPossiblestates, totalNoOfPossiblestates);

totalHamiltonian=kineticHamiltonian+potentialHamiltonian;



end