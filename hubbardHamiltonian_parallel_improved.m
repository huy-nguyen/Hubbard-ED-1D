function [ totalHamiltonian, kineticHamiltonian,  potentialHamiltonian] = hubbardHamiltonian_parallel_improved( t, U, noOfSites, noOfUp, noOfDn, NUM_CORES )

[ combinedBasis, totalNoOfPossiblestates,totalNoOfUpStates, totalNoOfDnStates, upStates, dnStates ] = generateBasis( noOfSites, noOfUp, noOfDn );

%% POTENTIAL HAMILTONIAN:
extracted_up_states_outside = combinedBasis(:,2); % for parfor
extracted_dn_states_outside = combinedBasis(:,3);
original_j = 1:totalNoOfPossiblestates;
j_assignment = splitvect(original_j, NUM_CORES);

potential_sparse_input = zeros(1, 3);

parfor core_counter_potential=1:NUM_CORES
    extracted_up_states = extracted_up_states_outside;
    extracted_dn_states = extracted_dn_states_outside;
    j_to_work_on = j_assignment{core_counter_potential};
    up_states_to_work_on = extracted_up_states(j_to_work_on);
    dn_states_to_work_on = extracted_dn_states(j_to_work_on);
    results_in_core_loop = zeros( length(j_to_work_on), 3);
    for j = 1:length(j_to_work_on)
        upSectorDec= up_states_to_work_on(j);
        dnSectorDec= dn_states_to_work_on(j);
        upSector = de2bi_modified(upSectorDec, noOfSites);
        dnSector= de2bi_modified(dnSectorDec, noOfSites); 
        doubleOccupancy=bitand(upSector,dnSector); % find all doubly-occupied sites
        potentialEnergy=sum(doubleOccupancy)*U; % sum up the number of doubly occupied sites and multiply by U
        results_in_core_loop(j, 3) = potentialEnergy;
        results_in_core_loop(j, 1) = j_to_work_on(j);
        results_in_core_loop(j, 2) = j_to_work_on(j);
    end
    potential_sparse_input = [potential_sparse_input; results_in_core_loop];
end

potentialHamiltonian = sparse( potential_sparse_input(2:end, 1), potential_sparse_input(2:end, 2), potential_sparse_input(2:end, 3),    totalNoOfPossiblestates, totalNoOfPossiblestates);
clearvars  extracted_up_states extracted_dn_states original j j_assisgnment upSectorDec dnSectorDec upSector dnSector;

%% KINETIC HAMILTONIAN:
% the number of electrons to be hopped over if the electrons hop around the lattice boundary (can easily see that this must be the case):
noOfUpInterior=noOfUp-1;
noOfDnInterior=noOfDn-1;
max_kinetic_num_non_zero_per_iteration = 2*max(1, (noOfSites - noOfUp)) + 2*max(1, (noOfSites - noOfDn) );
kinetic = zeros(1, 3);
original_m = 1:totalNoOfPossiblestates;
m_assignment = splitvect(original_m, NUM_CORES);
parfor core_counter_kinetic = 1:NUM_CORES % will be parfor    
    extracted_up_states = extracted_up_states_outside;
    extracted_dn_states = extracted_dn_states_outside;
    m_to_work_on = m_assignment{core_counter_kinetic};
    up_states_to_work_on = extracted_up_states( m_to_work_on);
    dn_states_to_work_on = extracted_dn_states( m_to_work_on);
    
    KINETIC_COUNTER = 0;
    kinetic_core_rows = zeros(ceil(max_kinetic_num_non_zero_per_iteration / NUM_CORES), 1);
    kinetic_core_cols = zeros( ceil( max_kinetic_num_non_zero_per_iteration / NUM_CORES), 1);
    kinetic_core_elems = zeros( ceil(max_kinetic_num_non_zero_per_iteration / NUM_CORES), 1);
    
    for internal_index = 1:length(m_to_work_on)    
%     for internal_index = m_to_work_on'
        m = m_to_work_on(internal_index);
        % save the unshifted spin up and spin down sectors:
        upSectorDec = up_states_to_work_on(internal_index);
        dnSectorDec = dn_states_to_work_on(internal_index);        
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
                   kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
                   kinetic_core_cols(KINETIC_COUNTER) = m;
                   kinetic_core_elems(KINETIC_COUNTER) = -t;
                   
               else % if the electron does hop around the boundary
                   if mod(noOfUpInterior,2)== 0 % if the number of electrons to be hopped over is even
                       KINETIC_COUNTER = KINETIC_COUNTER + 1;
                       kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
                       kinetic_core_cols(KINETIC_COUNTER) = m;
                       kinetic_core_elems(KINETIC_COUNTER) = -t;
                   else
                       KINETIC_COUNTER = KINETIC_COUNTER + 1;
                       kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
                       kinetic_core_cols(KINETIC_COUNTER) = m;
                       kinetic_core_elems(KINETIC_COUNTER) = +t;
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
                   kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
                   kinetic_core_cols(KINETIC_COUNTER) = m;
                   kinetic_core_elems(KINETIC_COUNTER) = -t;
               else
                   if mod(noOfUpInterior,2)== 0
                       KINETIC_COUNTER = KINETIC_COUNTER + 1;
                       kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
                       kinetic_core_cols(KINETIC_COUNTER) = m;
                       kinetic_core_elems(KINETIC_COUNTER) = -t;
                   else
                       KINETIC_COUNTER = KINETIC_COUNTER + 1;
                       kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
                       kinetic_core_cols(KINETIC_COUNTER) = m;
                       kinetic_core_elems(KINETIC_COUNTER) = +t;
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
                    kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
                    kinetic_core_cols(KINETIC_COUNTER) = m;
                    kinetic_core_elems(KINETIC_COUNTER) = -t;
                else % if the electron does hop around the boundary
                    if mod(noOfDnInterior,2)== 0 % if that number is even
                        KINETIC_COUNTER = KINETIC_COUNTER + 1;
                        kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
                        kinetic_core_cols(KINETIC_COUNTER) = m;
                        kinetic_core_elems(KINETIC_COUNTER) = -t;
                    else
                        KINETIC_COUNTER = KINETIC_COUNTER + 1;
                        kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
                        kinetic_core_cols(KINETIC_COUNTER) = m;
                        kinetic_core_elems(KINETIC_COUNTER) = +t;
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
                    kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
                    kinetic_core_cols(KINETIC_COUNTER) = m;
                    kinetic_core_elems(KINETIC_COUNTER) = -t;
                else % if the electron does hop around the boundary
                    if mod(noOfDnInterior,2)== 0 % if that number is even
                        KINETIC_COUNTER = KINETIC_COUNTER + 1;
                        kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
                        kinetic_core_cols(KINETIC_COUNTER) = m;
                        kinetic_core_elems(KINETIC_COUNTER) = -t;
                    else
                        KINETIC_COUNTER = KINETIC_COUNTER + 1;
                        kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
                        kinetic_core_cols(KINETIC_COUNTER) = m;
                        kinetic_core_elems(KINETIC_COUNTER) = +t;
                    end
                end
                
            end
        end
              
    end    
    kinetic_core_rows = kinetic_core_rows( kinetic_core_rows ~= 0);
    kinetic_core_cols = kinetic_core_cols( 1:length(kinetic_core_rows));
    kinetic_core_elems = kinetic_core_elems( 1:length(kinetic_core_rows));    
    kinetic_per_core = horzcat(kinetic_core_rows, kinetic_core_cols, kinetic_core_elems);      
    kinetic = [kinetic; kinetic_per_core];
end

kineticHamiltonian = sparse( kinetic(2:end,1), kinetic(2:end,2), kinetic(2:end,3), totalNoOfPossiblestates, totalNoOfPossiblestates);
%% TOTAL HAMILTONIAN:
totalHamiltonian=kineticHamiltonian+potentialHamiltonian;

end