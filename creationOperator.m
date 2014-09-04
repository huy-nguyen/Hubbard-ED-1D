function [ operatorMatrix ] = creationOperator( NUM_SITES, NUM_UP, NUM_DN , CREATION_INDEX, spin ) 
%inputs are for the space that is  to be acted on

if ( strcmp(spin,'up')==1 && NUM_UP < NUM_SITES ) || (strcmp(spin,'dn')==1 && NUM_DN < NUM_SITES)
    
    if strcmp(spin,'up')==1
        [ SMALL_BASIS, TOTAL_STATES_SMALL_BASIS,dummy1,dummy2,dummy3,dummy4 ] = generateBasis( NUM_SITES, NUM_UP, NUM_DN ); % original basis
        [ BIG_BASIS, TOTAL_STATES_BIG_BASIS,TOTAL_UP_STATES_BIG_BASIS, TOTAL_DN_STATES_BIG_BASIS, UP_BIG_BASIS, DN_BIG_BASIS ] = generateBasis( NUM_SITES, NUM_UP + 1, NUM_DN ); %basis for space with one more particle
    elseif strcmp(spin,'dn')==1
        [ SMALL_BASIS, TOTAL_STATES_SMALL_BASIS,dummy1,dummy2,dummy3,dummy4 ] = generateBasis( NUM_SITES, NUM_UP, NUM_DN );
        [ BIG_BASIS, TOTAL_STATES_BIG_BASIS,TOTAL_UP_STATES_BIG_BASIS, TOTAL_DN_STATES_BIG_BASIS, UP_BIG_BASIS, DN_BIG_BASIS ] = generateBasis( NUM_SITES, NUM_UP, NUM_DN + 1 );
    else
        disp('Error');
    end
    
    if strcmp(spin,'up')==1
        %the indices to extract the correct up/down state from the combined basis table
        INDEX_AFFECTED_SECTOR = 2;
        INDEX_OTHER_SECTOR = 3;
    elseif strcmp(spin,'dn')
        INDEX_AFFECTED_SECTOR = 3;
        INDEX_OTHER_SECTOR = 2;
    else
        disp('Error');
    end
    
    max_num_non_zero_elems = TOTAL_STATES_SMALL_BASIS;
    COUNTER = 0;
    output_rows = zeros(max_num_non_zero_elems, 1);
    output_cols = zeros(max_num_non_zero_elems, 1);
    output_elems = zeros(max_num_non_zero_elems, 1);
    
    for basisCounter=1:TOTAL_STATES_SMALL_BASIS %loop through all the columns of the operator matrix
        %apply the operator to the small basis
        currentStateDec = SMALL_BASIS(basisCounter, INDEX_AFFECTED_SECTOR);
        otherSectorDec = SMALL_BASIS(basisCounter, INDEX_OTHER_SECTOR);
        currentState= de2bi_modified(currentStateDec, NUM_SITES);
        otherSector= de2bi_modified(otherSectorDec, NUM_SITES); 
        
        if currentState(CREATION_INDEX)==0 % only need to act if that site is unoccupied
            resultantState=currentState;
            resultantState(CREATION_INDEX)=1; % create a particle on the correct site
            numElectronsInBetween=sum(currentState(1:CREATION_INDEX-1)==ones(1,CREATION_INDEX-1)); % number of electrons in between the create and destroy indices
            numElectronsToCommuteThrough=numElectronsInBetween;
            
            if mod(numElectronsToCommuteThrough,2)==0 % if the creation operator has to commute through an even number of creation operators to get to the correct position
                coefficient=1; % then the resulting state doesn't change sign
            else
                coefficient=-1; %otherwise there is a sign change
            end
            
            %then look up the resulting state in the big basis
            if strcmp(spin,'up')==1
                upIndexOfResultantState=  binaraysearchasc(UP_BIG_BASIS, bi2de_modified(resultantState) );
                dnIndexOfResultantState= binaraysearchasc(DN_BIG_BASIS, bi2de_modified(otherSector) );
            else
                dnIndexOfResultantState=  binaraysearchasc(DN_BIG_BASIS, bi2de_modified(resultantState) );
                upIndexOfResultantState= binaraysearchasc(UP_BIG_BASIS, bi2de_modified(otherSector) );
                
            end
            
            basisIndexOfResultantState=(upIndexOfResultantState-1)*TOTAL_DN_STATES_BIG_BASIS+dnIndexOfResultantState;
            
            COUNTER = COUNTER + 1;
            output_rows(COUNTER) = basisIndexOfResultantState;
            output_cols(COUNTER) = basisCounter;
            output_elems(COUNTER) = coefficient;
        end
        
    end

    output_rows = output_rows( output_rows ~= 0);
    output_cols = output_cols( 1:length(output_rows));
    output_elems = output_elems( 1:length(output_rows));
    
    operatorMatrix = sparse( output_rows, output_cols, output_elems, TOTAL_STATES_BIG_BASIS, TOTAL_STATES_SMALL_BASIS);
else
    error('Error: cannot apply creation operator when number of electrons = number of sites');
end
end

