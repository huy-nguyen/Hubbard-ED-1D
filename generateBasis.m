function [ combinedBasis, TOTAL_ALL_STATES,TOTAL_UP_STATES, TOTAL_DN_STATES, upStates, dnStates ] = generateBasis( NUM_SITES, NUM_UP, NUM_DN )
%CREATEBASIS Summary of this function goes here
%   Detailed explanation goes here

if (NUM_UP <= NUM_SITES) && (NUM_DN <= NUM_SITES)

[ TOTAL_ALL_STATES,TOTAL_UP_STATES, TOTAL_DN_STATES, temp_combinedBasis, temp_upStates, temp_dnStates, upStates, dnStates ] = helper_generateBasis( NUM_SITES, NUM_UP, NUM_DN );

% for compatibility:
upup = de2bi(temp_combinedBasis(:,2), 'left-msb');
dndn = de2bi(temp_combinedBasis(:,3), 'left-msb');
combinedBasis = horzcat( temp_combinedBasis(:,1), temp_combinedBasis(:,2),temp_combinedBasis(:,3), upup, dndn);

else
   error('Error: number of states greater than number of sites.'); 
end

end

function [ TOTAL_ALL_STATES,TOTAL_UP_STATES, TOTAL_DN_STATES, combinedBasis, upStates, dnStates, upupStates, dndnStates ] = helper_generateBasis( NUM_SITES, NUM_UP, NUM_DN )
% remove some of the output arguments after testing

TOTAL_UP_STATES=nchoosek(NUM_SITES,NUM_UP);
TOTAL_DN_STATES=nchoosek(NUM_SITES,NUM_DN);
TOTAL_ALL_STATES=TOTAL_UP_STATES*TOTAL_DN_STATES;


v_up = (dec2bin(0:(2^NUM_SITES)-1)=='1');
u_up =v_up(sum(v_up,2)==NUM_UP,:);
upStates(:,1) =  bin2dec(num2str(u_up));
% for compatibility:
upupStates = zeros( TOTAL_UP_STATES, 1 + NUM_SITES);
upupStates(:,1) =  bin2dec(num2str(u_up));
upupStates(:,2:end) = u_up;

v_dn = (dec2bin(0:(2^NUM_SITES)-1)=='1');
u_dn =v_dn(sum(v_dn,2)==NUM_DN,:);
dnStates(:,1) =  bin2dec(num2str(u_dn));
% for compatiblity:
dndnStates = zeros( TOTAL_DN_STATES, 1 + NUM_SITES);
dndnStates(:,1) =  bin2dec(num2str(u_dn));
dndnStates(:,2:end) = u_dn;

combinedBasis = repmat(dnStates, TOTAL_UP_STATES, 3);
combinedBasis(:,2) = reshape(repmat(upStates', TOTAL_DN_STATES, 1), TOTAL_ALL_STATES, []);
combinedBasis(:,1) = 1:TOTAL_ALL_STATES;
end