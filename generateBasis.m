function [ combinedBasis, TOTAL_ALL_STATES,TOTAL_UP_STATES, TOTAL_DN_STATES, upStates, dnStates ] = generateBasis( NUM_SITES, NUM_UP, NUM_DN )
%CREATEBASIS Summary of this function goes here
%   Detailed explanation goes here

if (NUM_UP <= NUM_SITES) && (NUM_DN <= NUM_SITES)

TOTAL_UP_STATES=nchoosek(NUM_SITES,NUM_UP);
TOTAL_DN_STATES=nchoosek(NUM_SITES,NUM_DN);
TOTAL_ALL_STATES=TOTAL_UP_STATES*TOTAL_DN_STATES;


v_up = (dec2bin(0:(2^NUM_SITES)-1)=='1');
u_up =v_up(sum(v_up,2)==NUM_UP,:);
upStates(:,1) =  bin2dec(num2str(u_up));

v_dn = (dec2bin(0:(2^NUM_SITES)-1)=='1');
u_dn =v_dn(sum(v_dn,2)==NUM_DN,:);
dnStates(:,1) =  bin2dec(num2str(u_dn));

combinedBasis = repmat(dnStates, TOTAL_UP_STATES, 3);
combinedBasis(:,2) = reshape(repmat(upStates', TOTAL_DN_STATES, 1), TOTAL_ALL_STATES, []);
combinedBasis(:,1) = 1:TOTAL_ALL_STATES;

else
   error('Error: number of states greater than number of sites.'); 
end

end