function [ combinedBasis, TOTAL_ALL_STATES,TOTAL_UP_STATES, TOTAL_DN_STATES, upStates, dnStates ] = generateBasis( NUM_SITES, NUM_UP, NUM_DN )
%CREATEBASIS Summary of this function goes here
%   Detailed explanation goes here

if (NUM_UP <= NUM_SITES) && (NUM_DN <= NUM_SITES)
% form the basis:
TOTAL_UP_STATES=nchoosek(NUM_SITES,NUM_UP);
TOTAL_DN_STATES=nchoosek(NUM_SITES,NUM_DN);
TOTAL_ALL_STATES=TOTAL_UP_STATES*TOTAL_DN_STATES;

tempUp=zeros(1,NUM_SITES);
tempUp(1:NUM_UP)=ones(1,NUM_UP);
upStates=zeros(TOTAL_UP_STATES,NUM_SITES+1);
upStates(:,2:NUM_SITES+1)=unique(perms(tempUp),'rows'); % state represented as a binary number
upStates(:,1)=bi2de(upStates(:,2:NUM_SITES+1),'left-msb'); % first column stores the decimal equivalent of the binary state
upStates=sortrows(upStates);


tempDn=zeros(1,NUM_SITES);
tempDn(1:NUM_DN)=ones(1,NUM_DN);
dnStates=zeros(TOTAL_DN_STATES,NUM_SITES+1);
dnStates(:,2:NUM_SITES+1)=unique(perms(tempDn),'rows');
dnStates(:,1)=bi2de(dnStates(:,2:NUM_SITES+1),'left-msb');
dnStates=sortrows(dnStates);

combinedBasis=zeros(TOTAL_ALL_STATES,3+NUM_SITES*2);

for i =1:TOTAL_ALL_STATES    
    dnIndex=mod( i-1,TOTAL_DN_STATES)+1;
    upIndex=floor((i-1)/TOTAL_DN_STATES)+1;
    
    combinedBasis(i,1)=i;
    combinedBasis(i,2)=upStates(upIndex,1);
    combinedBasis(i,3)=dnStates(dnIndex,1);
    
    combinedBasis(i,4:NUM_SITES+3)=upStates(upIndex,2:end);
    combinedBasis(i,NUM_SITES+4:end)=dnStates(dnIndex,2:end);       
end

else
   error('Error: number of states greater than number of sites.'); 
end

end

