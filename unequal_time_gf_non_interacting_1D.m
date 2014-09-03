function [ upGF, dnGF ] = unequal_time_gf_non_interacting_1D( noOfSites, noOfUp, noOfDn, t, tau )
% exact unequal time GF for 1-D lattice with odd number of sites, spin up and spin down electrons
if mod(noOfSites,2) ~= 1 || mod(noOfUp,2)~= 1 || mod(noOfDn,2)~= 1
    disp('Wrong inputs');
end

upGF=zeros(noOfSites,noOfSites);
dnGF=zeros(noOfSites,noOfSites);
for l=1:noOfSites
    for j=1:noOfSites
        elem=0;
        for n=((noOfUp+1)/2):((noOfSites-1)/2)
            elem=elem+ 2/noOfSites*cos(2*pi*n/noOfSites*(j-l))*exp(2*tau*t*cos(2*pi*n/noOfSites));
        end
        upGF(l,j)=elem;
    end
end

for l=1:noOfSites
    for j=1:noOfSites
        elem=0;
        for n=((noOfDn+1)/2):((noOfSites-1)/2)
            elem=elem+ 2/noOfSites*cos(2*pi*n/noOfSites*(j-l))*exp(2*tau*t*cos(2*pi*n/noOfSites));
        end
        dnGF(l,j)=elem;
    end
end

end

