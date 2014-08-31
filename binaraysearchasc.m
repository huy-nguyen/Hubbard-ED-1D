% % Binary search. 
% % Search 'sval' in sorted vector 'x', returns index of 'sval' in 'x'
% %  
% % INPUT:
% % x: vector of numeric values, x should already be sorted in ascending order
% %    (e.g. 2,7,20,...120)
% % sval: numeric value to be search in x
% %  
% % OUTPUT:
% % index: index of sval with respect to x. If sval is not found in x
% %        then index is empty.
       
function index = binaraysearchasc(x,sval)

index=[];
from=1;
to=length(x);

while from<=to
    mid = round((from + to)/2);    
    diff = x(mid)-sval;
    if diff==0
        index=mid;
        return
    elseif diff<0   % x(mid) < sval
        from=mid+1;
    else              % x(mid) > sval
        to=mid-1;			
    end
end

% % --------------------------------------------
% % Example code for Testing
% % x = sort(randint(1,1000,[0,400])); 
% % sval=12
% % index = binaraysearchasc(x,sval)
% % x(index)
% % --------------------------------------------

% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------


