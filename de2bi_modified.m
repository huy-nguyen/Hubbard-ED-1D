function b = de2bi_modified(d, n)
   
% Perform conversion.
b=rem(floor(d*pow2(1-n:0)),2);
