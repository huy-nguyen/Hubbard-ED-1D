function d = bi2de_modified(b2)
n = size(b2,2); 
pow2vector = 2.^(0:1:(n-1));
d = fliplr(b2)*pow2vector';

