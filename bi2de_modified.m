function d = bi2de_modified(b)
n = size(b,2);
   b2 = b;
   b = b2(:,n:-1:1);
pow2vector = 2.^(0:1:(size(b,2)-1));
d = b*pow2vector';

