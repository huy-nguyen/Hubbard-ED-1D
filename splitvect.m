function chuckCell = splitvect(v, n)
% Splits a vector into number of n chunks of the same size (if possible).
% In not possible the chunks are almost of equal size.
%
% based on http://code.activestate.com/recipes/425044/
chuckCell = {};
vectLength = numel(v);
splitsize = 1/n*vectLength;
for i = 1:n
%newVector(end + 1) =
idxs = [floor(round((i-1)*splitsize)):floor(round((i)*splitsize))-1]+1;
chuckCell{end + 1} = v(idxs);
end
assert(sum(cellfun(@numel, chuckCell)) == vectLength, ...
'More or less elements after split')