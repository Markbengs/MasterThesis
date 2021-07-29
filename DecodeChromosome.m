function x = DecodeChromosome(chromosome, nVariables, range)
nGenesInVariable = length(chromosome)/nVariables;
denominator = 1-2^-nGenesInVariable;
exponent = -(1:nGenesInVariable);
x = zeros(nVariables,1);
for i = 1:nVariables
    variableIndices = (i-1)*nGenesInVariable+1:i*nGenesInVariable;
    x(i) = -range(i,1) + 2*range(i,1)/denominator*(chromosome(variableIndices)*2.^exponent.')+range(i,2);
end
end