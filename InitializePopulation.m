function population = InitializePopulation(populationSize, nGenes)
population = binornd(1,1/2,populationSize, nGenes);
end