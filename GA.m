L = 30; N = 3; L0 = L/N; S = 1; VD = 1/5; tCircle = 30; wc = [0;0;0]; draw = 0;
tic
populationSize = 100; nGenes = 60; range = [7 18; 45 55; 1 1; 20 30];%[11 25; 10 100; 0 2; 10 50] [V,L,S,R0]
tournamentSize = 2; pTournament = 0.75; mutationProbability = 1/nGenes;
crossoverProbablility = 0.8;
nGenerations = 1;
 
population = InitializePopulation(populationSize, nGenes);
populationFitness = zeros(populationSize, 1);

%[M,varOrder,Sbs] = GAinitializeODE(L0,N,S,wc);

for j = 1:nGenerations
    for i = 1:populationSize
        chromosome = population(i,:);
        x = DecodeChromosome(chromosome, 4, range);
        %populationFitness(i) = EvaluateIndividual(M,x(1),L,N,x(2),x(3),VD,x(4),tCircle,varOrder,Sbs);
        populationFitness(i) = SimulationV1(x(1),x(2),N,x(3),x(4),VD,tCircle,wc,draw);
    end
    
    tempPopulation = population;
    
    for i=1:2:populationSize
        i1 = TournamentSelection(populationFitness, pTournament, tournamentSize);
        i2 = TournamentSelection(populationFitness, pTournament, tournamentSize);
        chromosome1 = population(i1,:);
        chromosome2 = population(i2,:);
        if rand < crossoverProbablility
            tempPopulation(i:i+1,:) = Cross(chromosome1, chromosome2);
        else
            tempPopulation(i:i+1,:) = [chromosome1; chromosome2];
        end
    end
    
    for i = 1:populationSize
        chromosome = tempPopulation(i,:);
        tempPopulation(i,:) = Mutate(chromosome, mutationProbability);
    end
                
    [bestIndividualFitness, bestIndividualIndex] = max(populationFitness);
    bestIndividual = population(bestIndividualIndex,:);
    tempPopulation = InsertBestIndividual(tempPopulation, bestIndividual, 1);
    population = tempPopulation;
end
disp('Minimum point')
disp(DecodeChromosome(bestIndividual, 4, range))
disp('Minimum function value')
disp(1/bestIndividualFitness)
toc