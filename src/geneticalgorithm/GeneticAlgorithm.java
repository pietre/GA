package geneticalgorithm;

import java.util.Arrays;

public class GeneticAlgorithm
{
  
    /*
    * Definiowanie rozmiaru populacji, chromosomu, maksymalnej ilości generacji,
    * prawdopodobieństwa mutacji i krzyżowania oraz precyzji rozwiązania (czyli
    * jaką część dopasowania najlepszego osobnika stanowi średnie dopasowanie
    * całej populacji).
    */
    final static int POP_SIZE = 30;
    final static int CHROM_SIZE = 16;
    final static int MAX_GENERATIONS = 300;
    final static double MUTATION_PROB = 0.2;
    final static double CROSSOVER_PROB = 0.8;
    final static double SOLUTION_CONVERGENCE = 0.99;
    
    /*
    * Precyzja rozwiązania (ile genów chromosomu jest przeznaczonych na
    * rozwinięcie dziesiętne fenotypu, np.: 4 oznacza precyzję 2^(-4)) oraz
    * początek i koniec przedziału poszukiwań.
    */
    final static int SOLUTION_PRECISION = 12;
    final static double INTERVAL_ORIGIN = -6;
    final static double INTERVAL_END =
            Math.pow(2, (CHROM_SIZE - SOLUTION_PRECISION)) + INTERVAL_ORIGIN;
    
    /*
    * Parametry aktualnej populacji: ilość krzyżowań, mutacji, numer generacji,
    * średnie dopasowanie osobników, suma przystosowań całej populacji, indeksy
    * najlepszego i najgorszego osobnika w danej populacji oraz zmienna
    * przechowująca genotyp najlepszego osobnika populacji.
    */
    static int crossovers;
    static int mutations;
    static int generation;
    static double averageFitness;
    static double fitnessSum;
    static int worstIndividualIndex;
    static int bestIndividualIndex;
    static int[] bestIndividual = new int[CHROM_SIZE];
       
    /*
    * Utworzenie populacji rodziców oraz potomnej o określonych parametrach.
    */
    static int[][] population = new int[POP_SIZE][CHROM_SIZE];
    static int[][] newPopulation = new int[POP_SIZE][CHROM_SIZE];
    
    /*
    * Utworzenie tablicy przechowującej przystosowanie wszystkich osobników
    * oraz tablicy przechowującej wynik losowania ruletki proporcjonalnej.
    */
    static double[] populationFitness = new double[POP_SIZE];
    static int[] rouletteResult = new int[POP_SIZE];
    
    /*
    * Funkcja do losowania liczby całkowitej z określonego przedziału.
    */
    static int randomInt(int min, int max)
    {
        return (int)(Math.random() * ((max - min) + 1) + min);
    }
    
    /*
    * Inicjalizacja - wylosowanie początkowych wartości wszystkich genów.
    */
    static void initialization()
    {
        for(int j = 0; j < POP_SIZE; j++)
            for(int i = 0; i < CHROM_SIZE; i++)
                population[j][i] = randomInt(0, 1);
    }
    
    /*
    * Pobieranie fenotypu, czyli odkodowanie parametrów osobnika (genotypu).
    */
    static double getPhenotype(int individual[])
    {
        double phenotype = 0;
        int n = individual.length;
        for(int i = 0; i < n; i++)
            phenotype += individual[i] * Math.pow(2, n - i - 1 - SOLUTION_PRECISION);
        phenotype+= INTERVAL_ORIGIN;
        return phenotype;
    }
            
    /*
    * Obliczenie wartości funkcji przystosowania wybranego osobnika.
    */
    static double fitness(int individual[])
    {
        double phenotype = getPhenotype(individual);
        return fitnessFunction(phenotype); 
    }  
    
    /*
    * Funkcja przystosowania, znalezienie jej maksimum na danym przedziale jest
    * celem działania algorytmu genetycznego.
    */
    static double fitnessFunction(double a)
    {       
        double function = -2 * (a + 6.2) * (a + 2.4) * (a + 2.4) * (a - 5.7)
                * (a - 5.7) * (a - 10);
        return function;
    }

    /*
    * Ocena przystosowania całej populacji.
    */
    static void evaluation()
    {
        for(int i = 0; i < POP_SIZE; i++)
            {
                populationFitness[i] = fitness(population[i]);
                fitnessSum += populationFitness[i];
            }
        averageFitness = fitnessSum / POP_SIZE;
    }
    
    /*
    * Ruletka proporcjonalna - określenie rozkładu przystosowania całej populacji
    * (kolejne sumy przystosowań), znormalizowanie do sumy 1 oraz losowanie
    * osobników które ulegną reprodukcji.
    * Wartości w tablicy rozkładu przystosowania wskazują na koniec przedziału
    * w którym jest wybierany dany osobnik.
    */
    static void rouletteWheelSelection()
    {
        double[] normalizedFitness = new double[POP_SIZE];
        
        normalizedFitness[0] = populationFitness[0];
        for(int i = 1; i < POP_SIZE; i++)
            normalizedFitness[i] = normalizedFitness[i-1] + populationFitness[i];
        
        for(int i = 0; i < POP_SIZE; i++)
            normalizedFitness[i] /= normalizedFitness[POP_SIZE-1];
        
        for(int i = 0; i < POP_SIZE; i++)
        {
            double randomDouble = Math.random();
            for(int j = 0; normalizedFitness[j] < randomDouble; j++)
                rouletteResult[i] = j + 1;        
        }  
    }
    
    /*
    * Przeniesienie do nowej populacji osobników wylosowanych metodą ruletki
    * proporcjonalnej.
    */
    static void replicateChosenIndividuals()
    {
        for(int i = 0; i < POP_SIZE; i++)
            newPopulation[i] = population[rouletteResult[i]];
    }

    /*
    * Krzyżowanie - losowanie liczby z przedziału <0; 1) dla każdego osobnika
    * w nowej populacji. W przypadku wylosowania liczby mniejszej od wybranego
    * prawdopodobieństwa krzyżowania, dla danego osobnika nowej populacji
    * losowany jest partner z populacji macierzystej, z którym się krzyżuje.
    * Proces polega na losowaniu miejsca krzyżowania, tj. pozycji od której
    * partnerzy wymieniają się bitami (genami). Do nowej populacji przenoszony
    * jest osobnik wylosowany z ruletki z cześcią genów partnera.
    */
    static void crossover()
    {
        for(int i = 0; i < POP_SIZE; i++)
        {
            if(Math.random() <= CROSSOVER_PROB)
            { 
                int crossoverPartnerIndex = 0;
                do
                    crossoverPartnerIndex = randomInt(0, POP_SIZE - 1);
                while(crossoverPartnerIndex == rouletteResult[i]);

                int crossoverPoint = randomInt(0, CHROM_SIZE - 1);

                System.arraycopy(
                        population[crossoverPartnerIndex], crossoverPoint,
                        newPopulation[i], crossoverPoint,
                        CHROM_SIZE - crossoverPoint);

                crossovers++;
            }
        }
    }
    
    /*
    * Mutacja - losowanie liczby z przedziału <0; 1) dla każdego osobnika.
    * W przypadku wylosowania liczby mniejszej od wybranego prawdopodobieństwa
    * mutacji, losowana jest pozycja (locus) bitu (genu), który jest negowany.
    */
    static void mutation()
    {
        for(int i = 0; i < POP_SIZE; i++)
        {
            if(Math.random() < MUTATION_PROB)
            {
                int position = randomInt(0, CHROM_SIZE - 1);
                newPopulation[i][position] = 1 - newPopulation[i][position];
                mutations++;
            }
        }
    }
    
    /*
    * Drukowanie informacji na temat bieżącej populacji.
    */
    static void printPopulationInfo()
    {
        System.out.println("------------------------------------------\n"
                + "------------------------------------------\n"
                + "------------------------------------------\n"
                + "Generacja: " + generation
                + "\tKrzyżowań: " + crossovers
                + "\tMutacji: " + mutations);
                
        System.out.println("\n------------------------------------------"
                + "\nOsobniki:\t\t\tPrzystosowanie:");
        for(int i = 0; i < POP_SIZE; i++)
            System.out.println(Arrays.toString(population[i]) + "\t\t\t" + populationFitness[i]);
        
        System.out.println("\t\t\tŚrednie przystosowanie: " + averageFitness
                + "\nIndeks najlepszego osobnika:\t\t" + bestIndividualIndex
                + "\nPrzystosowanie najlepszego osobnika:\t" + fitness(bestIndividual));
    }
    
    /*
    * Resetowanie parametrów określających daną populację.
    */
    static void resetValues()
    {
        averageFitness = 0;
        crossovers = 0;
        mutations = 0;
        fitnessSum = 0;
        bestIndividualIndex = 0;
        worstIndividualIndex = 0;
    }
    
    /*
    * Znajdowanie indeksu najlepszego osobnika populacji oraz kopiowanie jego
    * genotypu.
    */
    static void findBestIndividual()
    {
        for(int i = 0; i < POP_SIZE; i++)
            if(populationFitness[i] > populationFitness[bestIndividualIndex])
                bestIndividualIndex = i;
        bestIndividual = Arrays.copyOf(population[bestIndividualIndex], CHROM_SIZE);
    }
    
    /*
    * Znajdowanie indeksu najgorszego osobnika populacji.
    */
    static void findWorstIndex()
    {
        for(int i = 0; i < POP_SIZE; i++)
            if(populationFitness[i] < populationFitness[worstIndividualIndex])
                worstIndividualIndex = i;    
    }
    
    
    
    public static void main(String[] args)
    {
        System.out.println("Znajdowanie maksimum funkcji:"
                + "\nf(x) = -2(x+6.2)(x+2.4)^2(x-5.7)^2(x-10)"
                + "\nw przedziale: <" + INTERVAL_ORIGIN
                + ", " + INTERVAL_END + ">");
        
        initialization();
        evaluation();
        findBestIndividual();
        printPopulationInfo();
        
        while(generation < MAX_GENERATIONS && averageFitness <= SOLUTION_CONVERGENCE * fitness(bestIndividual))
        {
            resetValues();
            findBestIndividual();
            rouletteWheelSelection();
            replicateChosenIndividuals();
            crossover();
            mutation();
            for (int i = 0; i < POP_SIZE; i++)
                population[i] = Arrays.copyOf(newPopulation[i], CHROM_SIZE);
            evaluation();
            findWorstIndex();
            population[worstIndividualIndex] = bestIndividual;
            populationFitness[worstIndividualIndex] = fitness(bestIndividual);
            generation++;
            printPopulationInfo();
        }

        System.out.println("Genotyp najlepszego osobnika:\t\t"
                + Arrays.toString(bestIndividual)
                + "\nFenotp najlepszego osobnika:\t\t"
                + getPhenotype(bestIndividual));
    }
}