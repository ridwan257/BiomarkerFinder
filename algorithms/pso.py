import numpy as np
import warnings
import rutils as rt
import math

from sklearn.metrics import f1_score


def sigmoid(arr):
    return 1 / (1+np.exp(-arr))

class Agent:
    def __init__(self, total_feature, chromosome=None, velocity=None, rnd_engine=np.random.default_rng(None), make_valid=True):
        self.length = total_feature
        
        if chromosome is None:
            self.chr = rnd_engine.choice([True, False], size=self.length)
            self.velocity = rnd_engine.random(self.length)
        else : 
            self.chr = chromosome.copy()
            self.velocity = velocity.copy()

        # Ensuring VALID Chr
        if make_valid:
            self.make_valid(rnd_engine)
    

    def make_valid(self, rnd_engine):
        while self.chr.sum() == 0:
            warnings.warn(
                "Model Failed. Atleast one feature should have. Initilizing Random Chromosome...",
                stacklevel=2
            )
            self.chr = rnd_engine.choice([True, False], size=self.length)


    def copy(self):
        agent =  Agent(self.length, chromosome=self.chr, velocity=self.velocity, make_valid=False)
        if hasattr(self, 'fitness_'):
            agent.fitness_ = self.fitness_
            agent.frecord_ = self.frecord_.copy()
        
        return agent
    


    def calculate_boot_fitness(self, X, y, model, k, w1, w2, w3, n_jobs):
        ## objective function 1
        fn1 = 1 - self.chr.sum() / self.length

        ## Bootstrap Score Collection
        if n_jobs == 1:
            boot_scores, oob_scores = rt.evaluateBootstrap(X, y, self.chr, model, n_boot=100, seed=None)
        else:
            boot_scores, oob_scores = rt.evaluateBootstrapParallel(X, y, self.chr, model, n_boot=100, n_jobs=n_jobs)

        fn2 = oob_scores.mean() # mean score
        fn3 = math.exp(-k * oob_scores.var()) # variance score


        self.frecord_ = {
            'f1' : round(fn1, 4),
            'f2' : round(fn2, 4),
            'f3' : round(fn3, 4),
            'var_oob' : round(oob_scores.var(), 4),
            'mean_boot' : round(boot_scores.mean(), 4),
            'var_boot' : round(boot_scores.var(), 4)
        }

        self.fitness_ = w1*fn1 + w2*fn2 + w3*fn3
        return self.fitness_


    
    def calculate_loo_fitnes(self, X, y, model, w1, w2):
        ## objective function 1
        ## minimun number of features
        fn1 = 1 - self.chr.sum() / self.length

        ## LEAVE ONE OUT CROSSVALIDATION
        fn2 = rt.evaluateLOO(X, y, self.chr, model)

        self.frecord_ = {
            'f1' : round(fn1, 4),
            'f2' : round(fn2, 4),
        }

        self.fitness_ = w1*fn1 + w2*fn2

        return self.fitness_



class NBPSOModel:
    def __init__(
        self, X, y, model, total_pop=60, 
        inertia=0.9, C1=2, C2=1.8, Vmin=-6, Vmax=6,
        w1=0.05, w2=0.85, w3=0.10,
        kappa = 20, n_jobs = 1,
        random_state = None
    ):  
        self.rnd = np.random.default_rng(random_state)
        self.X = X
        self.y = y
        self.model = model
        self.n_samples, self.n_features = X.shape
        self.w1 = w1
        self.w2 = w2
        self.w3 = w3
        self.n_agent = total_pop
        self.generation = 0

        self.inertia = inertia
        self.C1 = C1
        self.C2 = C2
        self.Vmin = Vmin
        self.Vmax = Vmax

        self.kappa = kappa
        self.n_jobs = n_jobs

        self.best_local_agent = None
        self.best_local_score = 0
        self.best_global_agent = None
        self.best_global_score = 0

        self.fitness_ = np.empty(self.n_agent, dtype=np.float64)
        self.init_population()


    def fitness_evaluation(self):
        self.best_local_pop_ = None
        self.best_local_score_ = 0

        # w1, w2, w3 = self.get_weights()
        for i, pop in enumerate(self.population_):
            fitness_score = pop.calculate_boot_fitness(
                self.X, self.y,
                self.model,
                self.kappa, self.w1, self.w2, self.w3,
                self.n_jobs
            )
            # fitness_score = pop.calculate_loo_fitnes(
            #     self.X, self.y,
            #     self.model,
            #     self.w1, self.w2
            # )
            

            self.fitness_[i] = fitness_score

            # updating best local population
            if fitness_score >= self.best_local_score:
                self.best_local_score = fitness_score
                self.best_local_agent = pop.copy()

            # updating best global population
            if fitness_score >= self.best_global_score:
                self.best_global_score = fitness_score
                self.best_global_agent = pop.copy()
    
    def init_population(self):
        self.generation = 0
        self.population_ = [Agent(self.n_features, rnd_engine=self.rnd) for _ in range(self.n_agent)]
        self.fitness_evaluation()
    
    
    def update(self):
        for agent in self.population_:
            new_velocity = self.inertia*agent.velocity + \
            self.C1*self.rnd.random()*( np.int_(self.best_local_agent.chr) - np.int_(agent.chr) ) + \
            self.C2*self.rnd.random()*( np.int_(self.best_global_agent.chr) - np.int_(agent.chr) )
            
            agent.velocity = np.clip(new_velocity, self.Vmin, self.Vmax)
            vel_transform = 1 / (1 + np.exp(-agent.velocity))


            agent.chr = self.rnd.random(self.n_features) < vel_transform

    def calculate_rdf_cyclic(self):
        fitness_arr2 = np.roll(self.fitness_, -1)
        diff = np.abs(fitness_arr2 - self.fitness_)
        return diff.sum() / self.fitness_.sum()
    
    def calculate_population_diversity(self):
        # Compute pairwise Jaccard dissimilarity between all agents
        diversity = []
        for i in range(len(self.population_)):
            for j in range(i+1, len(self.population_)):
                s1 = set(np.where(self.population_[i].chr == True)[0])
                s2 = set(np.where(self.population_[j].chr == True)[0])
                jaccard = 1 - len(s1 & s2) / len(s1 | s2)
                diversity.append(jaccard)
        return np.mean(diversity)
    
    def converge_halt(self):
        rdf_value = self.calculate_rdf_cyclic()
        print(f'rdf values — {rdf_value}')
        if rdf_value < self.rnd.random():
            print('Randomizing 1/5 population.')
            break_point = round(self.n_agent / 5)
            new_population = [Agent(self.n_features, rnd_engine=self.rnd) for _ in range(break_point)]            
            self.population_ = new_population + self.population_[break_point:]

    
    def populate_next(self):
        self.fitness_evaluation()
        self.update()
        # self.converge_halt()
        d = self.calculate_population_diversity()
        print(f"Diversity IDX - {d}")
        self.generation += 1

        return self.best_local_agent

                      
                             


    


    




