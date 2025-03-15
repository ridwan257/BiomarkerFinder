import warnings
import numpy as np
from itertools import combinations
from sklearn.metrics import f1_score, roc_auc_score
from random import shuffle
import rutils as rt
import numbers
from copy import deepcopy
import math



def fit_IIC(X, y):
    unique_classes = np.unique(y)
    mean_dict = {}
    sd_dict = {}
    for label in unique_classes:
        data = X[y == label]
        mean_dict[label] = data.mean(axis=0)
        sd_dict[label] = data.std(axis=0, ddof=1)

    weight_list = np.zeros(X.shape[1])
    for i, j in combinations(unique_classes, 2):
        mu_num = np.abs(mean_dict[i] - mean_dict[j])
        mu_denom = sd_dict[i] + sd_dict[j]
        mu = mu_num / mu_denom

        sigma_num = sd_dict[i]**2 + sd_dict[j]**2
        sigma_denom = sd_dict[i] * sd_dict[j]
        sigma = np.log(sigma_num) - np.log(sigma_denom)

        weight_list += mu+sigma
    
    return weight_list


class Agent:
    def __init__(self, total_feature, chromosome=None, rnd_engine=np.random.default_rng(None), make_valid=True):
        self.length = total_feature

        if chromosome is None:
            self.chr = rnd_engine.choice([True, False], size=self.length)
        else : 
            self.chr = chromosome.copy()

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

    @classmethod
    def crossover(cls, mate1, mate2, rnd_engine):
        mask = rnd_engine.random(size=mate1.length) < 0.5
        child1_chr = np.where(mask, mate1.chr, mate2.chr)
        child2_chr = np.where(~mask, mate1.chr, mate2.chr)

        child1 = cls(mate1.length, child1_chr, rnd_engine)
        child2 = cls(mate2.length, child2_chr, rnd_engine)

        return child1, child2

    def copy(self):
        agent =  Agent(self.length, chromosome=self.chr, make_valid=False)
        if hasattr(self, 'fitness_'):
            agent.fitness_ = self.fitness_
            agent.frecord_ = self.frecord_.copy()
        return agent

    def mutate(self, mrate, swap_rate, rnd_engine):
        ## bit flip mutation
        flip_mask = rnd_engine.random(size=self.length) < mrate
        self.chr[flip_mask] = ~self.chr[flip_mask]

        ## swap mutation
        if rnd_engine.random() < swap_rate:
            on_indices = np.where(self.chr)[0]
            off_indices = np.where(~self.chr)[0]
            if len(on_indices) > 0 and len(off_indices) > 0:
                idx_on = rnd_engine.choice(on_indices)
                idx_off = rnd_engine.choice(off_indices)
                self.chr[idx_on], self.chr[idx_off] = False, True

        # Ensuring VALID chromosome
        self.make_valid(rnd_engine)

    def calculate_boot_fitness(self, X, y, model, k, w1, w2, w3, n_jobs, seed=None):
        ## objective function 1
        fn1 = 1 - self.chr.sum() / self.length

        ## Bootstrap Score Collection
        if n_jobs == 1:
            boot_scores, oob_scores = rt.evaluateBootstrap(
                X, y, model, feature_indices=self.chr, n_boot=100, seed=seed, metric=rt.pr_auc, use_proba=True
            )
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


## ================================================================================
## -------------------------------- GENETIC MODEL ---------------------------------
## ================================================================================
class GeneticModel:
    def __init__(
            self, X, y, model, 
            X_test = None, y_test = None,
            w1=0.5, w2=0.5, w3 = 0, kappa = 20,
            total_pop = 10, evolve_gen = 100, elite_factor = 0.1,
            elite_pop_var = 6,
            exploration_rate = 0.3, tounament_size = (3, 5),
            mutation_rate=0.05, swap_mrate = 0.1,
            n_jobs=1,
            random_state = None
    ):  
        self.rnd_engine = np.random.default_rng(random_state)
        self.n_samples, self.n_features = X.shape
        self.X = X
        self.X_test = X_test
        self.y = y
        self.y_test = y_test
        self.model = model
        self.w1 = w1
        self.w2 = w2
        self.w3 = w3
        self.kappa = kappa
        self.n_pop = total_pop
        self.n_jobs = n_jobs 

        ## hyperparameters, h_ -> hyperparam
        self.h_exploitation = 1 - exploration_rate
        self.h_tsize = tounament_size # TOurnament Size
        self.h_swap_rate = swap_mrate
        

        ## elitism factor
        if isinstance(elite_factor, numbers.Number):
            self.h_elitism = elite_factor
            self._adaptive_elite = False
        else:
            self.h_elitism = elite_factor[0]
            self.h_elitism_final = elite_factor[1]
            self._adaptive_elite = True
        self.elite_var = elite_pop_var
        
        
        ## Mutation Rate
        if isinstance(mutation_rate, numbers.Number):
            self.mrate = mutation_rate
            self._adaptive_mutation = False
        else:
            self.mrate = mutation_rate[0]
            self.mrate_final = mutation_rate[1]
            self._adaptive_mutation = True

        ## Stop changing hyper parameters
        self.n_unstable_generation = evolve_gen
        
        self.best_global_score = 0
        self.best_global_pop = None
        self.generation = 0
        self.population_log = []

        ## this will created later
        # self.population_
        # self.parents_
        # self.fitness_
        # self.best_local_pop_
        # self.best_local_score_
        # self.records_
        self.init_popupation()


    def init_popupation(self):
        self.generation = 0
        self.population_ = [Agent(self.n_features, rnd_engine=self.rnd_engine) for _ in range(self.n_pop)]

    def get_weights(self):
        return self.w1, self.w2, self.w3

    def fitness_evaluation(self):
        self.fitness_ = np.empty(self.n_pop, dtype=np.float64)
        self.best_local_pop_ = None
        self.best_local_score_ = -np.inf

        w1, w2, w3 = self.get_weights()
        for i, pop in enumerate(self.population_):
            fitness_score = pop.calculate_boot_fitness(
                self.X, self.y,
                self.model,
                self.kappa, w1, w2, w3,
                self.n_jobs, 
                seed=self.rnd_engine.integers(0, 2**30)
            )
            # fitness_score = pop.calculate_loo_fitnes(
            #     self.X, self.y,
            #     self.model,
            #     w1, w2
            # )

            self.fitness_[i] = fitness_score

            # updating best local population
            if fitness_score > self.best_local_score_:
                self.best_local_score_ = fitness_score
                self.best_local_pop_ = pop

            # updating best global population
            if fitness_score > self.best_global_score:
                self.best_global_score = fitness_score
                self.best_global_pop = pop


    def pick_parents(self, verbose=0):
        parents = []
        current_gen = min(self.generation, self.n_unstable_generation)
        t_size = round(self.h_tsize[0] + self.h_tsize[1] * (current_gen / self.n_unstable_generation))

        if verbose==1: print(f'tournament_size - {t_size}')
        ## aplying tournament selection
        for _ in range(self.n_pop):
            candidates = self.rnd_engine.choice(self.n_pop, size=t_size, replace=False)
            
            winner = candidates[np.argmax(self.fitness_[candidates])]
            parents.append(winner)
        
        parents = sorted(parents, key=lambda x: self.fitness_[x], reverse=True)

        ## now making mating pairs
        if self.rnd_engine.random() < self.h_exploitation:
            self.parents_ = np.array([parents[::2], parents[1::2]]).T
            if verbose==1: print('Exploiting best agents...')
        else:
            mid = self.n_pop // 2
            self.parents_ = np.array(list(zip( parents[:mid], parents[mid:] )))
            if verbose==1: print('Explorating agents...')
    

    def get_mutation_rate(self):
        ## Adaptive nutation
        if self._adaptive_mutation:
            current_gen = min(self.generation, self.n_unstable_generation)
            mrate = self.mrate + (self.mrate_final - self.mrate) * (current_gen / self.n_unstable_generation)
            # print('in get_mutation_rate mrate - ', mrate)
            return mrate
        
        return self.mrate
        

    def get_elite_pop(self):
        if self._adaptive_elite:
            current_gen = min(self.generation, self.n_unstable_generation)
            elite_factor = self.h_elitism + (self.h_elitism_final - self.h_elitism) * (current_gen / self.n_unstable_generation)
            n_elites = round(elite_factor * self.n_pop)
        else:
            n_elites = round(self.h_elitism * self.n_pop)
        
        n_elites += self.rnd_engine.binomial(n=self.elite_var, p=0.5) - self.elite_var//2
        elites = sorted(self.population_, key=lambda pop: pop.fitness_, reverse=True)[:n_elites]
        return n_elites, elites


    def evolve(self, verbose=0):
        # keeping some elites
        n_elites, elites = self.get_elite_pop()
        mrate = self.get_mutation_rate()

        if(verbose==1):
            print('N elites - ', n_elites)
            print('Mutation rate - ', mrate)

        new_pop = []

        ## crossover and mutaiton
        for m1, m2 in self.parents_:
            if len(new_pop) > (self.n_pop - n_elites):
                break  
            
            child1, child2 = Agent.crossover(
                mate1=self.population_[m1], mate2=self.population_[m2], 
                rnd_engine=self.rnd_engine
            )

            ## Mutation
            ## mutate(self, mrate, swap_rate, rnd_engine)
            child1.mutate(mrate, self.h_swap_rate, self.rnd_engine)
            child2.mutate(mrate, self.h_swap_rate, self.rnd_engine)

            new_pop.extend([child1, child2])
        
        # shuffle(new_pop)
        self.population_ = elites + new_pop[:self.n_pop - len(elites)]
        self.generation += 1



    def feature_importance(self):
        feature_mat = np.stack([pop.chr * pop.fitness_ for pop in self.population_log], axis=0)
        # feature_freq_mat = np.stack([pop.chr for pop in self.population_log], axis=0)
        feature_vec = feature_mat.sum(axis=0) / len(self.population_log)

        return feature_vec

    def populate_next(self, verbose=0):
        self.population_log.append(deepcopy(self.population_))

        self.fitness_evaluation()
        self.pick_parents(verbose)
        self.evolve(verbose)


        

        return self.best_local_pop_
    




    

        


## =========================================================================================================
## =========================================================================================================
## ----------------------------------- OLD FUNCTION ARCHIVES
## =========================================================================================================

## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## frist implemented method for crossover
## more simple methodss
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# def crossover(self):
#     chromosome_list = [None] * self.n_pop

#     for i, (mate1, mate2) in enumerate(self.parents_):
#         b_index = self.rnd_engine.integers(0, self.n_features)

#         chr1 = np.concatenate(
#             [self.population_[mate1].chr[:b_index], self.population_[mate2].chr[b_index:]]
#         )
#         chr2 = np.concatenate(
#             [self.population_[mate2].chr[:b_index], self.population_[mate1].chr[b_index:]]
#         )

#         mut_loc1 = self.rnd_engine.random(self.n_features) < self.mrate
#         mut_loc2 = self.rnd_engine.random(self.n_features) < self.mrate
#         chr1[mut_loc1] = ~chr1[mut_loc1]
#         chr2[mut_loc2] = ~chr2[mut_loc2]

#         chromosome_list[2*i] = Agent(self.n_features, chr1, self.rnd_engine)
#         chromosome_list[2*i+1] = Agent(self.n_features, chr2, self.rnd_engine)

#     return chromosome_list


## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## frist implemented method for pick parents
## more simple methodss
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# def pick_parents(self):
#     candidate_index = np.argsort(self.fitness_)[self.n_pop//2:]
#     fitness_density = self.fitness_[candidate_index] / np.sum(self.fitness_[candidate_index])
#     # fitness_density = self.fitness_ / self.fitness_.sum()

#     # self.parents_ = self.rnd_engine.choice(
#     #     self.n_pop, size=(self.n_pop//2, 2), replace=True, p=fitness_density
#     # )
#     # self.parents_ = self.rnd_engine.choice(
#     #     candidate_index, size=(self.n_pop//2, 2), replace=True, p=fitness_density
#     # )
#     parents = self.rnd_engine.choice(candidate_index, size=self.n_pop, replace=True, p=fitness_density)
#     self.rnd_engine.shuffle(parents)
#     self.parents_ = parents.reshape(self.n_pop//2, 2)