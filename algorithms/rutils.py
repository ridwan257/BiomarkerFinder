import numpy as np
from collections import defaultdict
from sklearn.metrics import confusion_matrix, f1_score, precision_recall_curve, auc
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from joblib import Parallel, delayed, dump, load
import matplotlib.pyplot as plt
import numbers

def echo():
    print("Hellow Babe How are you?")

def iterLog(i):
    print(f'\rIteration - {i}', end='', flush=True)

def writeln(lst):
    for e in lst:
        print(e)

def rank(arr, ascending=False):
    if ascending:
        sorted_indices = np.argsort(arr)
    else:
        sorted_indices = np.argsort(arr)[::-1]

    rank_array = np.empty_like(sorted_indices)
    rank_array[sorted_indices] = np.arange(1, len(arr) + 1)

    return rank_array

## ----------------
def set_seed(seed=None):
    return np.random.default_rng(seed)

def pr_auc(y_true, y_proba):
    precision, recall, _ = precision_recall_curve(y_true, y_proba)
    pr_auc = auc(recall, precision)
    return pr_auc

## ========================================================================
## ----------------------------- SOME PLOTING -----------------------------
## ========================================================================
def errorPlot(score_array, axis=0, vline=None, ddof=1):
    n = score_array.shape[1] if axis==0 else score_array.shape[0]

    if vline is not None:
        if isinstance(vline, numbers.Number):
            plt.axvline(vline, ls='--', alpha=0.8, c='black')
        else:
            for v in vline:
                plt.axvline(v, ls='--', alpha=0.8, c='black')

    plt.errorbar(
        x = range(n),
        y = score_array.mean(axis=axis), yerr=score_array.std(axis=axis, ddof=ddof),
        fmt='o', ecolor='red', alpha=0.4, capsize=3
    )
    plt.plot(range(n), score_array.mean(axis=axis), c='blue')

## ========================================================================
## ----------------------- BOOTSTRAP RESAMPLING ----------------------
## ========================================================================
def bootstrap_resample(X, y, shuffle=True, rnd_engine=None, oob=True, max_iters=100):
    if rnd_engine is None:
        rnd_engine = np.random.default_rng()

    index_map = defaultdict(list)
    all_indices = set(range(len(y)))
    n_y = len(set(y)) 

    # creating map of  
    for i, k in enumerate(y):
        index_map[k].append(i)

    counter = 0
    while True:
        bootstrap_index_list = []
        for indices in index_map.values():
            bootstrap_index_list.extend(rnd_engine.choice(indices, size=len(indices), replace=True))

        oob_sample_list = list(all_indices - set(bootstrap_index_list))

        # breaking the loop
        if not oob or (len(oob_sample_list) > 2 and np.unique(y[oob_sample_list]).size == n_y):
            break

        # stoping infinite loop
        if counter > max_iters:
            print(f"Warning: Failed to create valid OOB sample after {max_iters} iteraions.")
            break  
        counter += 1

    ## shuffling the bootstrap samples
    if shuffle:
        rnd_engine.shuffle(bootstrap_index_list)

    # Create bootstrap and OOB sets
    Xboot, yboot = X[bootstrap_index_list].copy(), y[bootstrap_index_list].copy()
    Xoob, yoob = X[oob_sample_list].copy(), y[oob_sample_list].copy()

    return Xboot, Xoob, yboot, yoob



## ========================================================================
## ----------------------- BOOTSTRAP CROSSVALIDATION ----------------------
## ========================================================================
def __each_iteration(X, y, model, metric, use_proba, seeds, i):
    rnd = np.random.default_rng(seeds[i])
    
    # Bootstrap sampling
    X_boot, X_oob, y_boot, y_oob = bootstrap_resample(X, y, rnd_engine=rnd)

    scaler = StandardScaler()
    X_boot_std = scaler.fit_transform(X_boot)
    X_oob_std = scaler.transform(X_oob)

    # Train model
    model.fit(X_boot_std, y_boot)

    # Get predictions
    y_boot_pred = model.predict(X_boot_std)
    y_oob_pred = model.predict(X_oob_std)

    if use_proba:
        y_boot_pred_proba = model.predict_proba(X_boot_std)[:, 1]
        y_oob_pred_proba = model.predict_proba(X_oob_std)[:, 1]
        boot_score = metric(y_boot, y_boot_pred_proba)
        oob_score = metric(y_oob, y_oob_pred_proba)
    else:
        boot_score = metric(y_boot, y_boot_pred)
        oob_score = metric(y_oob, y_oob_pred)

    return boot_score, oob_score

def evaluateBootstrapParallel(X, y, feature_indices, model, n_boot=100, metric=f1_score, use_proba=False, seed=None, n_jobs=3):
    X_subset = X[:, feature_indices]

    rnd = np.random.default_rng(seed)
    random_seeds = rnd.integers(low=0, high=2**31, size=n_boot)

    results = Parallel(n_jobs=n_jobs, backend="loky")(delayed(__each_iteration)(
        X_subset, y, model, metric, use_proba, random_seeds, i
    ) for i in range(n_boot))

    boot_scores, oob_scores = zip(*results)
    return np.array(boot_scores), np.array(oob_scores)

    

    
def evaluateBootstrap(X, y, model, feature_indices=None, n_boot=100, metric=f1_score, use_proba=False, seed=None, verbose=0):
    rnd = np.random.default_rng(seed)
    if feature_indices is None:
        X_subset = X.copy()
    else:
        X_subset = X[:, feature_indices]
    
    boot_scores = np.zeros(n_boot)
    oob_scores = np.zeros(n_boot)

    for i in range(n_boot):
        X_boot, X_oob, y_boot, y_oob = bootstrap_resample(X_subset, y, rnd_engine=rnd)

        scaler = StandardScaler()
        X_boot_std = scaler.fit_transform(X_boot)
        X_oob_std = scaler.transform(X_oob)

        model.fit(X_boot_std, y_boot)

        # Get predictions (either class labels or probabilities)
        y_boot_pred = model.predict(X_boot_std)
        y_oob_pred = model.predict(X_oob_std)

        if use_proba:
            y_boot_pred_proba = model.predict_proba(X_boot_std)[:, 1]
            y_oob_pred_proba = model.predict_proba(X_oob_std)[:, 1]
            boot_scores[i] = metric(y_boot, y_boot_pred_proba)
            oob_scores[i] = metric(y_oob, y_oob_pred_proba)
        else:
            boot_scores[i] = metric(y_boot, y_boot_pred)
            oob_scores[i] = metric(y_oob, y_oob_pred)


        if verbose == 1:
            print(f'\rIteration - {i}', end='', flush=True)

        if verbose == 2:
            print(f"Iteration {i} - Confusion Matrices:")
            print("Train:")
            print(confusion_matrix(y_boot, y_boot_pred))
            print("Test:")
            print(confusion_matrix(y_oob, y_oob_pred))

    if verbose == 1: print()
    return boot_scores, oob_scores


## ========================================================================
## -------------------- LEAVE-ONE-OUT CROSSVALIDATION ---------------------
## ========================================================================
def evaluateLOO(X, y, feature_indices, model, metric=f1_score, use_proba=False, verbose=0):
    X_subset = X[:, feature_indices]

    loo = LeaveOneOut()
    y_true, y_pred = [], []
    
    for i, (train_idx, test_idx) in enumerate(loo.split(X)):

        X_train, X_test = X_subset[train_idx], X_subset[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        model.fit(X_train, y_train)

        # TRue
        y_true.append(y_test[0])
        # Predict
        if use_proba:
            y_pred.append(model.predict_proba(X_test)[:, 1])
        else:
            y_pred.append(model.predict(X_test))

    y_true = np.array(y_true).flatten()
    y_pred = np.array(y_pred).flatten()

    if verbose == 1:
        print(confusion_matrix(y_true, y_pred))
    if verbose == 2 and not use_proba:
        print("True -", y_true)
        print("Pred -", y_pred)
        print('Confusion matrix -')
        print(confusion_matrix(y_true, y_pred))
    
    score = metric(y_true, y_pred)

    return score


## ========================================================================
## ------------------------ FISHER UNDERSAMPLING --------------------------
## ========================================================================
import numpy as np

def undersample_fisher(Xmat, y, total=None):
    count_y = np.bincount(y)
    u_label = np.argmin(count_y)  # Minority class
    o_label = np.argmax(count_y)  # Majority class

    u_idx = np.where(y == u_label)[0]
    o_idx = np.where(y == o_label)[0]
    
    bucket = []

    if total is None:
        total = count_y[u_label]

    # Compute Fisher score
    mean_squre = (Xmat[u_idx].mean(axis=0) - Xmat[o_idx].mean(axis=0))**2
    var_sum = Xmat[u_idx].var(axis=0, ddof=1) + Xmat[o_idx].var(axis=0, ddof=1)
    fisher_ratio = mean_squre / var_sum

    while len(bucket) < total:
        for idx in u_idx:
            x1 = Xmat[idx]
            # Compute distances efficiently using NumPy
            dist_arr = np.array([(j, np.sum(fisher_ratio * (x1 - Xmat[j])**2)) for j in o_idx])
            dist_arr = dist_arr[dist_arr[:, 1].argsort()]  # Sort by distance
            
            for key, _ in dist_arr:
                if key not in bucket:
                    bucket.append(key)
                    break

            if len(bucket) >= total:
                break
    
    complete_idx = sorted(np.hstack([bucket, u_idx]).astype(np.int_))

    return Xmat[complete_idx].copy(), y[complete_idx].copy()