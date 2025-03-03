import numba
from numba.typed import Dict
from numba import types
import numpy as np


def build_distance_matrix(data):
    sq_sum = np.sum(data**2, axis=1).reshape(-1, 1)
    dist_sq_matrix = sq_sum + sq_sum.T - 2 * np.dot(data, data.T)
    dist_sq_matrix[dist_sq_matrix < 0] = 0
    distance_matrix = np.sqrt(dist_sq_matrix)
    
    return distance_matrix

@numba.jit
def count_classes_dict(y):
    class_dict = Dict.empty(
        key_type = types.int64,
        value_type = types.int64
    )

    for val in y:
        if val in class_dict:
            class_dict[val] += 1
        else:
            class_dict[val] = 1

    return class_dict

@numba.jit
def class_prob_dict(y):
    class_dict = Dict.empty(
        key_type = types.int64,
        value_type = types.float64
    )

    for key in y:
        if key in class_dict:
            class_dict[key] += 1.0
        else:
            class_dict[key] = 1.0
    
    for k in class_dict.keys():
        class_dict[k] = class_dict[k] / y.size

    return class_dict

# np.unique(y, return_counts=True) is more effiecient
# but it does not work inside jit
@numba.jit
def count_classes_vec(y):
    unique_classes = np.unique(y).astype(np.int64)
    class_counter = np.zeros(unique_classes.size, dtype=np.int64)
    for i in range(y.size):
        for j in range(unique_classes.size):
            if y[i] == unique_classes[j]:
                class_counter[j] += 1
                break
    return unique_classes, class_counter 


def build_NNmatrix(dist_mat, y, k=-1):
    unique_classes, class_counter = np.unique(y, return_counts=True)

    if np.any(class_counter <= k):
        print("Unique  Class:", unique_classes)
        print("Class Counter:", class_counter)
        raise ValueError("Neighbour size is greater than total number of class")

    

    total_sample = dist_mat.shape[0]
    if k == -1:
        nearest_matrix = np.full((total_sample, total_sample), -1, dtype=np.int64)
        class_indices = np.repeat(unique_classes, class_counter)
        indices_locator = np.roll(np.cumsum(class_counter), 1)
        indices_locator[0] = 0
    else:
        nearest_matrix = np.full((total_sample, k * unique_classes.size), -1, dtype=np.int64)
        class_indices = np.repeat(unique_classes, k)
        indices_locator = np.arange(0, k * unique_classes.size, k)

    for i in range(total_sample):
        dist_arr_from_i = np.argsort(dist_mat[i])

        t_idx = indices_locator.copy()
        for j in dist_arr_from_i:
            if i == j: continue

            for r in range(unique_classes.size):
                if k != -1 and t_idx[r] >= indices_locator[r] + k: continue
                if unique_classes[r] == y[j]:
                    nearest_matrix[ i, t_idx[r] ] = j
                    t_idx[r] += 1
                    break

    return class_indices, nearest_matrix

## ====================================================================================
## ====================================================================================
@numba.jit
def neighbor_finder(class_indices, NN_vector, label, k):
    hits = np.empty(k, dtype=np.int64)
    j = 0
    for i in range(class_indices.size):
        if class_indices[i] == label:
            hits[j] = NN_vector[i]
            j += 1
    return hits
## ====================================================================================
## ====================================================================================
@numba.jit
def reliefF(X, y, class_indices, kNN_matrix, epoch=-1):
    total_sample, total_feature = X.shape
    weight_list = np.zeros(total_feature)
    k = np.sum(class_indices == class_indices[0])
    m = total_sample if epoch==-1 else epoch
    cls_prob = class_prob_dict(y)

    diffA = np.zeros(total_feature)
    for j in range(total_feature):
        diffA[j] = np.max(X[:, j]) - np.min(X[:, j])

    for i in range(m):
        idx = i if epoch==-1 else np.random.randint(0, total_sample)
        sample = X[idx]
        label = y[idx]

        ## hit processing
        hits = neighbor_finder(class_indices, kNN_matrix[idx], label, k)
        near_hits = X[hits]
        diff_sample_hits = np.abs(near_hits - sample)
        diff_sample_hits /= diffA

        # distance factor consideration
        dist_coef = diff_sample_hits.sum(axis=1) ** 2
        dist_coef /= dist_coef.sum()
        dist_coef = dist_coef.reshape((-1, 1))
        diff_sample_hits *= dist_coef

        DIST_same = diff_sample_hits.sum(axis=0) / m

        ## miss processing
        DIST_miss = np.zeros(total_feature)
        for other, prob in cls_prob.items():
            if other == label: continue
            misses = neighbor_finder(class_indices, kNN_matrix[idx], other, k)
            prob_factor = prob / (1 - cls_prob[label])
            near_misses = X[misses]
            diff_sample_misses = np.abs(near_misses - sample)
            diff_sample_misses /= diffA

            # distance factor consideration
            dist_coef = diff_sample_misses.sum(axis=1) ** 2
            dist_coef /= dist_coef.sum()
            dist_coef = dist_coef.reshape((-1, 1))
            diff_sample_misses *= dist_coef

            DIST_miss += prob_factor * diff_sample_misses.sum(axis=0) / m
        
        weight_list += DIST_miss - DIST_same

    return weight_list




        







## ====================================================================================
## ====================================================================================
def __init_jit():
    x_toy = np.random.random((10, 3))
    y_toy = np.repeat([1, 2, 3], [4, 5, 3])
    dist_mat_toy = build_distance_matrix(x_toy)
    c_idx, knn_mat = build_NNmatrix(dist_mat_toy, y_toy, k=2)
    reliefF(x_toy, y_toy, c_idx, knn_mat, epoch=2)

__init_jit()

    
