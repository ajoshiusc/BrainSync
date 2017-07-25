import scipy as sp
"""
BrainSync: An Orthogonal Transform for Syncronize the subject fMRI data to the reference fMRI data
Authors: Anand A Joshi, Minqi Chong, Jian Li, Richard M. Leahy
Based on 
AA Joshi, M Chong, RM Leahy, BrainSync: An Orthogonal Transformation for Synchronization of fMRI Data Across Subjects, Proc. MICCAI 2017.


"""

def normalizeData(pre_signal):
    """
     normed_signal, mean_vector, std_vector = normalizeData(pre_signal)
     This function normalizes the input signal to have 0 mean and unit
     norm in time.
     pre_signal: Time x Original Vertices data
     normed_signal: Normalized (Time x Vertices) signal
     mean_vector: Vertices x 1 mean for each time series
     std_vector : Vertices x 1 std dev for each time series
    """
    if sp.any(sp.isnan(pre_signal)):
        print('there are NaNs in the data matrix, synchronization\
may not work')

    pre_signal[sp.isnan(pre_signal)] = 0
    mean_vector = sp.mean(pre_signal, axis=0, keepdims=True)
    normed_signal = pre_signal - mean_vector
    std_vector = sp.linalg.norm(normed_signal, axis=0)  # , keepdims=True)
    std_vector[std_vector == 0] = 1e-116
    normed_signal = normed_signal / std_vector

    return normed_signal, mean_vector, std_vector


def brainSync(X, Y):
    """
   Input:
       X - Time series of the reference data (Time x Vertex)
       Y - Time series of the subject data (Time x Vertex)

   Output:
       Y2 - Synced subject data with respect to reference data (Time x Vertex)
       R - The orthogonal rotation matrix (Time x Time)

   Please cite the following publication:
       AA Joshi, M Chong, RM Leahy, BrainSync: An Orthogonal Transformation
       for Synchronization of fMRI Data Across Subjects, Proc. MICCAI 2017,
       in press.
       """
    if X.shape[0] > X.shape[1]:
        print('The input is possibly transposed. Please check to make sure \
that the input is time x vertices!')

    C = sp.dot(X, Y.T)
    U, _, V = sp.linalg.svd(C)
    R = sp.dot(U, V)
    Y2 = sp.dot(R, Y)
    return Y2, R
