#%AUM
# Shree Ganeshaya Namaha
import scipy as sp
import numpy as np
import os
from sklearn.model_selection import KFold

from sklearn.svm import SVC
import h5py
#from mvpa2.suite import *

data_dir = '/ImagePTE1/ajoshi/For_Anand_MICCAI/'
subid = '100307'
design = '/ImagePTE1/ajoshi/For_Anand_MICCAI/BD.mat'

fmri_file = os.path.join(data_dir, subid, 'tfMRI_LANGUAGE_LR.mat')
fmri_file = '/home/ajoshi/coding_ground/brainsync/src/dynamic_fmri/fmri_task_fitted_rest_wt.mat'

idNaN_file = os.path.join(data_dir, 'IdxNaN.mat')
f = h5py.File(idNaN_file, 'r')
idxNaNL = np.array(f['idxNaNL']).squeeze()
idxNaNR = np.array(f['idxNaNR']).squeeze()

idxNaN = np.concatenate((idxNaNL, idxNaNR))
f.close()

fmrif = h5py.File(fmri_file, 'r')
fmri_dataL = np.array(fmrif['dataL'])
fmri_dataR = np.array(fmrif['dataR'])

fmri_data = np.concatenate((fmri_dataL, fmri_dataR), axis=1)
fmrif.close()

f = h5py.File(design, 'r')
#print(f.keys())
lang = np.array(f['LANGUAGE']['LR']['X'])
#print(lang)
f.close()

lang = 1.0 * (lang > 0)
conditions = np.argmax(lang[2:, :], axis=0)
#print(conditions)

condition_mask = np.isin(conditions, [1, 2, 3, 4, 5, 6])
# We apply this mask in the sampe direction to restrict the
# classification to the face vs cat discrimination

fmri_masked = fmri_data[condition_mask, :]

conditions = conditions[condition_mask]

svc = SVC(kernel='linear')
#print(svc)
fmri_masked = fmri_masked[:, idxNaN == 0]
svc.fit(fmri_masked, conditions)

prediction = svc.predict(fmri_masked)
#print(prediction)

#print(LANGUAGE['LR'])
#print((prediction == conditions).sum() / float(len(conditions)))

svc.fit(fmri_masked[:-30], conditions[:-30])

prediction = svc.predict(fmri_masked[-30:])
#print((prediction == conditions[-30:]).sum() / float(len(conditions[-30:])))

from sklearn.model_selection import KFold

cv = KFold(n_splits=5)
cv_accuracy = []
# The "cv" object's split method can now accept data and create a
# generator which can yield the splits.
for train, test in cv.split(X=fmri_masked):
    conditions_masked = conditions[train]
    svc.fit(fmri_masked[train], conditions_masked)
    prediction = svc.predict(fmri_masked[test])
    cv_accuracy.append(
        (prediction == conditions[test]).sum() / float(len(conditions[test])))

cv_accuracy = np.array(cv_accuracy)
print('CV accuracy')
print(cv_accuracy)
print('******%s******' % subid)
print(
    'CV accuracy mean(std): %g (%g)' % (cv_accuracy.mean(), cv_accuracy.std()))
print('*****done*****')

#LANGUAGE.LR.X[:,:2]
