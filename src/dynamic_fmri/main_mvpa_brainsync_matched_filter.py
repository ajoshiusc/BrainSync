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

conditions = lang[0]
#print(conditions)

condition_mask = np.isin(conditions, [0, 1])
# We apply this mask in the sampe direction to restrict the
# classification to the face vs cat discrimination

fmri_masked = fmri_data[condition_mask, :]

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
cv_error = []
# The "cv" object's split method can now accept data and create a
# generator which can yield the splits.
for train, test in cv.split(X=fmri_masked):
    conditions_masked = conditions[train]
    svc.fit(fmri_masked[train], conditions_masked)
    prediction = svc.predict(fmri_masked[test])
    cv_error.append(
        (prediction == conditions[test]).sum() / float(len(conditions[test])))

cv_error = np.array(cv_error)
print('CV error')
print(cv_error)
print('******%s******' % subid)
print('CV error mean(std): %g (%g)' % (cv_error.mean(), cv_error.std()))
print('*****done*****')

#LANGUAGE.LR.X[:,:2]
