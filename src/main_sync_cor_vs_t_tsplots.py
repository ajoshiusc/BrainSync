# ||AUM||
import scipy.io
import scipy as sp
import matplotlib.pyplot as plt

import seaborn as sns
#sns.set(style="darkgrid")
sns.set_style("whitegrid", {'axes.grid' : True})


sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 1})
sns.set_style({'font.family':'serif', 'font.serif':'Times New Roman'})

IntV = range(10, 1210, 10)

a = sp.load('avg_corr.npz')
rho = a['rho']
rho_orig = a['rho_orig']

ind = rho[:, 1] != 0

p1=sns.tsplot(data=rho[ind, :], time=IntV, color=[0,0,1])
sns.tsplot(data=rho_orig[ind, :], time=IntV, color=[1,0,0])
sns.plt.ylim(-0.1, 1)
p1.figure.savefig('within_sub3.pdf')

p1.figure.clf()
sns.plt.close()

a = sp.load('avg_corr_sub.npz')
rho = a['rho']
rho_orig = a['rho_orig']
ind = rho[:, 1] != 0
p1=sns.tsplot(data=rho[ind, :], time=IntV, color=[0,0,1])
sns.tsplot(data=rho_orig[ind, :], time=IntV, color=[1,0,0])
sns.plt.ylim(-0.1, 1)
p1.figure.savefig('across_sub3.pdf')
