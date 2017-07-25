# ||AUM||
import scipy.io
import scipy as sp
import matplotlib.pyplot as plt
IntV = range(10, 1210, 10)

a = sp.load('avg_corr.npz')
rho = a['rho']
rho_orig = a['rho_orig']

ind = rho[:, 1] != 0

plt.plot(IntV, sp.mean(rho[ind, :], axis=0), label='synced wthin sub')
plt.plot(IntV, sp.mean(rho_orig[ind, :], axis=0), label='unsynced within sub')

plt.ylim(ymax=1, ymin=0)
plt.savefig('rho_sync_vs_len_same_sub3_v2.pdf')
plt.show()

a = sp.load('avg_corr_sub.npz')
rho = a['rho']
rho_orig = a['rho_orig']

ind = rho[:, 1] != 0

plt.plot(IntV, sp.mean(rho[ind, :], axis=0), label='synced across sub')
plt.plot(IntV, sp.mean(rho_orig[ind, :], axis=0), label='unsynced across sub')

plt.ylim(ymax=1, ymin=0)
plt.legend(loc=4)

plt.savefig('rho_sync_vs_len_sub_sub3_v2.pdf')

plt.show()
