import matplotlib.pylab as plt
from astropy.table import Table, join

import numpy as np

from scipy import stats
bad_sn = set(open('badsn.txt').read().splitlines()[1:])

def clean_table(table):
    to_remove, = np.where(np.vectorize(lambda x: x in bad_sn)(table['name']))
    print len(to_remove)	
    table.remove_rows(to_remove)


def report(x, y):
    valid = np.logical_and(np.isfinite(x), np.isfinite(y))
    valid *= np.logical_and(~x.mask, ~y.mask)
    x = np.array(x[valid])
    y = np.array(y[valid])

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    print "Spearmanr:", stats.spearmanr(x, y, axis=0)
    print "slope:", slope, "Intercept", intercept, "R-value", r_value, "P-value", p_value
    return slope, intercept


t0_lmax = Table.read('t0_lmax_forbumptime.dat', format='ascii')
csp_bumpt = Table.read('csp_bump_times_final.csv', format='ascii') 
cfa_bumpt = Table.read('cfa_bump_times_final.csv', format='ascii')

csp_types = Table.read('CSP_types.txt', format='ascii')
cfa_types = Table.read('CfA_types.txt', format='ascii')
sbv_csp = Table.read('sbvfile_csp.dat', format='ascii')
sbv_cfa = Table.read('sbv_cfa_optical.txt', format='ascii')
rflam_csp = Table.read('rflam_csp_final.dat', format='ascii')
rflam_cfa = Table.read('rflam_cfa_final.dat', format='ascii')

csp_table = join(csp_bumpt, sbv_csp, join_type='inner', keys='name')
csp_table = join(csp_table, csp_types, join_type='inner', keys='name')

cfa_table = join(cfa_bumpt, sbv_cfa, join_type='inner', keys='name')
cfa_table = join(cfa_table, cfa_types, join_type='inner', keys='name')

csp_table = join(csp_table, rflam_csp, join_type='outer', keys='name')
print csp_table.keys()

model_table = Table.read('rflam_models_final.dat', format='ascii')
print model_table.keys()

#print csp_table
#print cfa_table

clean_table(csp_table)
clean_table(cfa_table)
clean_table(t0_lmax)

plt.figure('bump_sbv')
norm = csp_table['SNID'] =="Normal"
pec = csp_table['SNID'] =="91T"
bg = csp_table['SNID'] =='91bg'

N = cfa_table['Wang'] == "N"
bgs = cfa_table['Wang'] == "91bg"
T = cfa_table['Wang'] == "91T"

plt.errorbar(csp_table['t_bump'][norm], csp_table['sbv'][norm], xerr=csp_table['t_bump_err'][norm], yerr=csp_table['sbv_err'][norm], fmt='o', label='CSP_{normal}', color='orange')
plt.errorbar(csp_table['t_bump'][pec], csp_table['sbv'][pec], xerr=csp_table['t_bump_err'][pec], yerr=csp_table['sbv_err'][pec], fmt='o', fillstyle='none',alpha=0.6, label='CSP_{91T}', color='black')
plt.errorbar(csp_table['t_bump'][bg], csp_table['sbv'][bg], xerr=csp_table['t_bump_err'][bg], yerr=csp_table['sbv_err'][bg], fmt='o', label='CSP_{91bg}', alpha=0.4, fillstyle='none', color='black')
plt.errorbar(cfa_table['t_bump'][N], cfa_table['sbv'][N], xerr=cfa_table['t_bump_err'][N], yerr=cfa_table['sbv_err'][N], color='r', fmt='s', label='CfA', alpha=0.6)

plt.errorbar(cfa_table['t_bump'][bgs], cfa_table['sbv'][bgs], xerr=cfa_table['t_bump_err'][bgs], yerr=cfa_table['sbv_err'][bgs], color='black', fmt='s', label='CfA_{91bg}', alpha=0.3)
plt.errorbar(cfa_table['t_bump'][T], cfa_table['sbv'][T], xerr=cfa_table['t_bump_err'][T], yerr=cfa_table['sbv_err'][T], color='black', fmt='s', label='CfA_{91T}', alpha=0.6)
plt.xlabel(r' $t_{bump,r}$')
plt.ylabel('$s_{BV}$')
plt.legend(prop={'size': 12})

condition = csp_table['sbv'] < 0.6

print csp_table['name'][condition]

print('t_bump and sbv relationship:')
slope, intercept = report(csp_table['t_bump'], csp_table['sbv'])
x = np.linspace(5, 30)
y = slope * x + intercept
print slope, intercept
plt.plot(x, y, color='k', linestyle='-.')

plt.show()
select = np.logical_and(csp_table['t_bump'] < 20, csp_table['sbv'] > 0.8)
select2 = np.logical_and(cfa_table['t_bump'] < 20, cfa_table['sbv'] > 0.8)
print csp_table['name'][select].tolist()
print cfa_table['name'][select2].tolist()

plt.savefig('plots/bump_sbv.pdf')

csp_table_lmax = join(csp_table, t0_lmax, join_type='inner', keys='name')
print csp_table_lmax.keys()

plt.figure('bump_t0')
plt.errorbar(csp_table_lmax['t_bump'], csp_table_lmax['t0.'], xerr=csp_table_lmax['t_bump_err'], yerr=csp_table_lmax['e_t0'], fmt='o', label='CSP')
plt.xlabel(r'$t_{bump,r}$')
plt.ylabel('t0')
plt.legend()
#plt.savefig('plots/bump_t0.pdf')

plt.figure('bump_lmax')
plt.errorbar(csp_table_lmax['t_bump'], csp_table_lmax['Lmax'], xerr=csp_table_lmax['t_bump_err'],yerr=csp_table_lmax['e_Lmax'], fmt='o', label='CSP')
plt.xlabel(r'$t_{bump,r}$')
plt.ylabel('L_{max}')
plt.legend()
plt.savefig('plots/bump_lmax.pdf')

'''
plt.figure('bump_lmax_hist')
N = csp_table_lmax['t_bump'].shape[0]
x = []
y = []
for i in range(1000):
    x.extend(csp_table_lmax['t_bump'] + csp_table_lmax['t_bump_err'] * np.random.randn(N))
    y.extend(csp_table_lmax['Lmax'] + csp_table_lmax['e_Lmax'] * np.random.randn(N))
plt.hexbin(x, y, bins=300, norm=mpl.colors.LogNorm())
plt.xlabel(r'$t_{bump,r}$')
plt.ylabel('Lmax')
plt.legend()
'''
#plt.savefig('plots/bump_lmax_hist.pdf')


plt.figure('sbv_lmax')
plt.errorbar(csp_table_lmax['Lmax'], csp_table_lmax['sbv'], yerr=csp_table_lmax['sbv_err'],xerr=csp_table_lmax['e_Lmax'], fmt='o', label='CSP')
plt.ylabel('$s_{BV}$')
plt.xlabel('$L_{max}$')
plt.legend()
plt.savefig('plots/sbv_lmax.pdf')

plt.figure('flux_t0')
plt.errorbar(csp_table_lmax['rflam'], csp_table_lmax['t0.'],yerr=csp_table_lmax['e_t0'], fmt='o', label='CSP')
plt.xlabel(r'$\langle f\rangle_{r,15-40}$')
plt.ylabel('t0')
plt.legend()

print 'rflam t0'

slope, intercept = report(csp_table_lmax['rflam'], csp_table_lmax['t0.'])

#Add plotting over the straight line.

x = np.linspace(0.18, 0.5)
y = slope * x + intercept
print slope, intercept
plt.plot(x, y, color='k', linestyle='-.')
plt.text(0.5, 0.5,'Slope')


plt.savefig('plots/flux_t0.pdf')

plt.figure('hist_rflam_ptf')


fluxes = np.genfromtxt('rflam_ptf_final.dat', delimiter=' ')[:, 1]
t0_ptf = slope * fluxes + intercept
plt.hist(t0_ptf, histtype='step', bins='auto', linewidth=2.0, color='black')
plt.xlabel('t0')
plt.ylabel('\# SNe')
plt.savefig('plots/hist_rflam_ptf.pdf')


names = np.genfromtxt('rflam_ptf_final.dat', delimiter=' ', dtype=None)
selected_names = [n[0] for n in names[t0_ptf < 26]]
selected_names.sort()
print selected_names


'''
bic = []
aic = []
component_list = range(1, 6)
for components in component_list:
    gmm = GaussianMixture(components)
    gmm.fit(t0_ptf)
    bic.append(gmm.bic(t0_ptf))
    aic.append(gmm.aic(t0_ptf))

    if components == 2:
        selected = gmm

x = np.linspace(t0_ptf.min(), t0_ptf.max(), num=1000)[:, None]        
logprob, log_responsibilities = selected._estimate_log_prob_resp(x)
pdf = np.exp(logprob)
responsibilities = np.exp(log_responsibilities)
pdf_individual = responsibilities * pdf[:, np.newaxis]        

'''
plt.figure('flux_lmax')
plt.errorbar(csp_table_lmax['rflam'], csp_table_lmax['Lmax'],yerr=csp_table_lmax['e_Lmax'], fmt='o', label='CSP')
plt.xlabel(r'$\langle f\rangle_{r,15-40}$')
plt.ylabel('$L_{max}$')
plt.legend()
plt.savefig('plots/flux_lmax.pdf')

csp_flam = join(rflam_csp, sbv_csp, join_type='inner', keys='name')
cfa_flam = join(rflam_cfa, sbv_cfa, join_type='inner', keys='name')


plt.figure('sbvflux')
plt.errorbar(cfa_flam['rflam'], cfa_flam['sbv'], yerr=cfa_flam['sbv_err'], fmt='o', label='CfA', color='r')
plt.errorbar(csp_flam['rflam'], csp_flam['sbv'], yerr=csp_flam['sbv_err'], fmt='o', label='CSP')

select = np.logical_and(cfa_flam['rflam'] < 0.25, cfa_flam['sbv'] > 0.6)

print
print cfa_flam['name'][select].tolist()
plt.xlabel(r'$\langle f\rangle_{r,15-40}$')
plt.ylabel('$s_{BV}$')
plt.legend()
plt.savefig('plots/sbvflux.pdf')


#plt.figure('bumpflux')

#plt.errorbar(cfa_flam['rflam'], cfa_table['t_bump'], yerr=cfa_table['t_bump_err'], fmt='o', label='CfA', color='r')
#plt.errorbar(csp_flam['rflam'], csp_flam['sb'], yerr=csp_flam['sbv_err'], fmt='o', label='CSP')
#plt.xlabel(r'$\langle f\rangle_{r,15-40}$')
#plt.ylabel('sbv')
#plt.legend()
plt.show()