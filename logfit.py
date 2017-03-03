import numpy as np
from lmfit import *

# Objective function to minimize when fitting log profile
def resid(params,z,u,eps_u):
    ustr = params['ustr']
    zo = params['zo']
    vk = 0.41
    model = (ustr/vk)*np.log(z/zo)
    return( u-model)/eps_u

def fit_log_profile(u,z,iverbose=1):
    # us lmfit to estimate non-linear least-squares fit
    params = Parameters()
    params.add('ustr',value = 0.01, min = 0., max = 10.)
    params.add('zo'  ,value = 0.001, min = 1.e-6, max = .5)
    eps_u = np.ones_like(u)

    out = minimize(resid,params,method='leastsq', args=(z,u,eps_u))
    ustr_fit = out.params['ustr'].value
    zo_fit = out.params['zo'].value
    ustr_stderr = out.params['ustr'].stderr
    zo_stderr = out.params['zo'].stderr
    if(iverbose):
       print('Best fit u*: {:.4f}; zo: {:.5f}'.format(ustr_fit,zo_fit))
       report_fit(out)
    return ustr_fit,ustr_stderr,zo_fit,zo_stderr

def main():
    # simple test
    ust = 0.2
    zot = 0.001
    zt = np.array((.1,.5,1,2))
    ut = (ust/.41)*np.log(zt/zot)
    us,use,zo,zoe = fit_log_profile(ut,zt,1)
    try:
        np.abs( ust-us )<1.e-5
        print('logfit seems to work')
    except:
        print('seems to be a problem with logfit')

if __name__ == "__main__":
   main()
