
import numpy as np

def floc_props(Df, Dp0, nf, rhos = 2650., rhow = 1025.0, nu = 0.001 ):
   g = 9.81
   rhof =rhow+(rhos-rhow)*(Dp0 /Df)**(3.0-nf)   # floc density
   wsf =g*(rhof-rhow)*Df**2.0 /(18.0*nu);       # settling velocity (Stokes)
   return rhof, wsf

def main():
   # test functions
   print("Testing floc props:\n")
   Dp0 = 4.0e-6
   Dpmax = 1500.e-6
   Df = 10.**np.linspace( np.log10(Dp0),np.log10(Dpmax),5)
   #Df = 200.e-6
   nf = 2.
   rhof, wsf = floc_props( Df, Dp0, nf )
   print("rhof: ",rhof)
   print("wsf : ",wsf)

if __name__ == "__main__":
   main()


