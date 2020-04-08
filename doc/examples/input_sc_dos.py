import numpy
from pytriqs.Base.Lattice.TightBinding import *
from pytriqs.Base.DOS.Hilbert_Transform import *
from pytriqs.Base.GF_Local import *

#
# SC, 1 site DMFT
#

N_Loops                    =  1
Beta                       =  5.0
N_Cycles                   =  50000
U_interact                 =  5.0
t                          = -1.0
Field_SC                   =  0.10
Field_SC_loc               =  0.00
mu                         =  U_interact/2.0

#
#  Solver
#
from pytriqs.Solvers.InteractionExpansion.Solver import Solver
class Solver_no_fit(Solver):
  def fitTails(self): pass

S = Solver_no_fit(Beta = Beta,
                  GFstruct = [('1',['p','m'])],
                  N_Cycles = N_Cycles,
                  Length_Cycle = 50,
                  N_Matsubara_Frequencies = 300)

# U[a,b][i,j]
S.Umatrix = [[ numpy.array([[0,U_interact],[U_interact,0]]) ]]

# alpha[a][i,s]
S.alpha = [ numpy.array([[0.5001, 0.6],[0.42, 0.4]]) ]

#  Init Self-Energy
S.Sigma << numpy.array([[mu,Field_SC],[Field_SC,mu]],numpy.complex)

#   Lattice
hop = { ( 1,0) : [[ t]],
        (-1,0) : [[ t]],
        (0, 1) : [[ t]],
        (0,-1) : [[ t]] }

BL = bravais_lattice(Units = [(1,0,0) , (0,1,0) ], Orbital_Positions= {"" :  (0,0,0)} )
TB = tight_binding (BL, hop)
d = dos (TB, nkpts= 500, neps = 101, Name = 'dos')[0]
HT = Hilbert_Transform(d)

# Green functions
G = GF( NameList = ['SC'],
        BlockList = [GFBloc_ImFreq(Indices = ['p','m'], Mesh = S.G.mesh)],
        Copy = False )
Sigma = G.copy()

#   My DMFT loop ...
for it in range(N_Loops):

     # Embedding
     Sigma['SC']['p','p']  = S.Sigma['1']['p','p']
     Sigma['SC']['p','m']  = S.Sigma['1']['p','m']
     Sigma['SC']['m','p']  = S.Sigma['1']['m','p']
     Sigma['SC']['m','m']  = S.Sigma['1']['m','m']

     def Epsilon_Hat_SC(eps) :
       a = numpy.zeros( (len(eps) if isSequenceType(eps) else 1,2,2) )
       a[:,0,0] = eps; a[:,1,1] = - numpy.array(eps)
       return a

     Dc = HT (Sigma['SC'], Field = numpy.array([[-mu ,Field_SC_loc],[Field_SC_loc, -mu]]), Epsilon_Hat= Epsilon_Hat_SC, Res = G['SC'] ).density()
     print("Total density = %f"%((Dc[0,0]-Dc[1,1]+1.0).real))

     # Extraction
     S.G['1']['p','p'] = G['SC']['p','p']
     S.G['1']['p','m'] = G['SC']['p','m']
     S.G['1']['m','p'] = G['SC']['m','p']
     S.G['1']['m','m'] = G['SC']['m','m']

     # Finally get S.G0
     S.G0 << inverse(S.Sigma + inverse(S.G))

     # Run and save
     S.G0.save("G0")
     S.G.save("G")
     S.Solve()

     g_tau = GFBloc_ImTime(Indices = ['p','m'], Beta = Beta, NTimeSlices = 10000)
     g_tau << InverseFourier(S.G['1'])
     g_tau.save("G_tau")

     Dc = S.G['1'].density()
     M = Dc[0,0]-Dc[1,1]+1.0
     MPI.report("Density = %.5f"%M.real)
