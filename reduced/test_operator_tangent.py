from __future__ import annotations

import unittest
import numpy as np

from reduced.causal_graph_circuit import CausalGraphCircuit
from reduced.operator_tangent import (
    compile_dense_passive_graph_tangent,
    run_operator_tangent,
)


def chain(theta):
    # theta scales intrinsic length with diameter held fixed:
    # membrane leak and capacitance scale with length; axial conductance 1/length.
    q=float(theta)
    gl=.02*q
    ax=.12/q
    G=np.array([
        [gl+ax,-ax,0,0],
        [-ax,gl+2*ax,-ax,0],
        [0,-ax,gl+2*ax,-ax],
        [0,0,-ax,gl+ax],
    ],dtype=float)
    C=np.array([.08,.07,.06,.08],dtype=float)*q
    # exact derivatives at q
    dax=-.12/(q*q)
    dG=np.array([
        [.02+dax,-dax,0,0],
        [-dax,.02+2*dax,-dax,0],
        [0,-dax,.02+2*dax,-dax],
        [0,0,-dax,.02+dax],
    ],dtype=float)
    dC=np.array([.08,.07,.06,.08],dtype=float)
    return G,C,dG,dC


def pulse(t,onset,tr=.3,td=5.0):
    x=t-onset
    y=np.where(x>0,np.exp(-x/td)-np.exp(-x/tr),0.0)
    m=float(np.max(y))
    return y/m if m>0 else y


def program(ntime=400,dt=.05):
    t=np.arange(ntime)*dt
    a=np.stack([pulse(t,4),pulse(t,5),pulse(t,6)])
    n=np.stack([pulse(t,4,.29,43),pulse(t,5,.29,43),pulse(t,6,.29,43)])
    return .0015*a,.006*n


class OperatorTangentTests(unittest.TestCase):
    def test_compiler_tangent_matches_finite_difference(self):
        q=1.13
        G,C,dG,dC=chain(q)
        pack=compile_dense_passive_graph_tangent(G,C,dG,dC,[1,2,3],dt_ms=.05)
        eps=1e-6
        Gp,Cp,_,_=chain(q+eps)
        Gm,Cm,_,_=chain(q-eps)
        pp=compile_dense_passive_graph_tangent(Gp,Cp,np.zeros_like(Gp),np.zeros_like(Cp),[1,2,3],dt_ms=.05)
        pm=compile_dense_passive_graph_tangent(Gm,Cm,np.zeros_like(Gm),np.zeros_like(Cm),[1,2,3],dt_ms=.05)
        numP=(pp.passive_step_matrix-pm.passive_step_matrix)/(2*eps)
        numX=(pp.input_step_matrix_mV_per_nA-pm.input_step_matrix_mV_per_nA)/(2*eps)
        self.assertTrue(np.allclose(pack.dpassive_step_dtheta[0],numP,rtol=3e-7,atol=3e-9))
        self.assertTrue(np.allclose(pack.dinput_step_dtheta_mV_per_nA[0],numX,rtol=3e-7,atol=3e-9))

    def test_full_causal_soma_tangent_matches_finite_difference(self):
        q=1.07
        G,C,dG,dC=chain(q)
        pack=compile_dense_passive_graph_tangent(G,C,dG,dC,[1,2,3],dt_ms=.05)
        circuit=CausalGraphCircuit(pack.passive_step_matrix,pack.input_step_matrix_mV_per_nA,[1,2,3],0)
        ga,gn=program()
        tan=run_operator_tangent(circuit,ga,gn,pack.dpassive_step_dtheta,pack.dinput_step_dtheta_mV_per_nA)
        self.assertTrue(tan.all_steps_converged)
        self.assertLessEqual(tan.max_newton_iterations,4)
        self.assertLess(tan.max_site_tangent_consistency_mV_per_unit,2e-9)

        eps=2e-5
        traces=[]
        for qq in (q+eps,q-eps):
            Gq,Cq,_,_=chain(qq)
            pq=compile_dense_passive_graph_tangent(Gq,Cq,np.zeros_like(Gq),np.zeros_like(Cq),[1,2,3],dt_ms=.05)
            cq=CausalGraphCircuit(pq.passive_step_matrix,pq.input_step_matrix_mV_per_nA,[1,2,3],0)
            traces.append(cq.run(ga,gn).soma_depolarization_mV)
        numeric=(traces[0]-traces[1])/(2*eps)
        self.assertTrue(np.allclose(tan.soma_tangent_mV_per_unit[0],numeric,rtol=2e-5,atol=2e-6))

    def test_two_parameter_stack(self):
        q=1.0
        G,C,dG,dC=chain(q)
        # second parameter scales all membrane capacitances only
        dGs=np.stack([dG,np.zeros_like(G)])
        dCs=np.stack([dC,C])
        pack=compile_dense_passive_graph_tangent(G,C,dGs,dCs,[1,2,3],dt_ms=.05)
        ga,gn=program(ntime=120)
        circuit=CausalGraphCircuit(pack.passive_step_matrix,pack.input_step_matrix_mV_per_nA,[1,2,3],0)
        tan=run_operator_tangent(circuit,ga,gn,pack.dpassive_step_dtheta,pack.dinput_step_dtheta_mV_per_nA)
        self.assertEqual(tan.soma_tangent_mV_per_unit.shape,(2,120))
        self.assertTrue(np.all(np.isfinite(tan.soma_tangent_mV_per_unit)))


if __name__ == '__main__':
    unittest.main()
