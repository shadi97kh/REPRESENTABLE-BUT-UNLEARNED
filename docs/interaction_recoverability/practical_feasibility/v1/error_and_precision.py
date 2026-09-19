"""Exact rational certificates and public precision formulas; no sampled data."""
from fractions import Fraction as F
import json
import math
from pathlib import Path


def exp_interval_positive(x,n=100):
    term=total=F(1)
    for k in range(1,n+1):
        term*=x/k;total+=term
    next_term=term*x/(n+1)
    return total,total+next_term/(1-x/(n+2))


def upper_float(q):
    return math.nextafter(float(q),math.inf)


def main():
    p=Path(__file__).with_name('error_and_precision.json')
    if p.exists():raise FileExistsError('Preserve previous certificates')
    lam=F(271,290);mu=F(129,152)
    checks={
      'sqrt829_upper':F(288,10)**2>829,
      'sqrt269_upper':F(1641,100)**2>269,
      'lower_root_gap_below_five':27*F(288,10)+13*F(1641,100)<1000,
      'sqrt884_lower':F(2973,100)**2<884,
      'sqrt244_lower':F(1562,100)**2<244,
      'upper_root_gap_above_five':28*F(2973,100)+12*F(1562,100)>1000,
      'entry_to_local_interval_after_32':lam**32*F(11,20)<F(1,10),
      'hp14_below_129':F(492,100)**2<F(29,10)**2*F(296,100),
      'hp26_above_152':F(1452,100)**2>F(52,10)**2*F(776,100),
      'local_denominator_above_007':mu<F(93,100)**2,
      'sqrt12_above_346':F(346,100)**2<12,
      'positive_warp_shift_below_051':F(14,10)/(1+F(1,5)*(10-F(14,11)))<F(51,100),
    }
    for x,bound in ((F(4),55),(F(5),149),(F(1,10),F(1106,1000)),(F(841,200),70)):
        checks['exp_'+str(x)+'_upper']=exp_interval_positive(x)[1]<bound
    assert all(checks.values())
    factor=mu**240
    a=F(209,100)
    etop=1/exp_interval_positive(F(348,100))[0]
    positive_tail=F(75)/F(346,100)*etop*(18/a+1/a**2)
    eneg=1/exp_interval_positive(F(48),200)[0]
    negative_tail=F(66)/F(346,100)*eneg*(F(17,7)+F(1,49))
    records=[]
    for N in (1000,10000,100000,1000000):
        K=2
        while lam**(K-2)>F(1,N*N):K+=1
        R=math.sqrt(1.5*math.log(N));M=math.ceil(4*R+30)+2
        H1=lambda x:1.1+.2*x
        L0=4+4*H1(R+M+4)
        JY=2*R*(4*H1(R)+2*M*L0)+36*K*(2*M+1)*L0
        cd=10350*(K+1)
        tau=1/(16*N*(1+cd+JY))
        records.append({'N':N,'K':K,'J_Y':JY,'C_data':cd,'tau_input':tau,
          'sufficient_fraction_bits':math.ceil(-math.log2(tau)),
          'cap_integer_magnitude_bits':math.ceil(math.log2(N*N+1)),
          'scope':'Public sufficient input allocation, not a full working-bit certificate.'})
    result={'status':'passed','rational_certificates':checks,
      'central_complete_influence_L2_remainder_K512_upper':upper_float(5000000000000*factor),
      'central_four_coefficient_vector_L2_remainder_upper':upper_float(200000000000*factor),
      'response_tail_L2_remainder_R12_upper':upper_float(positive_tail+negative_tail),
      'lattice_L2_remainder_M24_R12_upper':1e-45,
      'scope':'Analytic bounds for exact truncated influence, not a certificate of the floating-point quadrature.',
      'precision_allocations':records}
    p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
