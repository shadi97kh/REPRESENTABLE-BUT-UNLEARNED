"""Exact SMALL finite-state reference; standard bounded linear programming.

Not a new GNN layer, differentiable solver, end-to-end theorem or biological CI.
All arithmetic is rational. Inputs are integers or explicit decimal/rational strings.
Aborts without an answer when the enumeration budget is exceeded.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb

def rational(x):
 if isinstance(x,float):raise ValueError('Use explicit rational/decimal strings, not implicit floating-point inputs')
 return F(x)

def solve_square(A,b):
 n=len(b);a=[list(row)+[v] for row,v in zip(A,b)]
 for j in range(n):
  p=next((i for i in range(j,n) if a[i][j]),None)
  if p is None:return None
  a[j],a[p]=a[p],a[j];v=a[j][j];a[j]=[x/v for x in a[j]]
  for i in range(n):
   if i!=j:
    v=a[i][j];a[i]=[x-v*y for x,y in zip(a[i],a[j])]
 return tuple(row[-1] for row in a)

def dot(a,b):return sum(x*y for x,y in zip(a,b))

def exact_range(lower,upper,observations,query,classes=None,eta='0',budget=200000):
 """mu has one coordinate per UNIQUE measured state, including shared endpoints.

 observations=[(coefficient row, lower bound, upper bound)]. Encoding classes
 optionally impose |mu_i-mu_j|<=2 eta. These are EXTERNAL model assumptions.
 bounds and observations refer to population response means only if their
 input provenance supplies that interpretation; no sampling model is assumed.
 """
 n=len(lower);lo=list(map(rational,lower));hi=list(map(rational,upper));q=list(map(rational,query));e=rational(eta)
 if n==0 or len(hi)!=n or len(q)!=n or e<0:raise ValueError('invalid dimensions or approximation bound')
 if any(l>u for l,u in zip(lo,hi)):return {'status':'empty','reason':'contradictory coordinate bounds'}
 A=[];b=[]
 def add(row,rhs):A.append(tuple(map(rational,row)));b.append(rational(rhs))
 for i in range(n):
  row=[F(0)]*n;row[i]=F(1);add(row,hi[i]);add([-v for v in row],-lo[i])
 for row,l,u in observations:
  if len(row)!=n:raise ValueError('observation dimension')
  row=list(map(rational,row));l=rational(l);u=rational(u)
  add(row,u);add([-v for v in row],-l)
 if classes is not None:
  if len(classes)!=n:raise ValueError('class dimension')
  for i,j in combinations(range(n),2):
   if classes[i]==classes[j]:
    row=[F(0)]*n;row[i]=1;row[j]=-1;add(row,2*e);add([-v for v in row],2*e)
 # Exact deduplication avoids needless repeated active sets; retain strictest rhs.
 dedup={}
 for row,v in zip(A,b):dedup[row]=min(v,dedup.get(row,v))
 A=list(dedup);b=[dedup[row] for row in A];count=comb(len(A),n)
 if count>budget:return {'status':'unresolved','reason':'exact enumeration budget','active_sets':count,'budget':budget}
 vertices=set();max_violation=F(0)
 for ix in combinations(range(len(A)),n):
  x=solve_square([A[i] for i in ix],[b[i] for i in ix])
  if x is not None and all(dot(row,x)<=rhs for row,rhs in zip(A,b)):vertices.add(x)
 if not vertices:return {'status':'empty','active_sets':count,'reason':'bounded polytope has no vertex'}
 low=min(vertices,key=lambda x:(dot(q,x),x));high=max(vertices,key=lambda x:(dot(q,x),x))
 return {'status':'exact','lower':str(dot(q,low)),'upper':str(dot(q,high)),
         'lower_witness':list(map(str,low)),'upper_witness':list(map(str,high)),
         'vertices':len(vertices),'active_sets':count,'primal_violation':'0',
         'global_objective_gap':'0','dual_residual':None,
         'certificate':'exhaustive exact vertex enumeration; no dual solver',
         'arithmetic':'rational; no rounding allowance required for supplied rational problem'}

def union_range(results):
 """A finite union retains all admitted representation alternatives."""
 if any(r['status']=='unresolved' for r in results):return {'status':'unresolved','reason':'at least one representation unresolved'}
 valid=[r for r in results if r['status']=='exact']
 if not valid:return {'status':'empty'}
 return {'status':'exact','lower':str(min(F(r['lower']) for r in valid)),
         'upper':str(max(F(r['upper']) for r in valid)),
         'interpretation':'outer interval hull of union; gaps may exist'}

if __name__=='__main__':
 import json,sys
 print(json.dumps(exact_range(**json.load(open(sys.argv[1]))),indent=2))
