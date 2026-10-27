import math, numpy as np

class QuantumCircuit:
 def __init__(self,n,name="THEOS"):
  self.n=n
  self.name=name
  self.ops=[]
  self.t8=80
  self.lam=1.818e-10
  self.f=5.128
  self.state=np.zeros(2**n,dtype=complex)
  self.state[0]=1

 def _gate(self,g,q):
  mats=[]
  for i in range(self.n):
   mats.append(g if i==q else np.eye(2))
  f=mats[-1]
  for m in reversed(mats[:-1]):
   f=np.kron(f,m)
  return f

 def h(self,q):
  H=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
  self.state=self._gate(H,q)@self.state
  self.ops.append(f"H q{q}")

 def ry(self,th,q):
  c=math.cos(th/2)
  s=math.sin(th/2)
  R=np.array([[c,-s],[s,c]],complex)
  self.state=self._gate(R,q)@self.state
  self.ops.append(f"RY {th:.2f} q{q}")

 def cx(self,c,t):
  new=np.zeros_like(self.state)
  for i in range(len(self.state)):
   bit=(i>>c)&1
   if bit:
    j=i^(1<<t)
    new[j]=self.state[i]
   else:
    new[i]=self.state[i]
  self.state=new
  self.ops.append(f"CX {c}->{t}")

 def show(self):
  parts=[]
  for i,a in enumerate(self.state):
   if abs(a)>1e-6:
    b=format(i,'b').zfill(self.n)[::-1]
    parts.append(f"{a.real:.3f}|{b}>")
  return " + ".join(parts)

 def entropy(self, base='e'):
  p=np.abs(self.state)**2
  p=p[p>1e-12]
  if base=='2':
   return -np.sum(p*np.log2(p))
  return -np.sum(p*np.log(p))

 def von_neumann(self, q=None):
  if q is None:
   rho=np.outer(self.state, np.conj(self.state))
   eig=np.linalg.eigvalsh(rho)
  else:
   psi=self.state.reshape(2,2)
   if q==0:
    rho=psi @ psi.conj().T
   else:
    rho=psi.T.conj() @ psi.T
   eig=np.linalg.eigvalsh(rho)
  eig=eig[eig>1e-12]
  return -np.sum(eig*np.log(eig))

 def draw(self):
  print(f"THEOS {self.name}")
  print(f"T8 {self.t8}x {self.f}GHz")
  for o in self.ops:
   print(o)
  print("State:",self.show())
  print(f"Gibbs S/kB={self.entropy():.4f} Shannon={self.entropy(base='2'):.4f} bits")

def psi():
 q=QuantumCircuit(1,"0.6 0.8")
 q.ry(2*math.acos(0.6),0)
 return q

def bell():
 q=QuantumCircuit(2,"Bell")
 q.h(0)
 q.cx(0,1)
 return q

if __name__=="__main__":
 psi().draw()
 print()
 bell().draw()
 print(f"Reduced S={bell().von_neumann(q=0):.4f} = ln2")
