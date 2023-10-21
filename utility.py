
from hashlib import sha256
from math import log2, ceil, exp, pi, sqrt
from random import randint, normalvariate
from numpy import polymul,polyadd
from numpy.random import normal, binomial
from numpy import negative

n = 1024
mu= ceil(pow(log2(n), 2))
sigma=pow(n, 0.75)
q = 1187
coef_c_bound = ceil(log2(n))
deg_c_bound = n/2
coef_y_bound = ceil(pow(n, 1.5)*sigma*pow(log2(n), 3))
coef_z_bound = ceil((n-1)*pow(n, 0.5)*sigma*pow(log2(n), 3))

#SampU(a, b, N) generates N samples from {a, a+1, ..., b} and returns an array v
def SampU(a, b, N):
      v=[]
      for i in range(N):
          v.append(randint(a, b))
      return v

# Norm(v) computes the max_i |v_i|_q (where |v_i|_q=min(|v_i|, q-|v_i|)
def Norm(v):
      if len(v)==0:
         raise Exception("vector should have at least one element")

      max=0
      for i in range(len(v)):
          temp0=v[i]%q
          temp=min(temp0, q-temp0) 
          if max<temp:
             max=temp
      return max

# PolyMod(c) computes the representative of c in R_q=Z[x]/(q, x^n+1) with deg <n and coeff in Zq. 
def PolyMod(c):
      if len(c)<=0:
         raise Exception("invalid polynomial input")
      
      o=c
      while len(o)>n:
           A=o[-n:]
           o=polyadd(negative(o[:-n]), A)%q
      return o
      
#MULT(a, b) computes the mulitiplication of polynomial a and b in R_q=Z[x]/(q, x^n+1)
def MULT(a, b):
      if len(a)==0 or len(b)==0:
         raise Exception("invalid input for polynomial multiplication")
   
      c= polymul(a, b)%q
      #print("c[0]={}".format(c[0]))
      return PolyMod(c)

#MULT1tom(a, b) multiplies a polynomial to a vector of polynomials via the entry-wise multiplication 
def MULT1tom(a, b):
      temp=[]
      for i in range(len(b)):
            temp.append(MULT(a, b[i]))
      return temp


#AGG(v, P) computes the aggregation of polynomial vectors v and P: sum_i v[i]*P[i] over R_q. 
def AGG(v,P):
      if len(v)!=len(P):
         raise Exception("coeffient vector and polynomial vector do not have a matching length")
      if len(v)==0:
         raise Exception("aggregation requires at least one input element")
    
      sum=[0]
      for i in range(len(v)):
          sum=polyadd(sum, MULT(v[i], P[i]))%q
      return sum

#Hash(b) computes the hash of byte-stream b and output a n/2-sized array v with v_i in [-log2(n), log2(n)]
#the output is in the challenge set {\cal C} in the paper; 
def Hash(b):
      L=ceil(n/64)
      E=ceil(log2(n))
      A=b''
      for i in range(L):
          A=A+sha256(b+i.to_bytes(2, 'little')).digest() #create n/2 byte stream
          
      v=[]
      for j in range(ceil(n/2)):
          v.append(A[j]%(2*E)-E)  # Here we assume E<=128 (i.e. n<=2^{128}).

      return v
 

#SampD_R() implements D_{R, sigma}, using the simple rejection sampling for D_{Z,sigma} in GPV08
#the output is an array of size n from range [-sigma*6, sigma*6] 
#(6 here makes sure each sample before rejection decision lies in the interval except for probability ~ 10^(-50).  
def SampD_R():
       high=ceil(6*sigma) 
       low=-high
       
       v=[] 
       d=pi/pow(sigma, 2) 
       for i in range(n):
           s=randint(low, high)
           prob=exp(-pow(s, 2)*d)   #This is rho_sigma(s) for the rejection decision  

           while (binomial(n=1, p=prob)==0):
                 s=randint(low, high)
                 prob=exp(-pow(s, 2)*d)


           v.append(s)

       return v


#vecpolyadd adds the vector of polynomials a and b via the entry-wise addition. 
def vecpolyadd(a, b):
      temp=[]
      if len(a)!=len(b):
            raise Expection("inputs do not have the same length")
 
      for i in range(len(b)):
            temp.append(polyadd(a[i], b[i])%q)
      return temp

