#!/usr/bin/python3

from utility import *
import pickle
import time

M=b"this is the message to be signed"


#key generation: generate system parameter a and key pairs ((s1[i], s2[i]), u[i]), i=0, 1, 2 
print("key generation..")
a=SampU(0, q-1, n)
s1=[]
s2=[]
u=[]

for i in range(3): #generate key for signer i=0, 1, 2
    s1.append(SampD_R())
    s2.append(SampD_R())
    u.append(polyadd(MULT(a, s1[i]), s2[i])%q)

y1=[]
y2=[]
v=[]
r=[]

print("Round R-1 starts...")
start_time=time.time()  # the protocol start_time
#Round R-1:  signer i  generates v[i]=ay1[i]+y2[i]. 
for i in range(3):
    y1.append([])
    y2.append([])
    v.append([])    
    for j in range(mu):
        y1[i].append([])
        y2[i].append([])
        v[i].append([])
        y1[i][j]=SampU(-coef_y_bound, coef_y_bound, n)
        y2[i][j]=SampU(-coef_y_bound, coef_y_bound, n)
        v[i][j]=polyadd(MULT(a, y1[i][j]), y2[i][j])%q

    r.append(Hash(b'0'+pickle.dumps(v[i])+pickle.dumps(u[i])))

print("Round R-1 done and Round R-3 starts...") 

#Roud R-2:  #receiving all r[i], user needs to send v[i] to others; no computation needed.

#Round R-3:  #first, check if r[i]=Hash(0|v[i]|u[i])
for i in range(3):
    if r[i]!=Hash(b'0'+pickle.dumps(v[i])+pickle.dumps(u[i])):
        raise Exception("verification of r fails")
lambd=[]
for i in range(3):
    lambd.append(Hash(b'0'+pickle.dumps(u[i])+pickle.dumps(u)))

bar_u=MULT(lambd[0], u[0])
bar_v=MULT1tom(lambd[0],v[0])
for i in range(1, 3):
    bar_u=polyadd(bar_u,MULT(lambd[i], u[i]))%q
    bar_v=vecpolyadd(bar_v,MULT1tom(lambd[i], v[i]))

c=Hash(b'1'+pickle.dumps(bar_u)+pickle.dumps(bar_v)+M)

z1=MULT1tom(c, s1)
z2=MULT1tom(c, s2)
for i in range(3):
    for j in range(mu):
        z1[i]=polyadd(z1[i], y1[i][j])%q
        z2[i]=polyadd(z2[i], y2[i][j])%q

print("R-3 done and output generation...")

#output aggregation: aggregate z1[i], z2[i], i=0, 1, 2 and output bar_z1, bar_z2, bar_v
#aggregated public-key: bar_u|t (here t=3)
bar_z1=MULT(lambd[0], z1[0])
bar_z2=MULT(lambd[0], z2[0])
for i in range(1, 3):
    bar_z1=polyadd(bar_z1,MULT(lambd[i], z1[i]))%q
    bar_z2=polyadd(bar_z2,MULT(lambd[i], z2[i]))%q

#Verify: 
#check if ||bar_z1||_infty<eta_t and ||bar_z2||_infty<eta_t
#check if bar_v[1]+...bar_v[mu]=a*bar_z1+bar_z2-bar_u*c

print("output done and verification begin..")
if Norm(bar_z1)>eta_t or Norm(bar_z2)>eta_t: 
    raise Exception("signature verification fails: bar_z1 or bar_z2 is not short")

Left=[0]
for i in range(mu):
    Left=polyadd(Left, bar_v[i])%q

Right=MULT(a, bar_z1)
Right=polyadd(Right, bar_z2)%q
Right=polyadd(Right, negative(MULT(bar_u, c)))%q

LT=Left.tolist()
RT=Right.tolist()

if LT!=RT:
     flag=0
     print("The Multi-signature verification fails")
else:
     print("The Multi-signature passes the verification!")

end_time=time.time()
runtime=end_time-start_time
print("the total protocol execution time is {} seconds".format(runtime))

