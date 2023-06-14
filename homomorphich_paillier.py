from random import randint
import libnum
import sys
import random
import numpy as np ## in order to save the time

def gcd(a,b):
    """Compute the greatest common divisor of a and b"""
    while b > 0:
        a, b = b, a % b
    return a
    
def lcm(a, b):
    """Compute the lowest common multiple of a and b"""
    return a * b // gcd(a, b)

def L(x,n):
	return ((x-1)//n)


#%% Key generation (Schlüsselerzeugung)
find_primenumber=100

# Es funktioniert nicht, wenn z.b. 2**513 ist.
# sys.maxsize (36bits): 9223372036854775807
# Algorism of eratosthenes (Primzahlen finden )
a = [False,False] + [True]*(find_primenumber-1)
primes=[]

for i in range(2,find_primenumber+1):
  if a[i]:
    primes.append(i)
    for j in range(2*i, find_primenumber+1, i):
        a[j] = False

# print('the choosable prime numbers are : ',primes)
# p= random.sample(primes, 2)[0]
# q=random.sample(primes, 2)[1] ## two prime numbers are randomly choosed
# choose the two biggest prime numbers
# it must be two biggest prime numbers in oder to calculate stable results
# according to the p,q, the decrypted messages can be different.
p= primes[-1]
q= primes[-2]

# öffentliche Schlüssel public key (n,g)
n = p*q 
# um Schlüssellänge zu vergrößen, muss n noch größer werden
# sys.getsizeof(n)
# in case n = 2499995600000711 size of n is 32
g = n+1 #Pick a random integer g in the set Z∗n2 (integers between 1 and n2)
gLambda = lcm(p-1,q-1) ## in Matlab cipher = lambda

m1= int(input('Enter plaintext (123) :')) #must be string type

if (len(sys.argv)>1):
	m=int(sys.argv[1])

if (len(sys.argv)>2):
	p=int(sys.argv[2])

if (len(sys.argv)>3):
	q=int(sys.argv[3])

if (p==q):
	print("P and Q cannot be the same")
	sys.exit()
l = (pow(g, gLambda, n*n)-1)//n
gMu = libnum.invmod(l, n)    


if (gcd(g,n*n)==1):
	print("g is relatively prime to n*n")
else:
	print("WARNING: g is NOT relatively prime to n*n. Will not work!!!")
# geheime Schlüssel (private key)
r1=152


#%% Encryption
# changed: pow(a,b,n) for modular exponentiation a^b mod n
# --> much faster than the standard pow()
# shortened the code
k1 = pow(g, m1, n**2)
k2 = pow(r1, n, n**2)
c1=((k1*k2)%(n**2));


#%% Decryption
# changes: shortened the code, used pow(.,.,.)
x = pow(c1, gLambda, n**2)
l=L(x,n)
mess=((l * gMu) % n) ## this must be same with m_list

    
print("p=",p,"\tq=",q)
print("g=",g,"\tr=",r1)
print("================")
print("Mu:\t\t",gMu,"\tgLambda:\t",gLambda)
print("================")
print("Public key (n,g):\t\t",n,g)
print("Private key (lambda,mu):\t",gLambda,gMu)
print("================")
print("Message:\t",mess)

print("Cipher:\t\t",c1)
print("Decrypted Message:\t",mess)

print("================")
print("Now we will add a ciphered value of 2 to the encrypted value")

#%%Example, Homomorphic addition

m2= 37
r2= 999
k3 = pow(g, m2, n**2)
k4 = pow(r2, n, n**2)
c2=((k3*k4)%(n**2))
#Addition of two ciphertexts
#D.priv(E.pub(m1)⋅E.pub(m2) mod n^2)=m.1+m.2 mod n
c_add=((c1*c2)%(n**2))
x_add= pow(c_add, gLambda, n**2)
l=L(x_add,n)
m_add=((l * gMu) % n)
print("the added Message:\t",m2)
print("The decrypted message after the homomorphic addition", m_add)

#%% Example, Homomorphic multiplication by 1
# Multiplikation mit einem Plaintext
k=25
#More generally, a ciphertext raised to a constant k will decrypt to the product of the plaintext and the constant,
# E(m.1,r.1) = c1
c_mul=pow(c1, k, n**2)
x_mul= pow(c_mul, gLambda, n**2)
l = L(x_mul,n)
m_mul=((l * gMu) % n)


print("the multiplicated Message (k):\t",k)
print("The decrypted message after the homomorphic multiplication (m1*m3)", m_mul)
