%% Pailler_versuch3
clc;
clear;
close all;


%% Schlüsselerzeugung (Beispiel)
tic
% Since prime numbers are too big to find in Matlab, the found prime
% numbers from Python are used
% primes = readtable('two primenumbers from python.csv');
% p= table2array(primes(1,:));
% q= table2array(primes(2,:));
% since there is a conversion of number, it doesn't work correctly.
% n=5000;
% primes(n);
% primes= primes(n);
% 
% pq=maxk(primes,2);
% p= pq(1);
% q= pq(2);
p= 973049797;
q= 553957561;
disp('the chosen prime numbers from python are : ')
disp('p = : ')
disp(p)
disp('q = : ')
disp(q)
% öffentliche Schlüssel
n= sym(p)*sym(q);% n ist öffentliche Schlüssel (cipher)
g= n+1;%Pick a random integer g in the set Z*n2 (integers between 1 and n2).
% g= randi(n2); % wenn man g als eine randomisierte Zahl auswählen möchte
% geheime Schlüssel
r1=152;
p_1=p-1;
q_1=q-1;
lambda= lcm(sym((p_1)),sym((q_1))); %compute ambda as lcm(p-1,q-1)
cipher=lambda;
[~,mu]=gcd(Lx(powermod(g,lambda,n^2),n),n);
numbers = readtable('the list of numbers.csv');
m1= table2array(numbers(1,:));
disp('Taken time in key generation:')
toc

disp("g="), disp(g)
disp("r="), disp(r1)
disp("================")
disp("Mu:"), disp(mu)
disp("Lambda:"),disp(lambda)
disp("================")
disp("Public key (n,g):"), disp(n), disp(g)
disp("Private key (lambda,mu):"), disp(lambda), disp(mu)
disp("================")

%% Encryption
tic
for i = 1:99
    k1(i)= powermod(sym(g), m1(i), sym(n*n));
end

k2 = powermod(r1, n, sym(n*n)); % powermod is different from pow in python?

c1 = mod((k1*k2),sym(n*n)); % (k1*k2) is diffrent
disp('to show that the encryption functioned well')
disp('first element of Cipher as an example:');
disp(c1(1));
% Those c would be the encrypted numbers

disp('Taken time in encryption:')
toc
%% Decryption
tic
for i=1:99
   x(i)= powermod(c1(i), lambda, sym(n*n));
   l(i)=Lx(x(i),n); % here is diffrent
   mess(i)=mod((l(i)*mu),sym(n));
end

disp('decrypted Text (plain text):');
disp(mess);
disp('Taken time in decryption:')
toc
%% Homomorphic addition
% Multiplikation mit einem Plaintext
%More generally, a ciphertext raised to a constant k will decrypt to the product of the plaintext and the constant,
% E(m.1,r.1) = c1
disp('Now we start to homomorphic operations: k1*x1+k2*x2+k3*x3..+k100*x100')
tic
% choose the k
numbers = readtable('the list of numbers.csv');
k= table2array(numbers(2,:));

for i=1:99
     c_mul(i)= powermod(sym(c1(i)), k(i), sym(n*n));
     x_mul(i)= powermod(sym(c_mul(i)), lambda, sym(n*n));
     l_mul(i)=Lx(x_mul(i),sym(n));
     m_mul(i)=mod((l_mul(i)*mu),sym(n)); % plain text after multiplication
end

disp('to show that the encryption functioned well')
disp('first ciphered element of homomorphic multiplication:')
disp(c_mul(1))
disp("Message:") ,disp(mess)
%% Addition of two ciphertexts
% %D.priv(E.pub(m1)⋅E.pub(m2) mod n^2)=m.1+m.2 mod n
%left side of equation
c_add_prod=mod(sym(prod(c_mul,'all')), sym(n*n));
c_add= mod(sym(c_add_prod),sym(n*n));
x_add= powermod(sym(c_add), lambda, sym(n*n));
l=Lx(x_add,n);
m_add= mod((l * mu) , n);

%Addition of 100 ciphertexts
disp('The decrypted message after the homomorphic addition:')
disp(m_add)

disp('Taken time in Homomorphic Operation:')
toc

%% The result from the plain text
disp('This part will be the calculation without encryption (to compare the result)')

for i = 1:99
 summe_list(i) = k(i)*m1(i);
end
summe_from_original=mod(sum(summe_list),n);

disp('The multiplicated Message e.g.(k(1)):')
disp(k(1))
disp('The decrypted message after the homomorphic multiplication (k*m1)')
disp(m_mul)
disp('Formel: Summe = k1*x1+k2*x2+k3*x3...+k100*x100')
disp('The target result without encryption:')
disp(summe_from_original)
%%

function Lx = L(x,n)
    Lx = (x-1)/n;  
end