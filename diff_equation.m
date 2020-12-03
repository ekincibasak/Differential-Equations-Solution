
%//This matlab function calculates the flavor evolution of neutrinos in the Supernova.
%// Input:
%// Energy vector E. It should be given in MeV
%// Output:
%// r: vector of radius values. 
%// P: Matrix of survival probabilities.

clc
clear all
close all
%//Mlxlng angles and mass differance  from partlcle data group
%theta12=asin(sqrt(0.304))
theta13=0.1472;
delta_m_sqr31= 6.6845

%Stanndart parametrlzatlon of mlxlng matrlx
U=[cos(theta13) sin(theta13) ; -sin(theta13) cos(theta13)]
%uniter matrix for writing the different base.
UMUdagger=U*[0 0 ; 0 -6.6845]*ctranspose(U)


R=3.5
psi=[1;0];
E=10;  % in Mev



Hvacuum=(1/E)*UMUdagger
h=0.1;
l=1;
r(l)=0;
P(l)=abs(psi(1))^2;
%solution of the linear differential equation by fourth order runge-kutta
while r(l)<R
H=Hvacuum+0.189*expm(-10.54*r(l)/R)*[1 0 ; 0 0 ]; %in km^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1= h* (H*psi*(1/i));
k2=h*((Hvacuum+0.189*expm(-10.54*(r(l)+h/2)/R)*[1 0 ; 0 0 ])*(psi+k1/2)*(1/i));
k3=h*((Hvacuum+0.189*expm(-10.54*(r(l)+h/2)/R)*[1 0 ; 0 0 ])*(psi+k2/2)*(1/i));
k4=h*((Hvacuum+0.189*expm(-10.54*(r(l)+h)/R)*[1 0 ; 0 0 ])*(psi+k3)*(1/i));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi = psi +(k1+2*k2+2*k3+k4)/6;
l=l+1; 
 r(l)=r(l-1)+h;
 P(l)=abs(psi(1))^2;
end 

plot(r,P)
xlabel('Radius(km)')
ylabel('transition probability')
