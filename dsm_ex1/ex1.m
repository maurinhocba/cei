
close all
clear variables
clc

% conversion factors
feet2mm = 304.8;
ksi2MPa = 6.89476;
in22mm2 = 645.16;
in42mm4 = 416231.426;
kip2N   = 4448.22;

% nodes
n1=[ 0 20]*feet2mm;
n2=[20 20]*feet2mm;
n3=[20  0]*feet2mm;

% elements
E=29e3*ksi2MPa;
A=10*in22mm2;
I=500*in42mm4;

aux=n2-n1;
L1=norm(aux);
t1=aux/L1;
K1g = mrpp( A, E, L1, I, t1(1), t1(2) );

aux=n3-n2;
L2=norm(aux);
t2=aux/L2;
K2g = mrpp( A, E, L2, I, t2(1), t2(2) );

% assembled global matrix
K_full=zeros(9,9);
K_full(1:6,1:6)=K_full(1:6,1:6)+K1g;
K_full(4:9,4:9)=K_full(4:9,4:9)+K2g;

% assembled global load vector
L_full=zeros(9,1);
L_full(4)=5*kip2N;

% BCs
K_wBCs=K_full;
K_wBCs(:,[2 7:9])=[]; % delete rows
K_wBCs([2 7:9],:)=[]; % delete cols

L_wBCs=L_full;
L_wBCs([2 7:9])=[]; % delete rows

% disp.
U_wBCs = K_wBCs \ L_wBCs;

U_full=zeros(9,1);
U_full(1)=U_wBCs(1);
U_full(3:6)=U_wBCs(2:end);

% elemental load vectors
U1g=U_full(1:6);
L1g=K1g*U1g;
RT1=mRotTPP( t1(1), t1(2));
L1l=RT1*L1g;

U2g=U_full(4:9);
L2g=K2g*U2g;
RT2=mRotTPP( t2(1), t2(2));
L2l=RT2*L2g;
