
function K = mrpp( Area , E , L , I , gamma1 , gamma2 )

% Calcula la matriz de rigidez de una barra de pórtico plano, como se ve en Cálculo Estructural I.
% El vector de desplazamientos asociado es:
%             U_ij = [  u_i^x  u_i^y  phi_i^z  |  u_j^x  u_j^y  phi_j^z  ]^T


k  = Area*E/L;
k1 = 12*E*I/L^3;
k2 =  6*E*I/L^2;
k3 =  4*E*I/L  ;

A = gamma1^2*k + gamma2^2*k1;
B = gamma1*gamma2*(k-k1);
C = gamma2*k2;
D = gamma2^2*k + gamma1^2*k1;
E = gamma1*k2;

K = [  A  B -C     -A -B -C    ; ...
       B  D  E     -B -D  E    ; ...
      -C  E k3      C -E  k3/2 ; ...
%---------------------------------
      -A -B  C      A  B  C    ; ...
      -B -D -E      B  D -E    ; ...
      -C  E k3/2    C -E k3    ] ;
