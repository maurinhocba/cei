
function RT = mRotTPP( g1, g2)

% Calcula la matriz de rotaci�n de un vector de cargas o de desplazameinto
% de p�rtico plano, para ir de coords. globales a locales (t�picamente
% lo que se usa para resolver problemas en parciales) de a cuerdo a como
% se ve en C�lculo Estructural I.
% 
% El vector de desplazamientos asociado es:
%             U_ij = [  u_i^x  u_i^y  phi_i^z  |  u_j^x  u_j^y  phi_j^z  ]^T

R = [ g1 -g2 0;
      g2  g1 0;
       0   0 1]; % notar que transpone al final
O = zeros(3);
R=[ R O; O R];
RT=R';
