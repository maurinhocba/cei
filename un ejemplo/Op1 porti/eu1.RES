------------------------------------------------------------------------
    Programa PORTI      autor: Julio C. MASSA      versi�n: feb 2005

                      AN�LISIS DE P�RTICOS PLANOS

  Dto. Estructuras - F.C.E.F. y N. - Univ. Nacional de C�rdoba -  Arg.
------------------------------------------------------------------------
                                     fecha:  1/ 6/2024       hora: 12:14

EU1                                                                     

  4 NUDOS   3 BARRAS    4 APOYOS   0 INCLIN    2 SECC   1 MAT    1 HIPOT

BARRA          NUDO I     NUDO J        SECCI�N   MATERIAL      EXTREMOS
   1              1         2              1         1              0
   2              2         3              1         1              0
   3              3         4              2         1              3

 NUDO          COOR X         COOR Y
   1           0.0000         0.0000
   2         100.0000         0.0000
   3         500.0000         0.0000
   4         500.5000        -0.8660

SECC         AREA       INERCIA        Acorte         Hsuper      Hinfer
  1       100.000       833.330  ************          5.000       5.000
  2         0.000       833.330  ************          5.000       5.000

MATERIAL         E                 G               alfa             peso
  1          210000.0           81000.0         0.0000100         0.0000

HIPOTESIS DE CARGA NUMERO 1                                             

APOYOS -------->  codigo  ( 0 LIBRE   1 NULO   2 RESORTE   3 PREFIJADO )
NUDO      cod X  cod Y  cod Z           DATO X       DATO Y       DATO Z 
  1         0      1      0            0.00000      0.00000      0.00000
  2         0      1      0            0.00000      0.00000      0.00000
  3         3      3      0            1.00000      0.00000      0.00000
  4         1      1      0            0.00000      0.00000      0.00000
------------------------------------------------------------------------

RESULTADOS HIP�TESIS N�MERO  1

BARRA  NUDO        AXIAL         CORTE         FLECTOR           TENSION
   1     1         0.000         0.000            0.00               0.0
         2         0.000         0.000            0.00               0.0
   2     2         0.000         0.000            0.00               0.0
         3         0.000         0.000            0.00               0.0
   3     3        10.001         0.000            0.00          105004.6
         4       -10.001         0.000            0.00          105004.6

       D E S P L A Z A M I E N T O S      ---- R E A C C I O N E S -----
NUDO     DespX      DespY   Rotaci�n      FuerzaX   FuerzaY      MOMENTO
  1    1.00000    0.00000  0.0000000         0.00      0.00          0.0
  2    1.00000    0.00000  0.0000000         0.00      0.00          0.0
  3    1.00000    0.00000  0.0000000         5.00     -8.66          0.0
  4    0.00000    0.00000  0.0000000        -5.00      8.66          0.0
------------------------------------------------------------------------
