------------------------------------------------------------------------
    Programa PORTI      autor: Julio C. MASSA      versión: feb 2005

                      ANÁLISIS DE PÓRTICOS PLANOS

  Dto. Estructuras - F.C.E.F. y N. - Univ. Nacional de Córdoba -  Arg.
------------------------------------------------------------------------
                                     fecha:  1/ 6/2024       hora: 12:31

EU1                                                                     

  4 NUDOS   3 BARRAS    4 APOYOS   0 INCLIN    2 SECC   1 MAT    1 HIPOT

BARRA          NUDO I     NUDO J        SECCIÓN   MATERIAL      EXTREMOS
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
  3         3      3      0           -0.42750      1.00000      0.00000
  4         1      1      0            0.00000      0.00000      0.00000
------------------------------------------------------------------------

RESULTADOS HIPÓTESIS NÚMERO  1

BARRA  NUDO        AXIAL         CORTE         FLECTOR           TENSION
   1     1         0.000        26.250            0.00               0.0
         2         0.000       -26.250         2624.99              15.7
   2     2         0.000        -6.562        -2624.99              15.7
         3         0.000         6.562            0.00               0.0
   3     3       -21.596         0.000            0.00          226757.5
         4        21.596         0.000            0.00          226757.5

       D E S P L A Z A M I E N T O S      ---- R E A C C I O N E S -----
NUDO     DespX      DespY   Rotación      FuerzaX   FuerzaY      MOMENTO
  1   -0.42750    0.00000 -0.0002500         0.00     26.25          0.0
  2   -0.42750    0.00000  0.0005000         0.00    -32.81          0.0
  3   -0.42750    1.00000  0.0035000       -10.80     25.27          0.0
  4    0.00000    0.00000  0.0000000        10.80    -18.70          0.0
------------------------------------------------------------------------
