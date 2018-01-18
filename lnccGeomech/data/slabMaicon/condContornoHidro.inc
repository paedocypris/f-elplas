
*pressao_inicial_fluxo{
     6.0d7
}

*codigos_cond_contorno_fluxo{
         1        16         1         0
       481       496         1         1
         0
}
*hgNlvect{
         2
}

# Pressão prescrita nas bordas esquerda e direita
*valores_cond_contorno_fluxo{
         1         2   6.94D-5
                       6.94D-5
        15         1
         1         0   3.47D-5
        16         0   3.47D-5
         0
       481         2     1.0D7
                         1.0D7
        15         1
         0
}

*hgNltftn{
   2
}

*hgNptslf{
   7
}

*hgLoadTimeFunction{
           0.d0           1.d0
           1.d0           1.d0
           2.d0           1.d0
           3.d0           1.d0
           4.d0           1.d0
           5.d0           1.d0
           6.d0           1.d0
           0.D0           6.D0
           6.D0         3.58D0
          12.D0         2.07D0
          18.D0         1.40D0
          24.D0         1.14D0
          30.D0         1.05D0
          50.d0         1.00d0
}

*hgTimeIsTable{
   1
}

*tabelaTempos{
        3
        10           1.d0
         5           4.d0
        10         27.0D0
}
