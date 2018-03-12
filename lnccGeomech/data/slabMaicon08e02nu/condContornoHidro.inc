
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
         1         2   2.32D-5
                       2.32D-5
        15         1
         1         0   1.16D-5
        16         0   1.16D-5
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
           0.d0           0.d0
           1.d0           0.d0
           2.d0           0.d0
          30.d0           0.d0
       30.001d0           1.d0
          31.d0           1.d0
          32.d0           1.d0
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
        4
        10           1.d0
         5           4.d0
        10           1.D0
        13         20.0d0
}
