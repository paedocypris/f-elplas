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
         1         2     0.0d0
                         0.0d0
        15         1
         1         0     0.0d0
        16         0     0.0d0
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
   2
}

*hgLoadTimeFunction{
           0.d0           1.d0
           1.d0           1.d0
           0.D0           6.D0
          1.D-3           1.D0
}

*hgTimeIsTable{
   1
}

*hgTimeConvertion{
 86400.0d0
}

*tabelaTempos{
        1
        30          2.0d0
}
