*pressao_inicial_fluxo{
     0.0d0
}

*codigos_cond_contorno_fluxo{
         1        16         1         0
      1441      1456         1         1
         1      1441        16         0
        16      1456        16         0
         0
}
*hgNlvect{
         1
}

# Press�o prescrita nas bordas esquerda e direita
*valores_cond_contorno_fluxo{
      1441         2     0.0D0
                         0.0D0
        15         1
         0
}

*hgNltftn{
   0
}

*hgNptslf{
   2
}

*hgLoadTimeFunction{
           0.d0           1.d0
           1.d0           1.d0
}

*hgTimeIsTable{
   1
}

*hgTimeConvertion{
   1.0d0
}

*tabelaTempos{
       10
        1           1.0d3
        1           4.0d3
        1           5.0d3
        1           9.0d4
        1           4.0d5
        1           3.0d5
        1           2.0d5
        1           1.0d6
        1           1.0d6
        1           2.0d6
}