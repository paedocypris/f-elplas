*pressao_inicial_fluxo{
     0.0d0
}

*codigos_cond_contorno_fluxo{
         1         1         1         1
         0
}
*hgNlvect{
         1
}

# Pressão prescrita nas bordas esquerda e direita
*valores_cond_contorno_fluxo{
         1         2     0.0D0
                         0.0D0
         1         1
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
        1           1.5d3
        1           3.5d3
        1           5.0d3
        1           3.0d4
        1           6.0d4
        1           7.0d5
        1           2.0d5
        1           1.0d6
        1           1.0d6
        1           2.0d6
}
