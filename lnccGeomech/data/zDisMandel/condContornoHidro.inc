*pressao_inicial_fluxo{
     0.0d0
}

*codigos_cond_contorno_fluxo{
       181       186         1         1
         0
}
*hgNlvect{
         1
}

# Pressão prescrita nas bordas esquerda e direita
*valores_cond_contorno_fluxo{
       181         2     0.0D0
                         0.0D0
         5         1
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
         4
         1       1.00d+03
         1       9.00d+03
         1       9.00d+04
         1       9.00d+05
}
