*numeroPassosTempoP{
        50
}
*tempoTotal{
      12.0
}

*pressao_inicial_fluxo{
     5.0d7
}

*codigos_cond_contorno_fluxo{
         1        11         5         1
         5        15         5         1
         0
}
*hgNlvect{
         1
}

# Pressão prescrita nas bordas esquerda e direita
*valores_cond_contorno_fluxo{
         1         2     6.0d7
                         6.0d7
         2         5
       101         2     3.5d7
                         3.5d7
         2         5
         0
}

#impressao dos dados
*impressao_da_pressaoGalerkin{   
1
./out/
1
reservatorio  
}
