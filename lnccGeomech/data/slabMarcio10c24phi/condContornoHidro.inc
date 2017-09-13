*numeroPassosTempoP{
        10
}
*tempoTotal{
      36.0
}

*pressao_inicial_fluxo{
     5.0d7
}

*codigos_cond_contorno_fluxo{
         1      5051       101         1
       101      5151       101         1
         0
}
*hgNlvect{
         1
}

# Pressão prescrita nas bordas esquerda e direita
*valores_cond_contorno_fluxo{
         1         2     7.0d7
                         7.0d7
        50       101
       101         2     3.5d7
                         3.5d7
        50       101
         0
}
