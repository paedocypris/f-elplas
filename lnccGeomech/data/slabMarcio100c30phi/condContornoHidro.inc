*numeroPassosTempoP{
        10
}
*tempoTotal{
        1.0
}

*pressao_inicial_fluxo{
     5.0d7
}

*codigos_cond_contorno_fluxo{
         1        51         1         1
      5101      5151         1         1
         0
}
*hgNlvect{
         1
}

# Pressão prescrita nas bordas esquerda e direita
*valores_cond_contorno_fluxo{
         1         2     5.0d7
                         5.0d7
        50         1
      5101         2     2.0d7
                         2.0d7
        50         1
         0
}
