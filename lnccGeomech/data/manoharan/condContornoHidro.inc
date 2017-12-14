*pressao_inicial_fluxo{
     0.0d0
}
*hgNlvect{
         1
}

*codigos_cond_contorno_fluxo1{
         1        76        25         0
         1        25         1         0
        25      1225        25         0
      1201      1225         1         0
        76      1201        25         1
         0
}

# Pressão prescrita na parte sem carregamento
*valores_cond_contorno_fluxo1{
        76         2     0.0d0
                         0.0d0
        46         1
         0
}
