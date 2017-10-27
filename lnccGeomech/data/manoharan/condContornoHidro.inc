*pressao_inicial_fluxo{
     0.0d0
}

*codigos_cond_contorno_fluxo{
         1        76        25         0
         1        25         1         0
        25      1225        25         0
      1201      1225         1         0
        76      1201        25         1
         0
}
*hgNlvect{
         1
}

# Pressão prescrita na parte sem carregamento
*valores_cond_contorno_fluxo{
        76         2     0.0d0
                         0.0d0
        46         1
         0
}

*hgTimeIsTable{
   1
}

*tabelaTempos{
        19
        10        0.468D0
         5        0.936D0
         6        2.340D0
         5        4.680D0
         5        9.360D0
         6        2.340D1
         5        4.680D1
         5        9.360D1
         6        2.340D2
         5        4.680D2
         5        9.360D2
         6        2.340D3
         5        4.680D3
         5        9.360D3
         6        2.340D4
         5        4.680D4
         5        9.360D4
         6        2.340D5
         5        4.680D5
}
