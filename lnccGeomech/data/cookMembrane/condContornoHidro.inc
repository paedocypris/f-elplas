*numeroPassosTempoP{
        50
}
*tempoTotal{
      12.0
}
*hgCompressibilidade{
     1.d-7
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

# Press�o prescrita nas bordas esquerda e direita
*valores_cond_contorno_fluxo{
         1         2     6.0d7
                         6.0d7
        50       101
       101         2     3.5d7
                         3.5d7
        50       101
         0
}

#impressao dos dados
*impressao_da_pressaoGalerkin{   
1
./out/
1
reservatorio  
}