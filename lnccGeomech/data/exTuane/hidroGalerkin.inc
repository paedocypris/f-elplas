*numeroPassosTempoP{
        50
}
*tempoTotal{
0.00192901
}
*hgCompressibilidade{
     1.d-7
}

*codigos_cond_contorno_fluxo{
         1       111        11         1
        11       121        11         1
         0 
}
*hgNlvect{
         1
}

# Pressão prescrita nas bordas esquerda e direita
*valores_cond_contorno_fluxo{
         1         2    1.0d+6
                        1.0d+6
        10        11
        11         2    0.0d+0
                        0.0d+0
        10        11
         0
}

*pressao_inicial_fluxo{
     1.e+6
}

*impressao_da_pressaoGalerkin{   
1
./out/
1
reservatorio  
}
