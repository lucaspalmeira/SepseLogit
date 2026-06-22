# Relatório de interpretação dos candidatos a biomarcadores transcriptômicos

Este relatório interpreta os arquivos presentes em `results/biomarkers_multigrupo/` sem converter ProbeSet IDs em símbolos gênicos. A interpretação é conservadora: os probes são tratados como candidatos a biomarcadores ou marcadores transcriptômicos candidatos, não como biomarcadores validados.

As colunas originais disponíveis nos CSVs são `Alpha`, `DeltaMean_Log2`, `P_Welch`, `FDR_Welch`, `AUC_Separation` e `EvidenceScore`. Neste texto, elas são relatadas como `Alpha`, `DeltaMedia`, `pValue`, `FDR`, `AUC` e `ScoreIntegrado`.

## 1. Contraste principal: G5_vs_G1_infeccao_24h

O contraste `G5_vs_G1_infeccao_24h` compara a classe positiva G5, infectado 24 h sem tratamento, contra a classe negativa G1, controle não infectado T0. Assim, `Alpha` positivo indica associação com o grupo infectado 24 h e `Alpha` negativo indica associação com o controle. Como `DeltaMedia` foi calculado como média da classe positiva menos média da classe negativa, valores positivos indicam maior expressão média em G5, e valores negativos indicam maior expressão média em G1.

No Top50 integrado desse contraste, 46 probes têm `Alpha` positivo e 4 têm `Alpha` negativo. As probabilidades leave-one-out da assinatura alpha ficaram entre 0.999995 e 1.000000 para as amostras G5, e entre 4.368e-07 e 3.854e-06 para as amostras G1, compatíveis com separação forte neste conjunto pequeno.

### Principais candidatos no contraste G5 vs G1

| ProbeSetID | contraste | direção da associação | Alpha | DeltaMedia | abs(DeltaMedia) | pValue | FDR | AUC | ScoreIntegrado | interpretação curta |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1421262_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.02840 | 7.767 | 7.767 | 5.891e-05 | 0.005 | 1.00 | 3.543 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1418722_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.02462 | 6.111 | 6.111 | 7.326e-08 | 5.676e-04 | 1.00 | 3.416 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1427747_a_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.02193 | 6.496 | 6.496 | 7.940e-08 | 5.676e-04 | 1.00 | 3.371 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1450009_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.02252 | 5.231 | 5.231 | 1.489e-06 | 0.001 | 1.00 | 3.170 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1450188_s_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.02299 | 6.462 | 6.462 | 1.620e-04 | 0.008 | 1.00 | 3.131 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1437060_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.02106 | 5.286 | 5.286 | 1.337e-06 | 0.001 | 1.00 | 3.125 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1434046_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.02114 | 5.443 | 5.443 | 1.317e-05 | 0.002 | 1.00 | 3.061 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1425451_s_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.01684 | 5.105 | 5.105 | 3.092e-07 | 7.717e-04 | 1.00 | 2.981 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1419764_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.01658 | 5.117 | 5.117 | 1.352e-06 | 0.001 | 1.00 | 2.946 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1419532_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.01675 | 4.986 | 4.986 | 3.555e-06 | 0.001 | 1.00 | 2.906 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1440865_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.01734 | 4.339 | 4.339 | 2.394e-07 | 7.712e-04 | 1.00 | 2.900 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1417290_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.01480 | 4.928 | 4.928 | 1.259e-06 | 0.001 | 1.00 | 2.859 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1424509_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.01865 | 4.913 | 4.913 | 3.707e-05 | 0.004 | 1.00 | 2.852 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1419681_a_at | G5_vs_G1_infeccao_24h | associado ao infectado 24 h | 0.01716 | 5.479 | 5.479 | 4.142e-04 | 0.013 | 1.00 | 2.749 | associado a G5 infectado 24 h sem tratamento; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1422122_at | G5_vs_G1_infeccao_24h | associado ao controle | -0.01433 | -3.829 | 3.829 | 2.522e-04 | 0.011 | 1.00 | 2.462 | associado a G1 controle T0; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1440837_at | G5_vs_G1_infeccao_24h | associado ao controle | -0.01384 | -3.704 | 3.704 | 2.268e-04 | 0.010 | 1.00 | 2.435 | associado a G1 controle T0; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1442023_at | G5_vs_G1_infeccao_24h | associado ao controle | -0.01339 | -3.513 | 3.513 | 1.693e-04 | 0.008 | 1.00 | 2.411 | associado a G1 controle T0; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1450912_at | G5_vs_G1_infeccao_24h | associado ao controle | -0.01459 | -3.589 | 3.589 | 8.389e-04 | 0.019 | 1.00 | 2.378 | associado a G1 controle T0; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1455530_at | G5_vs_G1_infeccao_24h | associado ao controle | -0.01320 | -3.758 | 3.758 | 7.779e-04 | 0.019 | 1.00 | 2.354 | associado a G1 controle T0; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1429889_at | G5_vs_G1_infeccao_24h | associado ao controle | -0.01343 | -3.633 | 3.633 | 0.001 | 0.023 | 1.00 | 2.323 | associado a G1 controle T0; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1423226_at | G5_vs_G1_infeccao_24h | associado ao controle | -0.01231 | -3.137 | 3.137 | 4.617e-04 | 0.014 | 1.00 | 2.272 | associado a G1 controle T0; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1427860_at | G5_vs_G1_infeccao_24h | associado ao controle | -0.01303 | -2.884 | 2.884 | 0.002 | 0.030 | 1.00 | 2.189 | associado a G1 controle T0; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1419769_at | G5_vs_G1_infeccao_24h | associado ao controle | -0.01172 | -2.862 | 2.862 | 8.455e-04 | 0.019 | 1.00 | 2.183 | associado a G1 controle T0; DeltaMedia coerente; FDR baixo; AUC=1,00. |
| 1456328_at | G5_vs_G1_infeccao_24h | associado ao controle | -0.01138 | -2.834 | 2.834 | 9.290e-04 | 0.020 | 1.00 | 2.162 | associado a G1 controle T0; DeltaMedia coerente; FDR baixo; AUC=1,00. |


### Interpretação do contraste principal

Os candidatos mais fortes associados ao infectado 24 h combinam `Alpha` positivo, `DeltaMedia` positivo, FDR baixo e AUC de separação igual a 1,00. Dentro desse conjunto, destacam-se `1421262_at`, `1418722_at`, `1427747_a_at`, `1450009_at`, `1450188_s_at`, `1437060_at` e `1434046_at`. Os probes `1450912_at`, `1422122_at`, `1440837_at`, `1442023_at` e `1455530_at` aparecem como marcadores inversos associados ao controle: eles discriminam o contraste, mas seu sentido biológico é maior expressão média em G1, não em G5.

## 2. Ranking consenso entre contrastes

O arquivo `05_ranking_consenso_biomarcadores.csv` resume recorrência no Top50 entre contrastes, mas não contém a lista de contrastes por probe; por isso, essa lista foi reconstruída a partir dos `top50_biomarcadores_integrado.csv` de cada subpasta.

| ProbeSetID | Frequência Top50 | Contrastes | Melhor FDR | Melhor AUC | Maior abs(DeltaMedia) |
| --- | --- | --- | --- | --- | --- |
| 1418722_at | 5 | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina | 5.676e-04 | 1.00 | 6.111 |
| 1425451_s_at | 5 | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina; G3_vs_G1_vancomicina_sem_infeccao | 7.717e-04 | 1.00 | 5.105 |
| 1419764_at | 5 | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina; G3_vs_G1_vancomicina_sem_infeccao | 0.001 | 1.00 | 5.117 |
| 1419709_at | 4 | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina | 0.004 | 1.00 | 4.963 |
| 1450009_at | 4 | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina | 0.001 | 1.00 | 5.231 |
| 1434046_at | 4 | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida | 0.001 | 1.00 | 5.443 |
| 1419647_a_at | 4 | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G5_efeito_linezolida; G7_vs_G5_efeito_vancomicina | 0.002 | 1.00 | 3.595 |
| 1427747_a_at | 3 | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h | 5.676e-04 | 1.00 | 6.496 |
| 1421262_at | 3 | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida | 0.005 | 1.00 | 7.767 |
| 1434758_at | 3 | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G5_efeito_linezolida | 0.003 | 1.00 | 4.471 |
| 1419532_at | 3 | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G5_efeito_linezolida | 0.001 | 1.00 | 4.986 |
| 1419691_at | 3 | G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina | 0.089 | 1.00 | 2.570 |
| 1435761_at | 3 | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G7_linezolida_vs_vancomicina | 0.002 | 1.00 | 4.337 |
| 1435906_x_at | 3 | G5_vs_G4_progressao_2h_para_24h; G7_vs_G5_efeito_vancomicina; G6_vs_G7_linezolida_vs_vancomicina | 0.123 | 1.00 | 3.725 |
| 1437056_x_at | 3 | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G5_efeito_linezolida | 0.003 | 1.00 | 4.065 |


A recorrência alta deve ser lida com cautela: alguns probes reaparecem tanto em contrastes de infecção quanto em contrastes de tratamento, possivelmente refletindo o eixo principal de infecção/progressão e não um efeito terapêutico específico.

## 3. Candidatos por categoria

### Marcadores de infecção

| ProbeSetID | contrastes nos quais aparece | frequência/recorrência | melhor FDR | melhor AUC | maior abs(DeltaMedia) | direção predominante | justificativa curta |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1449366_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h | 2 | 5.105e-04 | 1.00 | 4.708 | associado aos grupos infectados (G4/G5) | melhor FDR=5.105e-04; melhor AUC=1.00; abs(DeltaMedia) máximo=4.708; frequência global Top50=2. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1427747_a_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h | 2 | 5.676e-04 | 1.00 | 6.496 | associado aos grupos infectados (G4/G5) | melhor FDR=5.676e-04; melhor AUC=1.00; abs(DeltaMedia) máximo=6.496; frequência global Top50=3. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1434758_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G5_efeito_linezolida | 2 | 0.003 | 1.00 | 4.471 | associado aos grupos infectados (G4/G5) | melhor FDR=0.003; melhor AUC=1.00; abs(DeltaMedia) máximo=4.471; frequência global Top50=3. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1419532_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G5_efeito_linezolida | 2 | 0.001 | 1.00 | 4.986 | associado aos grupos infectados (G4/G5) | melhor FDR=0.001; melhor AUC=1.00; abs(DeltaMedia) máximo=4.986; frequência global Top50=3. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1418722_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina | 2 | 5.676e-04 | 1.00 | 6.111 | associado aos grupos infectados (G4/G5) | melhor FDR=5.676e-04; melhor AUC=1.00; abs(DeltaMedia) máximo=6.111; frequência global Top50=5. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1435761_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G7_linezolida_vs_vancomicina | 2 | 0.002 | 1.00 | 4.337 | associado aos grupos infectados (G4/G5) | melhor FDR=0.002; melhor AUC=1.00; abs(DeltaMedia) máximo=4.337; frequência global Top50=3. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1434484_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h | 2 | 8.201e-04 | 1.00 | 4.979 | associado aos grupos infectados (G4/G5) | melhor FDR=8.201e-04; melhor AUC=1.00; abs(DeltaMedia) máximo=4.979; frequência global Top50=2. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1437056_x_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G5_efeito_linezolida | 2 | 0.003 | 1.00 | 4.065 | associado aos grupos infectados (G4/G5) | melhor FDR=0.003; melhor AUC=1.00; abs(DeltaMedia) máximo=4.065; frequência global Top50=3. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1437060_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G5_efeito_linezolida | 2 | 0.001 | 1.00 | 5.286 | associado aos grupos infectados (G4/G5) | melhor FDR=0.001; melhor AUC=1.00; abs(DeltaMedia) máximo=5.286; frequência global Top50=3. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1422953_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h | 2 | 8.201e-04 | 1.00 | 3.996 | associado aos grupos infectados (G4/G5) | melhor FDR=8.201e-04; melhor AUC=1.00; abs(DeltaMedia) máximo=3.996; frequência global Top50=2. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1434046_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida | 2 | 0.001 | 1.00 | 5.443 | associado aos grupos infectados (G4/G5) | melhor FDR=0.001; melhor AUC=1.00; abs(DeltaMedia) máximo=5.443; frequência global Top50=4. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1419647_a_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G5_efeito_linezolida; G7_vs_G5_efeito_vancomicina | 2 | 0.002 | 1.00 | 3.595 | associado aos grupos infectados (G4/G5) | melhor FDR=0.002; melhor AUC=1.00; abs(DeltaMedia) máximo=3.595; frequência global Top50=4. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1421262_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida | 1 | 0.005 | 1.00 | 7.767 | associado a G5 infectado 24 h sem tratamento | melhor FDR=0.005; melhor AUC=1.00; abs(DeltaMedia) máximo=7.767; frequência global Top50=3. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1450009_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina | 1 | 0.001 | 1.00 | 5.231 | associado a G5 infectado 24 h sem tratamento | melhor FDR=0.001; melhor AUC=1.00; abs(DeltaMedia) máximo=5.231; frequência global Top50=4. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |
| 1450188_s_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h | 1 | 0.008 | 1.00 | 6.462 | associado a G5 infectado 24 h sem tratamento | melhor FDR=0.008; melhor AUC=1.00; abs(DeltaMedia) máximo=6.462; frequência global Top50=2. Direção positiva indica maior expressão no grupo infectado; direção negativa indica marcador inverso associado ao controle. |


### Marcadores de progressão

| ProbeSetID | contrastes nos quais aparece | frequência/recorrência | melhor FDR | melhor AUC | maior abs(DeltaMedia) | direção predominante | justificativa curta |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1421262_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida | 1 | 0.026 | 1.00 | 7.514 | associado a G5 infectado 24 h sem tratamento | melhor FDR=0.026; melhor AUC=1.00; abs(DeltaMedia) máximo=7.514; frequência global Top50=3. Alpha positivo sugere aumento em 24 h; alpha negativo sugere maior expressão em 2 h. |
| 1450188_s_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h | 1 | 0.038 | 1.00 | 6.309 | associado a G5 infectado 24 h sem tratamento | melhor FDR=0.038; melhor AUC=1.00; abs(DeltaMedia) máximo=6.309; frequência global Top50=2. Alpha positivo sugere aumento em 24 h; alpha negativo sugere maior expressão em 2 h. |
| 1450009_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina | 1 | 0.015 | 1.00 | 5.122 | associado a G5 infectado 24 h sem tratamento | melhor FDR=0.015; melhor AUC=1.00; abs(DeltaMedia) máximo=5.122; frequência global Top50=4. Alpha positivo sugere aumento em 24 h; alpha negativo sugere maior expressão em 2 h. |
| 1425451_s_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina; G3_vs_G1_vancomicina_sem_infeccao | 1 | 0.006 | 1.00 | 4.475 | associado a G5 infectado 24 h sem tratamento | melhor FDR=0.006; melhor AUC=1.00; abs(DeltaMedia) máximo=4.475; frequência global Top50=5. Alpha positivo sugere aumento em 24 h; alpha negativo sugere maior expressão em 2 h. |
| 1419764_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina; G3_vs_G1_vancomicina_sem_infeccao | 1 | 0.012 | 1.00 | 4.598 | associado a G5 infectado 24 h sem tratamento | melhor FDR=0.012; melhor AUC=1.00; abs(DeltaMedia) máximo=4.598; frequência global Top50=5. Alpha positivo sugere aumento em 24 h; alpha negativo sugere maior expressão em 2 h. |
| 1419681_a_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h | 1 | 0.067 | 1.00 | 4.977 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.067); melhor AUC=1.00; abs(DeltaMedia) máximo=4.977; frequência global Top50=2. Alpha positivo sugere aumento em 24 h; alpha negativo sugere maior expressão em 2 h. |
| 1418722_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina | 1 | 0.012 | 1.00 | 3.392 | associado a G5 infectado 24 h sem tratamento | melhor FDR=0.012; melhor AUC=1.00; abs(DeltaMedia) máximo=3.392; frequência global Top50=5. Alpha positivo sugere aumento em 24 h; alpha negativo sugere maior expressão em 2 h. |
| 1436530_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G3_vs_G1_vancomicina_sem_infeccao | 1 | 0.049 | 1.00 | 4.664 | associado a G5 infectado 24 h sem tratamento | melhor FDR=0.049; melhor AUC=1.00; abs(DeltaMedia) máximo=4.664; frequência global Top50=3. Alpha positivo sugere aumento em 24 h; alpha negativo sugere maior expressão em 2 h. |
| 1450826_a_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h | 1 | 0.045 | 1.00 | 4.314 | associado a G5 infectado 24 h sem tratamento | melhor FDR=0.045; melhor AUC=1.00; abs(DeltaMedia) máximo=4.314; frequência global Top50=2. Alpha positivo sugere aumento em 24 h; alpha negativo sugere maior expressão em 2 h. |
| 1448562_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h | 1 | 0.071 | 1.00 | 4.567 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.071); melhor AUC=1.00; abs(DeltaMedia) máximo=4.567; frequência global Top50=2. Alpha positivo sugere aumento em 24 h; alpha negativo sugere maior expressão em 2 h. |
| 1450428_at | G4_vs_G1_infeccao_2h; G5_vs_G4_progressao_2h_para_24h | 1 | 0.021 | 1.00 | 2.690 | associado a G4 infectado 2 h | melhor FDR=0.021; melhor AUC=1.00; abs(DeltaMedia) máximo=2.690; frequência global Top50=2. Alpha positivo sugere aumento em 24 h; alpha negativo sugere maior expressão em 2 h. |
| 1429274_at | G5_vs_G4_progressao_2h_para_24h | 1 | 0.006 | 1.00 | 1.494 | associado a G4 infectado 2 h | melhor FDR=0.006; melhor AUC=1.00; abs(DeltaMedia) máximo=1.494; frequência global Top50=1. Alpha positivo sugere aumento em 24 h; alpha negativo sugere maior expressão em 2 h. |


### Marcadores de resposta à linezolida

| ProbeSetID | contrastes nos quais aparece | frequência/recorrência | melhor FDR | melhor AUC | maior abs(DeltaMedia) | direção predominante | justificativa curta |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1450009_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina | 1 | 0.628 | 0.84 | 2.342 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.628); melhor AUC=0.84; abs(DeltaMedia) máximo=2.342; frequência global Top50=4. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1419691_at | G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina | 1 | 0.626 | 0.92 | 2.185 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.626); melhor AUC=0.92; abs(DeltaMedia) máximo=2.185; frequência global Top50=3. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1419709_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina | 1 | 0.626 | 0.92 | 2.063 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.626); melhor AUC=0.92; abs(DeltaMedia) máximo=2.063; frequência global Top50=4. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1425451_s_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina; G3_vs_G1_vancomicina_sem_infeccao | 1 | 0.626 | 1.00 | 1.780 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.626); melhor AUC=1.00; abs(DeltaMedia) máximo=1.780; frequência global Top50=5. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1419532_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G5_efeito_linezolida | 1 | 0.626 | 1.00 | 1.975 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.626); melhor AUC=1.00; abs(DeltaMedia) máximo=1.975; frequência global Top50=3. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1419764_at | G5_vs_G1_infeccao_24h; G5_vs_G4_progressao_2h_para_24h; G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina; G3_vs_G1_vancomicina_sem_infeccao | 1 | 0.626 | 1.00 | 1.651 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.626); melhor AUC=1.00; abs(DeltaMedia) máximo=1.651; frequência global Top50=5. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1434758_at | G4_vs_G1_infeccao_2h; G5_vs_G1_infeccao_24h; G6_vs_G5_efeito_linezolida | 1 | 0.626 | 1.00 | 1.548 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.626); melhor AUC=1.00; abs(DeltaMedia) máximo=1.548; frequência global Top50=3. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1416125_at | G4_vs_G1_infeccao_2h; G6_vs_G5_efeito_linezolida | 1 | 0.626 | 1.00 | 1.510 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.626); melhor AUC=1.00; abs(DeltaMedia) máximo=1.510; frequência global Top50=2. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1425814_a_at | G6_vs_G5_efeito_linezolida; G7_vs_G5_efeito_vancomicina | 1 | 0.626 | 1.00 | 1.924 | associado a G6 infectado 24 h + linezolida | FDR elevado (melhor=0.626); melhor AUC=1.00; abs(DeltaMedia) máximo=1.924; frequência global Top50=2. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1449434_at | G6_vs_G5_efeito_linezolida; G6_vs_G7_linezolida_vs_vancomicina | 1 | 0.635 | 0.96 | 1.505 | associado a G6 infectado 24 h + linezolida | FDR elevado (melhor=0.635); melhor AUC=0.96; abs(DeltaMedia) máximo=1.505; frequência global Top50=2. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |


### Marcadores de resposta à vancomicina

| ProbeSetID | contrastes nos quais aparece | frequência/recorrência | melhor FDR | melhor AUC | maior abs(DeltaMedia) | direção predominante | justificativa curta |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1418930_at | G7_vs_G5_efeito_vancomicina; G6_vs_G7_linezolida_vs_vancomicina | 1 | 0.478 | 0.88 | 2.230 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.478); melhor AUC=0.88; abs(DeltaMedia) máximo=2.230; frequência global Top50=2. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1427381_at | G5_vs_G1_infeccao_24h; G7_vs_G5_efeito_vancomicina; G6_vs_G7_linezolida_vs_vancomicina | 1 | 0.452 | 0.84 | 2.153 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.452); melhor AUC=0.84; abs(DeltaMedia) máximo=2.153; frequência global Top50=3. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1452349_x_at | G7_vs_G5_efeito_vancomicina | 1 | 0.449 | 0.96 | 1.412 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.449); melhor AUC=0.96; abs(DeltaMedia) máximo=1.412; frequência global Top50=1. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1437226_x_at | G7_vs_G5_efeito_vancomicina | 1 | 0.449 | 0.96 | 1.698 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.449); melhor AUC=0.96; abs(DeltaMedia) máximo=1.698; frequência global Top50=1. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1435906_x_at | G5_vs_G4_progressao_2h_para_24h; G7_vs_G5_efeito_vancomicina; G6_vs_G7_linezolida_vs_vancomicina | 1 | 0.492 | 0.80 | 1.990 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.492); melhor AUC=0.80; abs(DeltaMedia) máximo=1.990; frequência global Top50=3. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1417300_at | G7_vs_G5_efeito_vancomicina | 1 | 0.479 | 0.88 | 1.569 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.479); melhor AUC=0.88; abs(DeltaMedia) máximo=1.569; frequência global Top50=1. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1428579_at | G7_vs_G5_efeito_vancomicina | 1 | 0.448 | 0.96 | 1.425 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.448); melhor AUC=0.96; abs(DeltaMedia) máximo=1.425; frequência global Top50=1. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1460282_at | G7_vs_G5_efeito_vancomicina | 1 | 0.448 | 0.96 | 1.623 | associado a G5 infectado 24 h sem tratamento | FDR elevado (melhor=0.448); melhor AUC=0.96; abs(DeltaMedia) máximo=1.623; frequência global Top50=1. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1425814_a_at | G6_vs_G5_efeito_linezolida; G7_vs_G5_efeito_vancomicina | 1 | 0.452 | 0.96 | 1.387 | associado a G7 infectado 24 h + vancomicina | FDR elevado (melhor=0.452); melhor AUC=0.96; abs(DeltaMedia) máximo=1.387; frequência global Top50=2. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |
| 1415682_at | G7_vs_G5_efeito_vancomicina; G3_vs_G1_vancomicina_sem_infeccao | 1 | 0.448 | 1.00 | 1.307 | associado a G7 infectado 24 h + vancomicina | FDR elevado (melhor=0.448); melhor AUC=1.00; abs(DeltaMedia) máximo=1.307; frequência global Top50=2. Candidato exploratório: neste contraste os FDRs dos top probes são altos, apesar de AUC/score úteis para ranqueamento. |


Nos contrastes de tratamento (`G6_vs_G5_efeito_linezolida` e `G7_vs_G5_efeito_vancomicina`), os principais probes têm FDRs elevados. Por isso, eles devem ser considerados sinais exploratórios de resposta terapêutica, não candidatos robustos no mesmo nível dos contrastes de infecção e progressão.

## 4. Qualidade dos classificadores

| Contraste | N | Acc completo | Acc reduzido | LOOCV alpha | KNN LOOCV | SVM LOOCV | Leitura conservadora |
| --- | --- | --- | --- | --- | --- | --- | --- |
| G4_vs_G1_infeccao_2h | 5+5 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | mais confiável dentro deste conjunto pequeno |
| G5_vs_G1_infeccao_24h | 5+5 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | mais confiável dentro deste conjunto pequeno |
| G5_vs_G4_progressao_2h_para_24h | 5+5 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | mais confiável dentro deste conjunto pequeno |
| G6_vs_G5_efeito_linezolida | 5+5 | 1.00 | 1.00 | 0.50 | 0.90 | 0.80 | interpretar com cautela; sinais de instabilidade entre validadores |
| G7_vs_G5_efeito_vancomicina | 5+5 | 1.00 | 1.00 | 0.70 | 0.80 | 0.90 | interpretar com cautela; sinais de instabilidade entre validadores |
| G6_vs_G7_linezolida_vs_vancomicina | 5+5 | 1.00 | 1.00 | 0.90 | 0.80 | 0.90 | interpretar com cautela; sinais de instabilidade entre validadores |
| G2_vs_G1_linezolida_sem_infeccao | 5+5 | 1.00 | 1.00 | 0.70 | 0.50 | 0.80 | interpretar com cautela; sinais de instabilidade entre validadores |
| G3_vs_G1_vancomicina_sem_infeccao | 5+5 | 1.00 | 1.00 | 0.60 | 0.90 | 1.00 | interpretar com cautela; sinais de instabilidade entre validadores |


Os contrastes `G4_vs_G1_infeccao_2h`, `G5_vs_G1_infeccao_24h` e `G5_vs_G4_progressao_2h_para_24h` são os mais consistentes dentro dos arquivos disponíveis: modelo completo, modelo reduzido e validações leave-one-out retornam acurácia 1,00. Mesmo assim, o resultado deve ser interpretado com cautela porque cada contraste usa apenas 10 amostras no total.

Os contrastes de tratamento apresentam maior instabilidade. Em `G6_vs_G5_efeito_linezolida`, a assinatura alpha em leave-one-out cai para 0,50, embora kNN e SVM fiquem acima disso. Em `G7_vs_G5_efeito_vancomicina`, a assinatura alpha fica em 0,70. Essa divergência entre classificadores e os FDRs altos indicam maior risco de overfitting e menor robustez estatística.

## 5. Metodologia dos scripts

A análise computacional foi conduzida a partir da matriz `data/gse_data_no_header.txt`, carregada como tabela tabulada sem cabeçalho explícito. A primeira coluna foi interpretada como identificador do gene/probe, aqui mantido como `ProbeSetID`, e as demais colunas foram convertidas para uma matriz numérica de expressão. A estrutura da matriz é de genes/probes nas linhas e amostras nas colunas, totalizando 35 amostras distribuídas em sete grupos experimentais com cinco amostras por grupo.

Os grupos experimentais foram definidos conforme o desenho do GSE38531: G1, controle não infectado T0; G2, não infectado tratado com linezolida; G3, não infectado tratado com vancomicina; G4, infectado por 2 h; G5, infectado por 24 h sem tratamento; G6, infectado por 24 h tratado com linezolida; e G7, infectado por 24 h tratado com vancomicina. Antes das análises, foram removidos genes/probes com pelo menos um valor `NaN` e genes/probes sem variação entre as amostras, evitando instabilidade numérica e testes sem informação discriminativa.

Para visualização global da estrutura dos dados, foi aplicada uma decomposição SVD/PCA sobre os genes/probes de maior variabilidade. A matriz foi centralizada por linha, de modo que os escores das componentes principais representassem padrões amostrais globais entre os sete grupos. Os escores foram salvos em `01_PCA_global_scores.csv`, e as figuras correspondentes foram geradas em `figures/biomarkers_multigrupo/`.

Em seguida, foi realizada ANOVA univariada entre os sete grupos para cada gene/probe. A estatística F, o p-valor e a correção por FDR pelo procedimento de Benjamini-Hochberg foram salvos em `02_ANOVA_7_grupos_todos_genes.csv`. Essa etapa fornece uma visão global de variação entre grupos, mas não substitui os contrastes biologicamente orientados.

Os scripts criaram contrastes binários definidos a priori: infecção precoce (`G4_vs_G1_infeccao_2h`), infecção de 24 h (`G5_vs_G1_infeccao_24h`), progressão de 2 h para 24 h (`G5_vs_G4_progressao_2h_para_24h`), resposta à linezolida (`G6_vs_G5_efeito_linezolida`), resposta à vancomicina (`G7_vs_G5_efeito_vancomicina`), comparação linezolida versus vancomicina (`G6_vs_G7_linezolida_vs_vancomicina`) e efeitos dos antibióticos em animais não infectados (`G2_vs_G1_linezolida_sem_infeccao` e `G3_vs_G1_vancomicina_sem_infeccao`). Em todos os contrastes, a classe positiva foi definida pelo primeiro grupo do nome do contraste e a classe negativa pelo segundo grupo.

Para cada contraste, foi ajustada uma regressão logística modificada. Os rótulos binários foram transformados em valores logit extremos, e os coeficientes `alpha` foram obtidos pela função `resolve` quando disponível no caminho do MATLAB/Octave; quando a função não estava disponível, o script utilizou a formulação equivalente por sistema linear aumentado. Um coeficiente `alpha` positivo indica associação com a classe positiva do contraste, enquanto um coeficiente `alpha` negativo indica associação com a classe negativa. Essa interpretação foi usada para separar probes associados ao grupo infectado, ao controle, à progressão temporal ou aos grupos tratados.

Os genes/probes com maiores coeficientes positivos e negativos foram selecionados para formar assinaturas reduzidas. Para cada contraste, o script salvou `top10_alpha_positivo.csv` e `top10_alpha_negativo.csv`, além de calcular probabilidades preditas para o modelo completo, para o modelo reduzido e para uma validação leave-one-out baseada na assinatura alpha. Classificadores kNN e SVM também foram avaliados sobre a assinatura reduzida, com acurácia aparente e leave-one-out reportadas em `03_resumo_contrastes_classificacao.csv`.

Paralelamente, para cada gene/probe e contraste, foram calculadas as médias de expressão nas classes positiva e negativa, a diferença média `DeltaMean_Log2` definida como média da classe positiva menos média da classe negativa, teste t de Welch bicaudal, p-valor, FDR de Benjamini-Hochberg, AUC univariada considerando a classe positiva mais alta e AUC de separação independente do sinal. Esses resultados foram reunidos em `todos_genes_ranqueados.csv` dentro de cada pasta de contraste.

A priorização integrada de candidatos combinou evidências normalizadas de magnitude de `alpha`, magnitude absoluta de `DeltaMean_Log2`, significância por FDR e capacidade de separação por AUC. O resultado foi o `EvidenceScore`, aqui interpretado como `ScoreIntegrado`. Os 50 probes de maior score em cada contraste foram salvos em `top50_biomarcadores_integrado.csv`.

Por fim, os scripts agregaram os 50 principais probes de todos os contrastes em `04_top_biomarcadores_todos_contrastes.csv` e construíram um ranking consenso em `05_ranking_consenso_biomarcadores.csv`, contendo frequência de recorrência no Top50, score médio, score máximo, maior `abs(alpha)`, maior `abs(DeltaMean_Log2)`, menor FDR e maior AUC de separação. A recorrência entre contrastes foi usada apenas como critério de priorização exploratória, não como validação independente.

Todos os resultados devem ser interpretados como geração de candidatos a biomarcadores transcriptômicos. O conjunto possui apenas cinco amostras por grupo, o que aumenta o risco de overfitting, sobretudo em modelos multivariados com muitos genes/probes. Portanto, os ProbeSet IDs priorizados requerem anotação posterior, validação externa e, idealmente, confirmação experimental antes de qualquer interpretação como biomarcadores definitivos.

## 6. Limitações

A principal limitação é o tamanho amostral: há apenas cinco amostras por grupo. A separação perfeita em alguns contrastes pode refletir sinal biológico forte, mas também pode ser amplificada por alta dimensionalidade e seleção de probes no mesmo conjunto. Além disso, o ranking usa ProbeSet IDs sem anotação gênica nesta etapa, e diferentes probes podem mapear para o mesmo gene ou para regiões com interpretação ambígua. Os candidatos precisam de anotação posterior, validação externa e confirmação experimental.

## 7. ProbeSet IDs prioritários para anotação posterior

Prioridade principal, associados à infecção/progressão com direção positiva em G5 quando disponível: `1421262_at`, `1418722_at`, `1427747_a_at`, `1450009_at`, `1450188_s_at`, `1437060_at`, `1434046_at`, `1419532_at`, `1449366_at`, `1434758_at`, `1425451_s_at`, `1419764_at`.

Marcadores inversos associados ao controle no contraste G5 vs G1: `1450912_at`, `1422122_at`, `1440837_at`, `1442023_at`, `1455530_at`.
