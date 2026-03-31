# SepseLogit

Repositório contendo scripts desenvolvidos na disciplina Algoritmos para Bioinformática I (UFMG), com foco na aplicação de técnicas de Aprendizado de Máquina em dados de expressão gênica provenientes de experimentos de microarray disponíveis no NCBI GEO.

Este repositório tem finalidade acadêmica e foi desenvolvido como parte das atividades práticas da disciplina. Também o repositório tem por objetivo organizar o ambiente computacional em MATLAB e consolidar scripts para:

- processamento de dados de microarray
- construção de modelos de classificação
- seleção de atributos (feature selection)
- análise de dados biológicos em formato vector space model

### Dataset utilizado

O conjunto de dados utilizado nas atividades é o:

- GSE38531

Este dataset contém dados de expressão gênica obtidos por microarray, permitindo estudar diferenças de expressão entre condições biológicas distintas.

O problema é modelado como um vector space model, onde:

- Entidades → amostras biológicas (pacientes ou condições experimentais)
- Atributos → níveis de expressão gênica
- Objetivo → classificar as entidades em categorias (ex.: controle vs condição patológica)

### Estrutura do repositório

```
.
├── data/
│   └── GSE38531_series_matrix.txt
│
├── scripts/
│   ├── load_matrix.m
│   ├── preprocess_matrix.m
│   ├── logistic_model.m
│   ├── feature_selection.m
│   └── reduced_model.m
│
├── figures/
│
└── README.md
```