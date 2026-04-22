%% =========================================================
% GSE38531 - Classificação binária com regressão logística modificada
% Disciplina: Algoritmos para Bioinformática I
%
% MATRIZ:
%   A -> linhas = genes/atributos
%        colunas = amostras/entidades
%
% COMPARAÇÃO ESCOLHIDA:
%   Classe 1 = infectado 24 h sem tratamento
%   Classe 0 = controle não infectado
%
% IMPORTANTE:
% 1) Se A já estiver no workspace, o script usa A diretamente
% 2) Se A não existir, o script carrega do arquivo de entrada
% 3) A primeira coluna do arquivo é usada apenas para gene_ids
% 4) As colunas da segunda em diante formam a matriz A
%
% ORDEM REAL DAS 35 AMOSTRAS:
% Para cada réplica biológica:
%   1 = T0 controle
%   2 = T2h infectado
%   3 = T0 + linezolida
%   4 = T0 + vancomicina
%   5 = T24h infectado
%   6 = T24h infectado + linezolida
%   7 = T24h infectado + vancomicina
%
% Repetido para rep1, rep2, rep3, rep4, rep5
%% =========================================================

clearvars -except A gene_ids
clc
close all

%% ---------------------------------------------------------
% 0) Verificar se A existe
%    Se não existir, carregar do arquivo:
%    - coluna 1 = nomes/IDs dos genes
%    - colunas 2:end = matriz numérica A
%% ---------------------------------------------------------
if ~exist('A', 'var')

    % Caminho do arquivo de entrada
    inputFile = '~/Documentos/algoritmos/gse_no_head_data.txt';

    % Expandir "~" para a HOME do usuário
    if startsWith(inputFile, '~/')
        inputFile = fullfile(getenv('HOME'), inputFile(3:end));
    end

    fprintf('A variável A não está no workspace.\n');
    fprintf('Carregando matriz a partir do arquivo:\n%s\n', inputFile);

    if ~isfile(inputFile)
        error('Arquivo de entrada não encontrado: %s', inputFile);
    end

    % Ler arquivo tabulado, sem cabeçalho
    T = readtable(inputFile, ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'ReadVariableNames', false, ...
        'VariableNamingRule', 'preserve');

    % Primeira coluna = nomes/IDs dos genes
    gene_ids = T{:,1};

    % Garantir formato string
    if iscell(gene_ids)
        gene_ids = string(gene_ids);
    elseif isnumeric(gene_ids)
        gene_ids = string(gene_ids);
    elseif ischar(gene_ids)
        gene_ids = string(cellstr(gene_ids));
    end

    % Colunas 2:end = matriz numérica
    A = T{:,2:end};

    fprintf('Matriz A carregada do arquivo com sucesso.\n');
    fprintf('Número de genes lidos: %d\n', size(A,1));
    fprintf('Número de amostras lidas: %d\n', size(A,2));

else
    fprintf('A variável A já existe no workspace. Usando matriz já carregada.\n');

    if ~exist('gene_ids', 'var')
        warning(['A matriz A já existe, mas gene_ids não existe no workspace. ' ...
                 'O script poderá identificar os índices dos genes, ' ...
                 'mas não mostrará os nomes/IDs reais.']);
    end
end

%% ---------------------------------------------------------
% 1) Conferir dimensão da matriz
%% ---------------------------------------------------------
[m, n] = size(A);
fprintf('Dimensão original da matriz A: %d genes x %d amostras\n', m, n);

if n ~= 35
    warning(['Foram encontradas %d amostras. O desenho experimental do GSE38531 ' ...
             'descreve 35 amostras. Revise a matriz antes de seguir.'], n);
end

%% ---------------------------------------------------------
% 2) Remover linhas com NaN (sem imputação)
%% ---------------------------------------------------------
linhasComNaN = any(isnan(A), 2);
numNaN = sum(linhasComNaN);

fprintf('Número de genes com pelo menos um NaN: %d\n', numNaN);

A = A(~linhasComNaN, :);

% Filtrar gene_ids junto com A
if exist('gene_ids', 'var')
    gene_ids = gene_ids(~linhasComNaN);
end

[m, n] = size(A);
fprintf('Dimensão após remover genes com NaN: %d genes x %d amostras\n', m, n);

%% ---------------------------------------------------------
% 3) Definição dos grupos experimentais
%
% Ordem REAL das amostras:
% rep1: 1:7
% rep2: 8:14
% rep3: 15:21
% rep4: 22:28
% rep5: 29:35
%
% Dentro de cada réplica:
% 1 = T0 controle
% 2 = T2h infectado
% 3 = T0 + linezolida
% 4 = T0 + vancomicina
% 5 = T24h infectado
% 6 = T24h infectado + linezolida
% 7 = T24h infectado + vancomicina
%% ---------------------------------------------------------
if n ~= 35
    error(['Este script está configurado para 35 amostras. ' ...
           'Como sua matriz atual não tem 35 colunas, primeiro revise a montagem da A.']);
end

idx_G1 = [1 8 15 22 29];   % controle não infectado (T0)
idx_G2 = [3 10 17 24 31];  % não infectado + linezolida (T0)
idx_G3 = [4 11 18 25 32];  % não infectado + vancomicina (T0)
idx_G4 = [2 9 16 23 30];   % infectado 2h sem tratamento
idx_G5 = [5 12 19 26 33];  % infectado 24h sem tratamento
idx_G6 = [6 13 20 27 34];  % infectado 24h + linezolida
idx_G7 = [7 14 21 28 35];  % infectado 24h + vancomicina

Indicadores = zeros(7, n);
Indicadores(1, idx_G1) = 1;
Indicadores(2, idx_G2) = 1;
Indicadores(3, idx_G3) = 1;
Indicadores(4, idx_G4) = 1;
Indicadores(5, idx_G5) = 1;
Indicadores(6, idx_G6) = 1;
Indicadores(7, idx_G7) = 1;

if any(sum(Indicadores,1) ~= 1)
    error('Há amostras sem grupo definido ou associadas a mais de um grupo.');
end

%% ---------------------------------------------------------
% 4) Escolher as duas classes do classificador binário
%
% Classe 1 = infectado 24 h sem tratamento
% Classe 0 = controle não infectado
%% ---------------------------------------------------------
idxClasse1 = idx_G5;  % casos
idxClasse0 = idx_G1;  % controles

idxSel = [idxClasse1 idxClasse0];
A2 = A(:, idxSel);

nSel = size(A2, 2);
fprintf('Número de amostras usadas no classificador binário: %d\n', nSel);

iClasse1 = 1:length(idxClasse1);
iClasse0 = length(idxClasse1) + (1:length(idxClasse0));

%% ---------------------------------------------------------
% 5) Montar vetor b usando transformação logit
%
% Evitar 0 e 1 exatos para não gerar infinito no logit
%% ---------------------------------------------------------
p1 = 0.999999;
p0 = 0.000001;

lgch1 = log(p1 / (1 - p1));
lgch0 = log(p0 / (1 - p0));

b = zeros(nSel, 1);
b(iClasse1) = lgch1;
b(iClasse0) = lgch0;

%% ---------------------------------------------------------
% 6) Resolver para alpha no modelo completo
%
% Professor:
% alpha = resolve(M', b)
%
% Aqui M = A2
%% ---------------------------------------------------------
if exist('resolve', 'file') == 2
    alpha = resolve(A2', b);
    metodo_alpha = 'resolve(M'', b)';
else
    warning(['A função resolve() não foi encontrada. ' ...
             'Será usada uma solução regularizada como aproximação.']);
    lambda = 1e-6;
    alpha = A2 * ((A2' * A2 + lambda * eye(nSel)) \ b);
    metodo_alpha = 'aproximação regularizada';
end

fprintf('Método usado para alpha (modelo completo): %s\n', metodo_alpha);

%% ---------------------------------------------------------
% 7) Calcular p(x) no modelo completo
%% ---------------------------------------------------------
aux = A2' * alpha;
num = exp(aux);
p = num ./ (1 + num);

figure;
plot(1:nSel, p, '*');
grid on;
ylim([-0.05 1.05]);
title('Probabilidades p(x) - modelo completo');
xlabel('Amostras');
ylabel('p(x)');

disp('Probabilidades do modelo completo:');
disp(p);

%% ---------------------------------------------------------
% 8a) Seleção de atributos
%
% 10 alphas mais negativos + 10 mais positivos
%% ---------------------------------------------------------
[val, pos] = sort(alpha);

nAtributos = 10;

indxMin = pos(1:nAtributos);
indxMax = pos(end-nAtributos+1:end);

fprintf('\nAtributos mais negativos:\n');
disp(indxMin');

fprintf('Atributos mais positivos:\n');
disp(indxMax');

figure;
plot(alpha, '*');
hold on;
plot(indxMin, alpha(indxMin), 'or');
plot(indxMax, alpha(indxMax), 'og');
grid on;
title('Pesos alpha e atributos selecionados');
xlabel('Índice do atributo (gene)');
ylabel('alpha');
legend('Todos os alpha', 'Mais negativos', 'Mais positivos');
hold off;

%% ---------------------------------------------------------
% 8b) Guardar os 10 alphas mais negativos e os 10 mais positivos
%% ---------------------------------------------------------
alphaNegValores = alpha(indxMin);
alphaPosValores = alpha(indxMax);

if exist('gene_ids', 'var')
    top10Neg = table(indxMin(:), gene_ids(indxMin), alphaNegValores(:), ...
        'VariableNames', {'IndiceGene', 'GeneID', 'Alpha'});

    top10Pos = table(indxMax(:), gene_ids(indxMax), alphaPosValores(:), ...
        'VariableNames', {'IndiceGene', 'GeneID', 'Alpha'});
else
    top10Neg = table(indxMin(:), alphaNegValores(:), ...
        'VariableNames', {'IndiceGene', 'Alpha'});

    top10Pos = table(indxMax(:), alphaPosValores(:), ...
        'VariableNames', {'IndiceGene', 'Alpha'});
end

disp('Top 10 genes com alpha mais negativo:');
disp(top10Neg);

disp('Top 10 genes com alpha mais positivo:');
disp(top10Pos);

%% ---------------------------------------------------------
% 9) Montar matriz reduzida
%% ---------------------------------------------------------
Ar = [A2(indxMin, :); A2(indxMax, :)];

fprintf('Dimensão da matriz reduzida Ar: %d atributos x %d amostras\n', size(Ar,1), size(Ar,2));

%% ---------------------------------------------------------
% 10) Resolver alpha do modelo reduzido
%% ---------------------------------------------------------
alphar = Ar' \ b;

%% ---------------------------------------------------------
% 11) Calcular p(x) no modelo reduzido
%% ---------------------------------------------------------
auxr = Ar' * alphar;
numr = exp(auxr);
pr = numr ./ (1 + numr);

figure;
plot(1:nSel, pr, '*');
grid on;
ylim([-0.05 1.05]);
title('Probabilidades p(x) - modelo reduzido');
xlabel('Amostras');
ylabel('pr(x)');

disp('Probabilidades do modelo reduzido:');
disp(pr);

%% ---------------------------------------------------------
% 12) Verificar acurácia do classificador
%% ---------------------------------------------------------
y_true = [ones(5,1); zeros(5,1)];

y_pred_full = p >= 0.5;
acc_full = mean(y_pred_full == y_true);

y_pred_red = pr >= 0.5;
acc_red = mean(y_pred_red == y_true);

fprintf('Acurácia modelo completo: %.2f%%\n', acc_full*100);
fprintf('Acurácia modelo reduzido: %.2f%%\n', acc_red*100);

%% ---------------------------------------------------------
% 13) Salvar resultados
%% ---------------------------------------------------------
if exist('gene_ids', 'var')
    save('resultado_GSE38531_classificador.mat', ...
        'alpha', 'alphar', 'p', 'pr', ...
        'indxMin', 'indxMax', ...
        'alphaNegValores', 'alphaPosValores', ...
        'top10Neg', 'top10Pos', ...
        'idxClasse1', 'idxClasse0', ...
        'A2', 'Ar', 'Indicadores', 'b', ...
        'gene_ids');
else
    save('resultado_GSE38531_classificador.mat', ...
        'alpha', 'alphar', 'p', 'pr', ...
        'indxMin', 'indxMax', ...
        'alphaNegValores', 'alphaPosValores', ...
        'top10Neg', 'top10Pos', ...
        'idxClasse1', 'idxClasse0', ...
        'A2', 'Ar', 'Indicadores', 'b');
end

fprintf('\nResultados salvos em: resultado_GSE38531_classificador.mat\n');

%% ---------------------------------------------------------
% 14) Decomposição SVD e gráfico dos valores singulares
%% ---------------------------------------------------------
Msvd = A2;

[T, S, D] = svd(Msvd, 'econ');

A_reconstruida = T * S * D';
erro_reconstrucao = norm(Msvd - A_reconstruida, 'fro');
fprintf('Erro de reconstrução pela SVD: %.6e\n', erro_reconstrucao);

s = diag(S);

figure;
plot(1:length(s), s, '-o');
grid on;
title('Valores singulares de A2');
xlabel('Índice');
ylabel('Valor singular');

s_rel = s / sum(s);

figure;
plot(1:length(s_rel), s_rel, '-o');
grid on;
title('Valores singulares relativos de A2');
xlabel('Índice');
ylabel('Peso relativo');

disp('Valores singulares absolutos:');
disp(s);

disp('Valores singulares relativos:');
disp(s_rel);

s_acum = cumsum(s_rel);

figure;
plot(1:length(s_acum), s_acum, '-o');
grid on;
title('Soma acumulada dos valores singulares relativos');
xlabel('Número de componentes');
ylabel('Proporção acumulada');
ylim([0 1.05]);

disp('Soma acumulada dos valores singulares relativos:');
disp(s_acum);

%% ---------------------------------------------------------
% 15) Projeção nos 3 primeiros componentes principais
%% ---------------------------------------------------------
Aux = S * D';

x = Aux(1,:);
y = Aux(2,:);
z = Aux(3,:);

figure;
plot3(x, y, z, 'o');
grid on;

title('Projeção nos 3 primeiros componentes principais (SVD)');
xlabel('Componente 1');
ylabel('Componente 2');
zlabel('Componente 3');
