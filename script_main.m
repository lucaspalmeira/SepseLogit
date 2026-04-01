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
% 1) Antes de rodar, a variável A deve estar carregada no workspace
% 2) Confirme se A possui 35 colunas e se a ordem das amostras está correta
%% =========================================================

clearvars -except A
clc
close all

%% ---------------------------------------------------------
% 0) Verificar se A existe
%% ---------------------------------------------------------
if ~exist('A', 'var')
    error(['A variável A não está no workspace. ' ...
           'Carregue primeiro a matriz de expressão gênica no MATLAB.']);
end

%% ---------------------------------------------------------
% 1) Conferir dimensão da matriz
%% ---------------------------------------------------------
[m, n] = size(A);
fprintf('Dimensão original da matriz A: %d genes x %d amostras\n', m, n);

if n ~= 35
    warning(['Foram encontradas %d amostras. O desenho experimental do GSE38531 ' ...
             'descreve 35 amostras (7 grupos de 5). Revise a matriz antes de seguir.'], n);
end

%% ---------------------------------------------------------
% 2) Remover linhas com NaN (sem imputação)
%% ---------------------------------------------------------
linhasComNaN = any(isnan(A), 2);
numNaN = sum(linhasComNaN);

fprintf('Número de genes com pelo menos um NaN: %d\n', numNaN);

A = A(~linhasComNaN, :);

[m, n] = size(A);
fprintf('Dimensão após remover genes com NaN: %d genes x %d amostras\n', m, n);

%% ---------------------------------------------------------
% 3) Definição dos grupos experimentais
%
% Ordem assumida:
% G1 = controle não infectado
% G2 = não infectado + linezolida
% G3 = não infectado + vancomicina
% G4 = infectado 2h sem tratamento
% G5 = infectado 24h sem tratamento
% G6 = infectado 24h + linezolida
% G7 = infectado 24h + vancomicina
%
% Cada grupo com 5 amostras
%% ---------------------------------------------------------
if n ~= 35
    error(['Este script está configurado para 35 amostras distribuídas em 7 grupos de 5. ' ...
           'Como sua matriz atual não tem 35 colunas, primeiro revise a montagem da A.']);
end

idx_G1 = 1:5;
idx_G2 = 6:10;
idx_G3 = 11:15;
idx_G4 = 16:20;
idx_G5 = 21:25;
idx_G6 = 26:30;
idx_G7 = 31:35;

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
% 8) Seleção de atributos
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
% 9) Montar matriz reduzida
%% ---------------------------------------------------------
Ar = [A2(indxMin, :); A2(indxMax, :)];

fprintf('Dimensão da matriz reduzida Ar: %d atributos x %d amostras\n', size(Ar,1), size(Ar,2));

%% ---------------------------------------------------------
% 10) Resolver alpha do modelo reduzido
%
% Conforme o espírito do enunciado, agora temos poucos atributos.
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
save('resultado_GSE38531_classificador.mat', ...
    'alpha', 'alphar', 'p', 'pr', ...
    'indxMin', 'indxMax', ...
    'idxClasse1', 'idxClasse0', ...
    'A2', 'Ar', 'Indicadores', 'b');

fprintf('\nResultados salvos em: resultado_GSE38531_classificador.mat\n');
