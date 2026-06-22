%% =========================================================
% GSE38531 - Analise multigrupo para candidatos a biomarcadores
% Repositorio: SepseLogit
% Local sugerido: scripts/script_biomarcadores_multigrupos.m
%
% Este script foi adaptado para a estrutura do repositorio:
% .
% ├── data/
% │   └── gse_data_no_header.txt
% ├── scripts/
% │   ├── script_main.m
% │   ├── resolve.m
% │   └── script_biomarcadores_multigrupos.m
% ├── figures/
% │   └── biomarkers_multigrupo/
% └── results/
%     └── biomarkers_multigrupo/
%
% Objetivos:
% 1) Reproduzir a logica do script_main.m: classificacao binaria,
%    coeficientes alpha, selecao dos 10 pesos negativos e 10 positivos,
%    modelo reduzido, PCA/SVD, kNN e SVM.
% 2) Aplicar essa logica a varios contrastes biologicos do GSE38531.
% 3) Gerar tabelas de genes/probes candidatos a biomarcadores usando:
%    - coeficiente alpha;
%    - diferenca media de expressao;
%    - teste t de Welch;
%    - FDR Benjamini-Hochberg;
%    - AUC univariada;
%    - escore integrado de evidencia;
%    - frequencia do probe entre contrastes.
%
% Observacao importante:
% Este script gera candidatos a biomarcadores, mas nao valida biomarcadores.
% O dataset possui apenas 5 amostras por grupo. Para artigo, recomenda-se
% validacao externa, validacao experimental e comparacao com metodos como
% limma, elastic net, random forest e SVM-RFE.
%% =========================================================

clearvars
clc
close all

%% ---------------------------------------------------------
% 0) Definir diretorios do projeto seguindo o repositorio
%% ---------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

% O script deve ficar em SepseLogit/scripts/.
% Assim, o diretorio raiz do projeto e a pasta acima de scripts/.
projectRoot = fileparts(scriptDir);

% Caso o usuario execute uma copia do script a partir da raiz do projeto,
% ajusta automaticamente.
[~, currentFolderName] = fileparts(scriptDir);
if ~strcmp(currentFolderName, 'scripts') && isfolder(fullfile(scriptDir, 'data'))
    projectRoot = scriptDir;
    scriptDir = fullfile(projectRoot, 'scripts');
end

% Garante que scripts/ esteja no path, para encontrar resolve.m.
if isfolder(scriptDir)
    addpath(scriptDir);
end

inputFile = fullfile(projectRoot, 'data', 'gse_data_no_header.txt');
figRoot   = fullfile(projectRoot, 'figures');
figDir    = fullfile(figRoot, 'biomarkers_multigrupo');
outRoot   = fullfile(projectRoot, 'results');
outDir    = fullfile(outRoot, 'biomarkers_multigrupo');

if ~isfile(inputFile)
    error(['Arquivo de entrada nao encontrado: %s\n', ...
           'Coloque gse_data_no_header.txt em: %s'], inputFile, fullfile(projectRoot, 'data'));
end

if ~exist(figRoot, 'dir'); mkdir(figRoot); end
if ~exist(figDir,  'dir'); mkdir(figDir);  end
if ~exist(outRoot, 'dir'); mkdir(outRoot); end
if ~exist(outDir,  'dir'); mkdir(outDir);  end

fprintf('Diretorio do projeto: %s\n', projectRoot);
fprintf('Arquivo de entrada: %s\n', inputFile);
fprintf('Figuras serao salvas em: %s\n', figDir);
fprintf('Resultados serao salvos em: %s\n', outDir);

% Parametros principais
nTopAlpha = 10;              % 10 negativos + 10 positivos, como no script original
nTopBiomarkers = 50;         % top candidatos integrados por contraste
nTopVariableGenesPCA = 2000; % PCA global usando genes mais variaveis

%% ---------------------------------------------------------
% 1) Carregar matriz de expressao
%% ---------------------------------------------------------
T = readtable(inputFile, ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'ReadVariableNames', false, ...
    'VariableNamingRule', 'preserve');

gene_ids = T{:,1};
if iscell(gene_ids)
    gene_ids = string(gene_ids);
elseif isnumeric(gene_ids)
    gene_ids = string(gene_ids);
elseif ischar(gene_ids)
    gene_ids = string(cellstr(gene_ids));
else
    gene_ids = string(gene_ids);
end

A = double(T{:,2:end});
clear T

[m, n] = size(A);
fprintf('Dimensao original da matriz A: %d probes/genes x %d amostras\n', m, n);

if n ~= 35
    error('Este script espera uma matriz com 35 amostras. Verifique o arquivo data/gse_data_no_header.txt.');
end

%% ---------------------------------------------------------
% 2) Remover genes com NaN e genes sem variacao
%% ---------------------------------------------------------
linhasComNaN = any(isnan(A), 2);
linhasSemVariacao = std(A, 0, 2) == 0;
remover = linhasComNaN | linhasSemVariacao;

fprintf('Genes com pelo menos um NaN: %d\n', sum(linhasComNaN));
fprintf('Genes sem variacao: %d\n', sum(linhasSemVariacao));

A = A(~remover, :);
gene_ids = gene_ids(~remover);
[m, n] = size(A);
fprintf('Dimensao apos filtros: %d probes/genes x %d amostras\n', m, n);

%% ---------------------------------------------------------
% 3) Definir os 7 grupos experimentais
%% ---------------------------------------------------------
% Os indices abaixo seguem o desenho usado no script_main.m.
groupNames = [ ...
    "G1_Control_T0", ...
    "G2_Control_Linezolid", ...
    "G3_Control_Vancomycin", ...
    "G4_Infected_2h", ...
    "G5_Infected_24h_Untreated", ...
    "G6_Infected_24h_Linezolid", ...
    "G7_Infected_24h_Vancomycin"];

groupDesc = [ ...
    "controle nao infectado T0", ...
    "nao infectado + linezolida", ...
    "nao infectado + vancomicina", ...
    "infectado 2h", ...
    "infectado 24h sem tratamento", ...
    "infectado 24h + linezolida", ...
    "infectado 24h + vancomicina"];

groupIdx = { ...
    [1 8 15 22 29], ...
    [3 10 17 24 31], ...
    [4 11 18 25 32], ...
    [2 9 16 23 30], ...
    [5 12 19 26 33], ...
    [6 13 20 27 34], ...
    [7 14 21 28 35]};

Indicadores = zeros(7, n);
for g = 1:7
    Indicadores(g, groupIdx{g}) = 1;
end
if any(sum(Indicadores,1) ~= 1)
    error('Ha amostras sem grupo definido ou com mais de um grupo.');
end

indicesAmostras = string(cellfun(@mat2str, groupIdx, 'UniformOutput', false));
indicesAmostras = indicesAmostras(:);
groupInfo = table((1:7)', groupNames(:), groupDesc(:), indicesAmostras, ...
    'VariableNames', {'Grupo', 'NomeCurto', 'Descricao', 'IndicesAmostras'});
writetable(groupInfo, fullfile(outDir, '00_grupos_experimentais.csv'));

%% ---------------------------------------------------------
% 4) PCA/SVD global com todos os grupos
%% ---------------------------------------------------------
% Para visualizacao global, usa genes mais variaveis para reduzir ruido.
geneVar = var(A, 0, 2);
[~, idxVar] = sort(geneVar, 'descend');
idxPCA = idxVar(1:min(nTopVariableGenesPCA, length(idxVar)));
A_pca = A(idxPCA, :);
A_pca_centered = A_pca - mean(A_pca, 2);
[~, Sglobal, Dglobal] = svd(A_pca_centered, 'econ');
scoresGlobal = Sglobal * Dglobal';

sampleGroup = zeros(n,1);
for g = 1:7
    sampleGroup(groupIdx{g}) = g;
end

% Garante que todas as variaveis da tabela tenham formato coluna (n x 1).
sampleID = (1:n)';
sampleGroup = sampleGroup(:);
sampleGroupName = reshape(groupNames(sampleGroup), [], 1);
PC1 = scoresGlobal(1,:)';
PC2 = scoresGlobal(2,:)';
PC3 = scoresGlobal(3,:)';

pcaTable = table(sampleID, sampleGroup, sampleGroupName, PC1, PC2, PC3, ...
    'VariableNames', {'Amostra', 'GrupoID', 'Grupo', 'PC1', 'PC2', 'PC3'});
writetable(pcaTable, fullfile(outDir, '01_PCA_global_scores.csv'));

fig = figure('Visible','off');
hold on
for g = 1:7
    idx = sampleGroup == g;
    scatter(scoresGlobal(1,idx), scoresGlobal(2,idx), 70, 'filled');
end
grid on
xlabel('PC1'); ylabel('PC2');
title('PCA/SVD global - todos os grupos');
legend(groupNames, 'Interpreter', 'none', 'Location', 'bestoutside');
hold off
save_figure(fig, fullfile(figDir, '01_PCA_global_PC1_PC2.png'));
close(fig)

%% ---------------------------------------------------------
% 5) ANOVA univariada entre os 7 grupos
%% ---------------------------------------------------------
[anovaF, anovaP, anovaFDR] = one_way_anova_matrix(A, groupIdx);
anovaTable = table((1:m)', gene_ids(:), anovaF(:), anovaP(:), anovaFDR(:), ...
    'VariableNames', {'IndiceGene', 'ProbeSetID', 'F_ANOVA_7grupos', 'P_ANOVA_7grupos', 'FDR_ANOVA_7grupos'});
anovaTable = sortrows(anovaTable, {'FDR_ANOVA_7grupos', 'P_ANOVA_7grupos', 'F_ANOVA_7grupos'}, {'ascend', 'ascend', 'descend'});
writetable(anovaTable, fullfile(outDir, '02_ANOVA_7_grupos_todos_genes.csv'));

%% ---------------------------------------------------------
% 6) Contrastes biologicos principais
%% ---------------------------------------------------------
% Em cada contraste, a classe positiva e a classe interpretada como "caso".
% Alpha positivo = associado a classe positiva.
% Alpha negativo = associado a classe negativa.
% Inicializa a struct diretamente com o primeiro contraste.
% Isso evita o erro "Subscripted assignment between dissimilar structures"
% em algumas versoes do MATLAB ao fazer contrasts = struct([]) e depois
% contrasts(end+1) = make_contrast(...).
contrasts = make_contrast('G4_vs_G1_infeccao_2h',               4, 1, groupNames, groupIdx, 'infectado 2h vs controle T0');
contrasts(end+1) = make_contrast('G5_vs_G1_infeccao_24h',              5, 1, groupNames, groupIdx, 'infectado 24h sem tratamento vs controle T0');
contrasts(end+1) = make_contrast('G5_vs_G4_progressao_2h_para_24h',    5, 4, groupNames, groupIdx, 'progressao da infeccao: 24h vs 2h');
contrasts(end+1) = make_contrast('G6_vs_G5_efeito_linezolida',         6, 5, groupNames, groupIdx, 'efeito da linezolida em animais infectados 24h');
contrasts(end+1) = make_contrast('G7_vs_G5_efeito_vancomicina',        7, 5, groupNames, groupIdx, 'efeito da vancomicina em animais infectados 24h');
contrasts(end+1) = make_contrast('G6_vs_G7_linezolida_vs_vancomicina', 6, 7, groupNames, groupIdx, 'linezolida vs vancomicina em animais infectados 24h');
contrasts(end+1) = make_contrast('G2_vs_G1_linezolida_sem_infeccao',   2, 1, groupNames, groupIdx, 'linezolida em animal nao infectado vs controle T0');
contrasts(end+1) = make_contrast('G3_vs_G1_vancomicina_sem_infeccao',  3, 1, groupNames, groupIdx, 'vancomicina em animal nao infectado vs controle T0');

summaryRows = cell(numel(contrasts), 1);
topAll = table();

for c = 1:numel(contrasts)
    fprintf('\n=========================================================\n');
    fprintf('Contraste %d/%d: %s\n', c, numel(contrasts), contrasts(c).name);
    fprintf('%s\n', contrasts(c).description);

    result = run_binary_contrast(A, gene_ids, contrasts(c), nTopAlpha, nTopBiomarkers, outDir, figDir);
    summaryRows{c} = result.summary;

    topTmp = result.topBiomarkers;
    topTmp.Contraste = repmat(string(contrasts(c).name), height(topTmp), 1);
    topTmp.ClassePositiva = repmat(string(contrasts(c).positiveName), height(topTmp), 1);
    topTmp.ClasseNegativa = repmat(string(contrasts(c).negativeName), height(topTmp), 1);
    if isempty(topAll)
        topAll = topTmp;
    else
        topAll = [topAll; topTmp]; %#ok<AGROW>
    end
end

summaryTable = vertcat(summaryRows{:});
writetable(summaryTable, fullfile(outDir, '03_resumo_contrastes_classificacao.csv'));

%% ---------------------------------------------------------
% 7) Ranking consenso de candidatos a biomarcadores
%% ---------------------------------------------------------
if ~isempty(topAll)
    writetable(topAll, fullfile(outDir, '04_top_biomarcadores_todos_contrastes.csv'));

    [uniqueProbe, ~, ic] = unique(topAll.ProbeSetID);
    freq = accumarray(ic, 1);
    meanScore = accumarray(ic, topAll.EvidenceScore, [], @mean);
    maxScore = accumarray(ic, topAll.EvidenceScore, [], @max);
    maxAbsAlpha = accumarray(ic, abs(topAll.Alpha), [], @max);
    maxAbsDelta = accumarray(ic, abs(topAll.DeltaMean_Log2), [], @max);
    minFDR = accumarray(ic, topAll.FDR_Welch, [], @min);
    maxAUCsep = accumarray(ic, topAll.AUC_Separation, [], @max);

    consensus = table(uniqueProbe, freq, meanScore, maxScore, maxAbsAlpha, maxAbsDelta, minFDR, maxAUCsep, ...
        'VariableNames', {'ProbeSetID', 'FrequenciaTop50', 'MeanEvidenceScore', 'MaxEvidenceScore', ...
        'MaxAbsAlpha', 'MaxAbsDeltaMean_Log2', 'MinFDR_Welch', 'MaxAUC_Separation'});
    consensus = sortrows(consensus, {'FrequenciaTop50', 'MaxEvidenceScore', 'MaxAUC_Separation'}, {'descend', 'descend', 'descend'});
    writetable(consensus, fullfile(outDir, '05_ranking_consenso_biomarcadores.csv'));

    % Heatmap dos 50 candidatos mais recorrentes/fortes
    nHeat = min(50, height(consensus));
    probesHeat = consensus.ProbeSetID(1:nHeat);
    [~, idxHeat] = ismember(probesHeat, gene_ids);
    idxHeat = idxHeat(idxHeat > 0);
    if ~isempty(idxHeat)
        Xheat = zscore_rows(A(idxHeat, :));
        fig = figure('Visible','off');
        imagesc(Xheat);
        colorbar;
        title('Top candidatos consenso - z-score por probe');
        xlabel('Amostras'); ylabel('Probes');
        yticks(1:numel(idxHeat));
        yticklabels(gene_ids(idxHeat));
        set(gca, 'TickLabelInterpreter', 'none');
        save_figure(fig, fullfile(figDir, '02_heatmap_top_consenso.png'));
        close(fig)
    end
end

%% ---------------------------------------------------------
% 8) Salvar workspace reduzido
%% ---------------------------------------------------------
save(fullfile(outDir, 'workspace_GSE38531_biomarcadores.mat'), ...
    'gene_ids', 'groupNames', 'groupDesc', 'groupIdx', 'Indicadores', ...
    'summaryTable', 'anovaTable', 'topAll');

fprintf('\nAnalise concluida.\n');
fprintf('Tabelas: %s\n', outDir);
fprintf('Figuras: %s\n', figDir);

%% =========================================================
% FUNCOES LOCAIS
%% =========================================================

function C = make_contrast(name, posGroup, negGroup, groupNames, groupIdx, description)
    C.name = string(name);
    C.positiveGroup = posGroup;
    C.negativeGroup = negGroup;
    C.positiveName = groupNames(posGroup);
    C.negativeName = groupNames(negGroup);
    C.positiveIdx = groupIdx{posGroup};
    C.negativeIdx = groupIdx{negGroup};
    C.description = string(description);
end

function result = run_binary_contrast(A, gene_ids, C, nTopAlpha, nTopBiomarkers, outDir, figDir)
    safeName = sanitize_filename(C.name);
    contrastDir = fullfile(outDir, char(safeName));
    if ~exist(contrastDir, 'dir'); mkdir(contrastDir); end

    idxSel = [C.positiveIdx(:); C.negativeIdx(:)]';
    A2 = A(:, idxSel);
    y = [ones(numel(C.positiveIdx),1); zeros(numel(C.negativeIdx),1)];
    nSel = numel(y);

    b = labels_to_logit(y);

    % Modelo completo: alpha para todos os genes/probes.
    alpha = solve_alpha_modified_logistic(A2', b);
    probFull = sigmoid(A2' * alpha);
    predFull = probFull >= 0.5;
    accFull = mean(predFull == y);

    % Selecionar top negativos e positivos por alpha.
    [~, posSort] = sort(alpha, 'ascend');
    idxNeg = posSort(1:nTopAlpha);
    idxPos = flipud(posSort(end-nTopAlpha+1:end));
    idxReduced = [idxNeg(:); idxPos(:)];

    topNeg = build_alpha_table(idxNeg, gene_ids, alpha, C.negativeName, C.positiveName, "negativo");
    topPos = build_alpha_table(idxPos, gene_ids, alpha, C.negativeName, C.positiveName, "positivo");
    writetable(topNeg, fullfile(contrastDir, 'top10_alpha_negativo.csv'));
    writetable(topPos, fullfile(contrastDir, 'top10_alpha_positivo.csv'));

    % Modelo reduzido com os 20 genes/probes selecionados.
    Ar = A2(idxReduced, :);
    alphar = Ar' \ b;
    probReduced = sigmoid(Ar' * alphar);
    predReduced = probReduced >= 0.5;
    accReduced = mean(predReduced == y);

    % Validacao leave-one-out com selecao feita dentro de cada fold.
    [accLOOCV_alpha, predLOOCV_alpha, probLOOCV_alpha] = loocv_alpha_signature(A2, y, nTopAlpha);

    % kNN e SVM no conjunto reduzido. Acuracia aparente e LOOCV.
    Xred = Ar';
    accKNN_app = NaN; accKNN_loocv = NaN;
    accSVM_app = NaN; accSVM_loocv = NaN;
    try
        mdlKNN = fitcknn(Xred, y, 'NumNeighbors', 3, 'Standardize', true);
        predKNN = predict(mdlKNN, Xred);
        accKNN_app = mean(predKNN == y);
        cvKNN = crossval(mdlKNN, 'Leaveout', 'on');
        accKNN_loocv = 1 - kfoldLoss(cvKNN);
    catch ME
        warning('kNN nao executado para %s: %s', C.name, ME.message);
    end
    try
        mdlSVM = fitcsvm(Xred, y, 'KernelFunction', 'linear', 'Standardize', true);
        predSVM = predict(mdlSVM, Xred);
        accSVM_app = mean(predSVM == y);
        cvSVM = crossval(mdlSVM, 'Leaveout', 'on');
        accSVM_loocv = 1 - kfoldLoss(cvSVM);
    catch ME
        warning('SVM nao executado para %s: %s', C.name, ME.message);
    end

    % Estatisticas univariadas por gene/probe.
    posData = A(:, C.positiveIdx);
    negData = A(:, C.negativeIdx);
    meanPos = mean(posData, 2);
    meanNeg = mean(negData, 2);
    deltaMean = meanPos - meanNeg; % dados de microarray normalmente em escala log2; aproximacao de log2FC.
    [pWelch, tWelch] = welch_ttest_matrix(posData, negData);
    fdrWelch = bh_fdr(pWelch);
    [aucPosHigher, aucSep] = auc_matrix(posData, negData);

    scoreAlpha = minmax01(abs(alpha));
    scoreDelta = minmax01(abs(deltaMean));
    scoreFDR = minmax01(-log10(fdrWelch + realmin));
    scoreAUC = minmax01(aucSep);
    evidenceScore = scoreAlpha + scoreDelta + scoreFDR + scoreAUC;

    allGenes = table((1:numel(gene_ids))', gene_ids(:), alpha(:), meanPos(:), meanNeg(:), deltaMean(:), ...
        tWelch(:), pWelch(:), fdrWelch(:), aucPosHigher(:), aucSep(:), evidenceScore(:), ...
        'VariableNames', {'IndiceGene', 'ProbeSetID', 'Alpha', 'MeanPositive', 'MeanNegative', ...
        'DeltaMean_Log2', 'T_Welch', 'P_Welch', 'FDR_Welch', 'AUC_PositiveHigher', ...
        'AUC_Separation', 'EvidenceScore'});

    allGenes = sortrows(allGenes, {'EvidenceScore', 'AUC_Separation', 'FDR_Welch'}, {'descend', 'descend', 'ascend'});
    writetable(allGenes, fullfile(contrastDir, 'todos_genes_ranqueados.csv'));

    topBiomarkers = allGenes(1:min(nTopBiomarkers, height(allGenes)), :);
    writetable(topBiomarkers, fullfile(contrastDir, 'top50_biomarcadores_integrado.csv'));

    % Tabela com probabilidades por amostra.
    probTable = table((1:nSel)', idxSel(:), y(:), probFull(:), predFull(:), probReduced(:), predReduced(:), ...
        probLOOCV_alpha(:), predLOOCV_alpha(:), ...
        'VariableNames', {'AmostraLocal', 'AmostraOriginal', 'ClasseReal', 'ProbModeloCompleto', ...
        'PredModeloCompleto', 'ProbModeloReduzido', 'PredModeloReduzido', 'ProbLOOCV_AlphaSignature', ...
        'PredLOOCV_AlphaSignature'});
    writetable(probTable, fullfile(contrastDir, 'probabilidades_amostras.csv'));

    % Figuras do contraste.
    fig = figure('Visible','off');
    plot(1:nSel, probFull, '*', 'MarkerSize', 8); hold on
    plot(1:nSel, probReduced, 'o', 'MarkerSize', 8);
    yline(0.5, '--');
    grid on; ylim([-0.05 1.05]);
    title(sprintf('Probabilidades - %s', C.name), 'Interpreter', 'none');
    xlabel('Amostras'); ylabel('Probabilidade da classe positiva');
    legend('Modelo completo', 'Modelo reduzido', 'Limiar 0.5', 'Location', 'best');
    save_figure(fig, fullfile(figDir, char(safeName + "_probabilidades.png")));
    close(fig)

    fig = figure('Visible','off');
    plot(alpha, '*'); hold on
    plot(idxNeg, alpha(idxNeg), 'or');
    plot(idxPos, alpha(idxPos), 'og');
    grid on
    title(sprintf('Coeficientes alpha - %s', C.name), 'Interpreter', 'none');
    xlabel('Indice do probe/gene'); ylabel('Alpha');
    legend('Todos', 'Top negativos', 'Top positivos', 'Location', 'best');
    save_figure(fig, fullfile(figDir, char(safeName + "_alpha.png")));
    close(fig)

    fig = figure('Visible','off');
    scatter(deltaMean, -log10(fdrWelch + realmin), 10, 'filled');
    grid on
    xlabel('Delta medio de expressao: positivo - negativo');
    ylabel('-log10(FDR)');
    title(sprintf('Volcano aproximado - %s', C.name), 'Interpreter', 'none');
    save_figure(fig, fullfile(figDir, char(safeName + "_volcano.png")));
    close(fig)

    % Resumo do contraste.
    result.summary = table(string(C.name), string(C.description), string(C.positiveName), string(C.negativeName), ...
        numel(C.positiveIdx), numel(C.negativeIdx), accFull, accReduced, accLOOCV_alpha, accKNN_app, accKNN_loocv, accSVM_app, accSVM_loocv, ...
        'VariableNames', {'Contraste', 'Descricao', 'ClassePositiva', 'ClasseNegativa', ...
        'N_Positive', 'N_Negative', 'AccFull_Apparent', 'AccReduced_Apparent', ...
        'AccAlphaSignature_LOOCV', 'AccKNN_Apparent', 'AccKNN_LOOCV', 'AccSVM_Apparent', 'AccSVM_LOOCV'});

    result.topBiomarkers = topBiomarkers;
    result.topNeg = topNeg;
    result.topPos = topPos;
end

function tbl = build_alpha_table(idx, gene_ids, alpha, negativeName, positiveName, direction)
    if direction == "positivo"
        assoc = repmat(string(positiveName), numel(idx), 1);
    else
        assoc = repmat(string(negativeName), numel(idx), 1);
    end
    tbl = table(idx(:), gene_ids(idx(:)), alpha(idx(:)), assoc, ...
        'VariableNames', {'IndiceGene', 'ProbeSetID', 'Alpha', 'ClasseAssociada'});
end

function b = labels_to_logit(y)
    p1 = 0.999999;
    p0 = 0.000001;
    lgch1 = log(p1 / (1 - p1));
    lgch0 = log(p0 / (1 - p0));
    b = zeros(numel(y), 1);
    b(y == 1) = lgch1;
    b(y == 0) = lgch0;
end

function alpha = solve_alpha_modified_logistic(A_samples_by_genes, b)
    % Usa resolve.m quando ele estiver disponivel no path; caso contrario,
    % executa internamente a mesma formulacao por sistema aumentado.
    if exist('resolve', 'file') == 2
        alpha = resolve(A_samples_by_genes, b);
        return
    end

    [m, n] = size(A_samples_by_genes);
    Im = speye(m);
    In = speye(n);
    As = sparse(A_samples_by_genes);
    M = [Im, -As; As', In];
    nb = zeros(m+n, 1);
    nb(1:m) = -b;
    x = M \ nb;
    alpha = x(m+1:end);
end

function p = sigmoid(x)
    % Funcao logistica numericamente estavel.
    p = zeros(size(x));
    idx = x >= 0;
    p(idx) = 1 ./ (1 + exp(-x(idx)));
    ex = exp(x(~idx));
    p(~idx) = ex ./ (1 + ex);
end

function [acc, pred, prob] = loocv_alpha_signature(A2, y, nTopAlpha)
    n = numel(y);
    pred = false(n,1);
    prob = NaN(n,1);
    for i = 1:n
        trainIdx = true(n,1);
        trainIdx(i) = false;
        yTrain = y(trainIdx);
        bTrain = labels_to_logit(yTrain);

        Atrain = A2(:, trainIdx);
        alphaTrain = solve_alpha_modified_logistic(Atrain', bTrain);
        [~, posSort] = sort(alphaTrain, 'ascend');
        idxUse = [posSort(1:nTopAlpha); flipud(posSort(end-nTopAlpha+1:end))];

        ArTrain = Atrain(idxUse, :);
        alpharTrain = ArTrain' \ bTrain;
        prob(i) = sigmoid(A2(idxUse, i)' * alpharTrain);
        pred(i) = prob(i) >= 0.5;
    end
    acc = mean(pred == y);
end

function [p, tstat] = welch_ttest_matrix(X1, X0)
    n1 = size(X1, 2);
    n0 = size(X0, 2);
    mu1 = mean(X1, 2);
    mu0 = mean(X0, 2);
    v1 = var(X1, 0, 2);
    v0 = var(X0, 0, 2);
    se = sqrt(v1/n1 + v0/n0);
    tstat = (mu1 - mu0) ./ se;
    df = (v1/n1 + v0/n0).^2 ./ ((v1.^2)/((n1^2)*(n1-1)) + (v0.^2)/((n0^2)*(n0-1)));
    invalid = se == 0 | isnan(se) | isnan(df) | df <= 0;
    tstat(invalid) = 0;
    df(invalid) = 1;

    % p-valor bicaudal com fallback sem Statistics Toolbox.
    p = 2 * t_upper_tail_local(abs(tstat), df);
    p(invalid) = 1;
    p(~isfinite(p)) = 1;
end

function tail = t_upper_tail_local(t, v)
    % P(T >= t) para t >= 0, T ~ t Student(v).
    x = v ./ (v + t.^2);
    tail = 0.5 * betainc(x, v/2, 0.5);
end

function q = bh_fdr(p)
    q = NaN(size(p));
    valid = ~isnan(p);
    pValid = p(valid);
    [pSort, order] = sort(pValid(:), 'ascend');
    m = numel(pSort);
    if m == 0
        return
    end
    qSort = pSort .* m ./ (1:m)';
    qSort = min(1, qSort);
    qSort = flipud(cummin(flipud(qSort)));
    qValid = NaN(size(pValid));
    qValid(order) = qSort;
    q(valid) = qValid;
end

function [aucPosHigher, aucSep] = auc_matrix(posData, negData)
    m = size(posData, 1);
    nPos = size(posData, 2);
    nNeg = size(negData, 2);
    aucPosHigher = NaN(m,1);
    for i = 1:m
        scores = [posData(i,:), negData(i,:)]';
        labels = [ones(nPos,1); zeros(nNeg,1)];
        r = tied_ranks_local(scores);
        rankPos = sum(r(labels == 1));
        aucPosHigher(i) = (rankPos - nPos*(nPos+1)/2) / (nPos*nNeg);
    end
    aucSep = max(aucPosHigher, 1 - aucPosHigher);
end

function r = tied_ranks_local(x)
    [xs, order] = sort(x(:));
    rSorted = zeros(size(xs));
    i = 1;
    n = numel(xs);
    while i <= n
        j = i;
        while j < n && xs(j+1) == xs(i)
            j = j + 1;
        end
        rSorted(i:j) = (i + j) / 2;
        i = j + 1;
    end
    r = zeros(size(x(:)));
    r(order) = rSorted;
end

function y = minmax01(x)
    x = x(:);
    finiteIdx = isfinite(x);
    y = zeros(size(x));
    if ~any(finiteIdx)
        return
    end
    xmin = min(x(finiteIdx));
    xmax = max(x(finiteIdx));
    if xmax == xmin
        y(finiteIdx) = 0;
    else
        y(finiteIdx) = (x(finiteIdx) - xmin) ./ (xmax - xmin);
    end
end

function [F, p, fdr] = one_way_anova_matrix(A, groupIdx)
    K = numel(groupIdx);
    N = size(A, 2);
    grandMean = mean(A, 2);
    SSB = zeros(size(A,1), 1);
    SSW = zeros(size(A,1), 1);
    for k = 1:K
        Xk = A(:, groupIdx{k});
        nk = size(Xk, 2);
        muk = mean(Xk, 2);
        SSB = SSB + nk * (muk - grandMean).^2;
        SSW = SSW + sum((Xk - muk).^2, 2);
    end
    dfB = K - 1;
    dfW = N - K;
    F = (SSB / dfB) ./ (SSW / dfW);
    F(~isfinite(F)) = 0;

    % p = P(F_dist >= F). Implementado via betainc para nao depender de fcdf.
    x = (dfB .* F) ./ (dfB .* F + dfW);
    cdfVal = betainc(x, dfB/2, dfW/2);
    p = 1 - cdfVal;
    p(~isfinite(p)) = 1;
    fdr = bh_fdr(p);
end

function z = zscore_rows(X)
    mu = mean(X, 2);
    sd = std(X, 0, 2);
    sd(sd == 0) = 1;
    z = (X - mu) ./ sd;
end

function s = sanitize_filename(name)
    s = regexprep(string(name), '[^a-zA-Z0-9_\-]', '_');
end

function save_figure(figHandle, filePath)
    try
        exportgraphics(figHandle, filePath, 'Resolution', 300);
    catch
        saveas(figHandle, filePath);
    end
end
