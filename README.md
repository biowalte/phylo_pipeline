Projeto de Pipeline Filogenético com Docker
Este projeto utiliza um script Python (run_phylogeny.py) para orquestrar uma série de Biocontainers (Docker) e construir árvores filogenéticas a partir de um arquivo multi-fasta.

O objetivo é garantir uma análise 100% reprodutível, onde qualquer pessoa com Docker e Python possa rodar a mesma análise e obter os mesmos resultados, sem se preocupar com a instalação de dependências de bioinformática.

Ferramentas Utilizadas
Orquestrador: Python 3 (subprocess, pathlib)

Dependências do Host: Docker, Python 3

Biocontainers:

quay.io/biocontainers/mafft: Para alinhamento de sequências.

quay.io/biocontainers/trimal: Para filtragem/trimagem do alinhamento.

quay.io/biocontainers/biopython: Para conversão de formatos e análise de Neighbor-Joining.

quay.io/biocontainers/iqtree: Para análise de Máxima Verossimilhança.

quay.io/biocontainers/mrbayes: Para análise de Inferência Bayesiana.

Como Usar
Certifique-se de que você tem o Docker e o Python 3 instalados em seu sistema.

Coloque seu arquivo multi-fasta (ex: sequence.fasta) neste diretório.

Execute o script orquestrador:

    Bash

    python3 run_phylogeny.py
    O script pedirá o nome do seu arquivo de entrada (ex: sequence.fasta).

O script executará o alinhamento e a trimagem automaticamente.

O script pedirá que você escolha o método de construção da árvore (NJ, ML ou BI).

O Fluxo do Pipeline (Parâmetros)
O script executa a seguinte cadeia de comandos:

1. Input
O script solicita o nome do arquivo de entrada.

Input: sequence.fasta

2. Passo 1: Alinhamento (MAFFT)
O script chama o container mafft para alinhar as sequências.

Comando executado: mafft --auto /data/sequence.fasta > /data/sequence.aligned.fasta

Parâmetro: --auto seleciona automaticamente a melhor estratégia de alinhamento (L-INS-i, G-INS-i, etc.) com base no tamanho dos dados.

Resultado: sequence.aligned.fasta

3. Passo 2: Trimagem (trimAl)
O script chama o container trimal para remover posições do alinhamento que contêm muitos gaps ou são pouco informativas.

Comando executado: trimal -in ... -out ... -automated1

Parâmetro: -automated1 utiliza uma heurística otimizada para filogenia, geralmente mais "agressiva" na limpeza do que outros métodos.

Resultado: sequence.trimmed.fasta

4. Passo 3: Escolha do Método
Opção 1: Neighbor-Joining (NJ)
Ferramenta: biopython

Parâmetros: O script calcula uma matriz de distância (DistanceCalculator('identity')) e usa o algoritmo Neighbor-Joining (DistanceTreeConstructor('nj')) para construir a árvore.

Resultado: sequence.nj_tree.newick

Opção 2: Máxima Verossimilhança (ML)
Ferramenta: iqtree

Parâmetros:

-s: O alinhamento trimado (sequence.trimmed.fasta).

-m MFP: "ModelFinder Plus". O IQ-TREE testa automaticamente centenas de modelos de substituição e escolhe o melhor (ex: GTR+G+I) com base no BIC (Critério de Informação Bayesiano).

-B 1000: "Ultrafast Bootstrap". Executa 1000 réplicas de bootstrap para calcular os valores de suporte (confiança) para cada ramo da árvore.

Resultado: sequence.trimmed.fasta.treefile (e outros arquivos de log).

Opção 3: Inferência Bayesiana (BI)
Ferramenta: biopython (para conversão) e mrbayes (para análise).

Etapa A (Conversão): O script converte o alinhamento de FASTA para NEXUS.

Etapa B (Análise Interativa): O script inicia o MrBayes em modo interativo.

[!WARNING] Atenção: Ação Manual Requerida

Ao escolher a Inferência Bayesiana, o script não roda a análise sozinho. Ele apenas te coloca dentro do prompt MrBayes >.

Você deve digitar os comandos manualmente para iniciar a análise. Os comandos sugeridos pelo script (e que você pode alterar) são:

execute /data/sequence.trimmed.nex

(Carrega seu arquivo de alinhamento)

mcmc ngen=100000 samplefreq=100

(Inicia a análise MCMC com 100.000 gerações, salvando uma árvore a cada 100 gerações. Para uma publicação, este número ngen geralmente precisa ser muito maior, na casa dos milhões, até que a convergência seja atingida).

sumt

(Quando a análise terminar, este comando "summarize trees" descarta o "burn-in" inicial e calcula a árvore de consenso com as probabilidades posteriores).

quit

(Sai do MrBayes).

Resultado: sequence.trimmed.nex.con.tre (a árvore de consenso Bayesiana).