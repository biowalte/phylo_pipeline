# phylo_pipeline
Um pipeline filogenético automatizado que usa Python para orquestrar Biocontainers (Docker). Garante 100% de reprodutibilidade. Recebe um multi-fasta, alinha com MAFFT, filtra com trimAl e constrói árvores usando Neighbor-Joining (Biopython), Máxima Verossimilhança (IQ-TREE) ou Inferência Bayesiana (MrBayes)
