# Arquivo: run_phylogeny.py
# VERSÃO CORRIGIDA

import os
import subprocess
import sys
from pathlib import Path

# --- Configuração dos Nomes dos Containers ---
CONTAINERS = {
    "mafft": "quay.io/biocontainers/mafft:7.525--h031d066_0",
    "trimal": "quay.io/biocontainers/trimal:1.5",
    "iqtree": "quay.io/biocontainers/iqtree:3.0.1--h503566f_0",
    "biopython": "quay.io/biocontainers/biopython:1.79",
    "mrbayes": "quay.io/biocontainers/mrbayes:3.2.7--hd0d793b_7"
}

# --- Função Auxiliar para Rodar Comandos Docker ---
# MUDANÇA AQUI: 'command' agora é 'list'
def run_docker_command(container_name: str, command: list, interactive: bool = False):
    """
    Executa um comando dentro de um container Docker,
    montando o diretório de trabalho atual.
    """
    work_dir = Path.cwd().resolve()
    
    docker_base = [
        "docker", "run", "--rm",
        "-v", f"{work_dir}:/data",
        "-w", "/data"
    ]
    
    if interactive:
        docker_base.extend(["-it"])
        
    # MUDANÇA AQUI: 'command' já é uma lista, não precisa de .split()
    full_command = docker_base + [CONTAINERS[container_name]] + command

    # Usamos 'shlex.join' para imprimir o comando de forma legível
    # (import shlex no topo do arquivo seria ideal, mas ' '.join funciona para log)
    print(f"\n--- Executando: {' '.join(full_command)} ---")
    try:
        subprocess.run(full_command, text=True, check=True)
        print("--- Comando concluído com sucesso! ---")
    except subprocess.CalledProcessError as e:
        print(f"ERRO ao executar o comando: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nExecução interrompida pelo usuário.", file=sys.stderr)
        sys.exit(1)


# --- Função Principal do Pipeline ---
def main():
    print("Iniciando o Pipeline de Filogenia com Containers")
    
    # --- PASSO 0: INPUT ---
    try:
        input_file = input("Por favor, digite o nome do seu arquivo multi-fasta: ")
    except EOFError:
        print("\nEntrada cancelada.", file=sys.stderr)
        return
        
    input_path = Path(input_file)
    
    if not input_path.exists():
        print(f"ERRO: Arquivo '{input_file}' não encontrado.", file=sys.stderr)
        return

    aligned_file = f"{input_path.stem}.aligned.fasta"
    trimmed_file = f"{input_path.stem}.trimmed.fasta"
    
    # --- PASSO 1: ALINHAMENTO (MAFFT) ---
    print(f"\n[PASSO 1/3] Alinhando '{input_file}' com MAFFT...")
    # MUDANÇA AQUI: O comando agora é uma lista.
    # 'sh -c' recebe o comando com aspas como um *único argumento*.
    mafft_cmd = ['sh', '-c', f"mafft --auto /data/{input_path.name} > /data/{aligned_file}"]
    run_docker_command("mafft", mafft_cmd)

    # --- PASSO 2: TRIMAGEM (trimAl) ---
    print(f"\n[PASSO 2/3] Filtrando alinhamento com trimAl...")
    # MUDANÇA AQUI: Comando como lista
    trimal_cmd = ['trimal', '-in', f'/data/{aligned_file}', '-out', f'/data/{trimmed_file}', '-automated1']
    run_docker_command("trimal", trimal_cmd)

    # --- PASSO 3: ESCOLHA DO MÉTODO ---
    print(f"\n[PASSO 3/3] Alinhamento filtrado salvo em '{trimmed_file}'.")
    print("Qual método filogenético você deseja usar?")
    print("[1] Neighbor-Joining (NJ) - Rápido, para uma primeira olhada")
    print("[2] Máxima Verossimilhança (ML) - Robusto, padrão-ouro (IQ-TREE)")
    print("[3] Inferência Bayesiana (BI) - Robusto, mas interativo (MrBayes)")
    
    try:
        choice = input("Digite sua escolha (1, 2, ou 3): ")
    except EOFError:
        print("\nEntrada cancelada.", file=sys.stderr)
        return

    if choice == '1':
        # --- Rota 1: Neighbor-Joining (Biopython) ---
        print("\nExecutando Neighbor-Joining (NJ)...")
        nj_script_name = "_temp_nj_script.py"
        nj_output_file = f"{input_path.stem}.nj_tree.newick"
        
        nj_script_content = f"""
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import sys
sys.dont_write_bytecode = True
print('Lendo alinhamento...')
aln = AlignIO.read('/data/{trimmed_file}', 'fasta')
print('Calculando matriz de distância...')
calculator = DistanceCalculator('identity')
print('Construindo árvore NJ...')
constructor = DistanceTreeConstructor(calculator, 'nj')
tree = constructor.build_tree(aln)
tree.ladderize()
Phylo.write(tree, '/data/{nj_output_file}', 'newick')
print(f'Árvore NJ salva em: {nj_output_file}')
print('\\n--- Árvore ASCII ---')
Phylo.draw_ascii(tree)
"""
        with open(nj_script_name, "w") as f:
            f.write(nj_script_content)
        
        # MUDANÇA AQUI: Comando como lista
        biopython_cmd = ['python', f'/data/{nj_script_name}']
        run_docker_command("biopython", biopython_cmd)
        
        os.remove(nj_script_name)

    elif choice == '2':
        # --- Rota 2: Máxima Verossimilhança (IQ-TREE) ---
        print("\nExecutando Máxima Verossimilhança (ML)...")
        print("Isso pode demorar. O IQ-TREE fará o ModelFinder e 1000 bootstraps.")
        # MUDANÇA AQUI: Comando como lista
        iqtree_cmd = ['iqtree', '-s', f'/data/{trimmed_file}', '-m', 'MFP', '-B', '1000']
        run_docker_command("iqtree", iqtree_cmd)
        print(f"Árvore ML salva em: {trimmed_file}.treefile")

    elif choice == '3':
        # --- Rota 3: Inferência Bayesiana (MrBayes) ---
        print("\nIniciando setup para Inferência Bayesiana (BI)...")
        nexus_file = f"{input_path.stem}.trimmed.nex"
        
        print("Convertendo FASTA para NEXUS (usando Biopython)...")
        # MUDANÇA AQUI: Comando como lista
        biopython_convert_cmd = [
            'python', '-c', 
            f"from Bio import AlignIO; AlignIO.convert('/data/{trimmed_file}', 'fasta', '/data/{nexus_file}', 'nexus', molecule_type='DNA')"
        ]
        run_docker_command("biopython", biopython_convert_cmd)
        
        print(f"\nArquivo '{nexus_file}' criado.")
        print("Iniciando o MrBayes em modo interativo.")
        print("--- DENTRO DO MRBAYES, DIGITE OS SEGUINTES COMANDOS: ---")
        print(f"1. execute /data/{nexus_file}")
        print(f"2. mcmc ngen=100000 samplefreq=100")
        print(f"3. sumt")
        print(f"4. quit")
        print("---------------------------------------------------------")
        
        # MUDANÇA AQUI: Comando como lista
        mrbayes_cmd = ['mb']
        run_docker_command("mrbayes", mrbayes_cmd, interactive=True)
        print(f"Arquivos do MrBayes (incluindo .con.tre) devem estar no diretório.")

    else:
        print(f"Escolha '{choice}' inválida. Saindo.", file=sys.stderr)

    print("\nProjeto de Filogenia concluído!")


if __name__ == "__main__":
    main()