import re
import os
from Bio import Entrez, SeqIO
from collections import defaultdict

# --- CONFIGURAÇÃO ---
Entrez.email = "waltergomesbio@gmail.com"
INPUT_FILE = "raw_table.txt"
OUTPUT_DIR = "dataset_lamiaceae"
BATCH_SIZE = 100 

def extract_ids(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    # Padrão GenBank: 2 letras + 6 números ou 1 letra + 5 números
    pattern = r'\b([A-Z]{2}\d{6}|[A-Z]{1}\d{5})\b'
    ids = list(set(re.findall(pattern, content)))
    print(f"--> {len(ids)} IDs únicos encontrados em {filepath}")
    return ids

def download_and_sort(id_list):
    gene_buckets = defaultdict(list)
    total = len(id_list)
    
    print(f"--> Iniciando download de {total} sequências via NCBI Entrez...")

    for i in range(0, total, BATCH_SIZE):
        batch = id_list[i:i+BATCH_SIZE]
        try:
            handle = Entrez.efetch(db="nucleotide", id=",".join(batch), rettype="fasta", retmode="text")
            records = list(SeqIO.parse(handle, "fasta"))
            handle.close()

            for record in records:
                desc = record.description.lower()
                # Triagem baseada nos 5 genes do artigo (matK, ndhF, rbcL, rps16, trnL-F)
                if "rbcl" in desc or "ribulose" in desc:
                    gene_buckets["rbcL"].append(record)
                elif "matk" in desc or "maturase" in desc:
                    gene_buckets["matK"].append(record)
                elif "ndhf" in desc or "nadh" in desc:
                    gene_buckets["ndhF"].append(record)
                elif "rps16" in desc:
                    gene_buckets["rps16"].append(record)
                elif "trnl" in desc or "trnf" in desc or "intergenic" in desc:
                    gene_buckets["trnL-F"].append(record)
                else:
                    gene_buckets["outros"].append(record)
            
            print(f"   Progresso: {min(i+BATCH_SIZE, total)}/{total} baixados...")
        except Exception as e:
            print(f"   Erro no lote {i}: {e}")

    return gene_buckets

def save_files(buckets):
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    
    for gene, records in buckets.items():
        if records:
            path = os.path.join(OUTPUT_DIR, f"{gene}.fasta")
            SeqIO.write(records, path, "fasta")
            print(f"--> Salvo: {path} ({len(records)} sequências)")

if __name__ == "__main__":
    if not os.path.exists(INPUT_FILE):
        print(f"Erro: Arquivo {INPUT_FILE} não encontrado!")
    else:
        ids = extract_ids(INPUT_FILE)
        if ids:
            buckets = download_and_sort(ids)
            save_files(buckets)
            print("\nConcluído! Seus arquivos FASTA estão em case_study/dataset_lamiaceae/")