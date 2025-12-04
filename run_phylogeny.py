import os
import subprocess
import sys
from pathlib import Path
import threading
import customtkinter as ctk
from tkinter import filedialog, messagebox

# ==============================================================================
# --- 1. CONFIGURAÇÕES E FUNÇÕES DE AJUDA (Adaptação do seu código) ---
# ==============================================================================

# --- Configuração dos Nomes dos Containers ---
CONTAINERS = {
    "mafft": "quay.io/biocontainers/mafft:7.525--h031d066_0",
    "trimal": "quay.io/biocontainers/trimal:1.5",
    "iqtree": "quay.io/biocontainers/iqtree:3.0.1--h503566f_0",
    "biopython": "quay.io/biocontainers/biopython:1.79",
    "mrbayes": "quay.io/biocontainers/mrbayes:3.2.7--hd0d793b_7"
}

# --- Função Auxiliar para Rodar Comandos Docker (ADAPTADA para GUI) ---
def run_docker_command(container_name: str, command: list, log_callback, interactive: bool = False):
    """
    Executa um comando dentro de um container Docker e reporta o output completo 
    para a GUI via log_callback. Levanta erro em caso de falha.
    """
    work_dir = Path.cwd().resolve()
    
    docker_base = [
        "docker", "run", "--rm",
        "-v", f"{work_dir}:/data",
        "-w", "/data"
    ]
    
    if interactive:
        # Modo interativo (para MrBayes)
        docker_base.extend(["-it"]) 
        
    full_command = docker_base + [CONTAINERS[container_name]] + command

    log_callback(f"\n--- Executando: {' '.join(full_command)} ---", is_highlight=True)
    
    try:
        # Comando bloqueante, mas executado em Thread separada.
        process = subprocess.run(
            full_command, 
            text=True, 
            check=True, 
            capture_output=True
        )
        
        # Reporta o output completo (STDOUT)
        if process.stdout:
             log_callback("--- SAÍDA DO CONTAINER (STDOUT) ---")
             log_callback(process.stdout)
        
        log_callback("--- Comando concluído com sucesso! ---", is_success=True)
        return True
        
    except subprocess.CalledProcessError as e:
        # Em caso de erro do comando
        log_callback(f"ERRO DE EXECUÇÃO: Comando {container_name} falhou.", is_error=True)
        log_callback(f"STATUS CODE: {e.returncode}", is_error=True)
        if e.stdout: log_callback(f"STDOUT:\n{e.stdout}", is_error=True)
        if e.stderr: log_callback(f"STDERR (ERRO):\n{e.stderr}", is_error=True)
        
        # Levanta a exceção para que a função chamadora interrompa o pipeline
        raise 
        
    except Exception as e:
        # Outros erros (ex: Docker não instalado)
        log_callback(f"ERRO DESCONHECIDO NO DOCKER: {e}", is_error=True)
        raise

# ==============================================================================
# --- 2. CLASSE DA APLICAÇÃO GUI (CustomTkinter) ---
# ==============================================================================

class PhylogenyApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        # --- Configuração Básica da Janela ---
        self.title("Pipeline de Filogenia (Docker) - por Walter")
        self.geometry("800x650")
        ctk.set_appearance_mode("System")

        # Variáveis de Controle
        self.input_file_path = ctk.StringVar(value="")
        self.choice_var = ctk.StringVar(value="2") # ML (IQ-TREE) como padrão
        
        # Estrutura do Layout (Grid)
        self.grid_columnconfigure(0, weight=1) 
        self.grid_rowconfigure(2, weight=1)   

        # 1. Frame de Entrada (Input Frame)
        self.input_frame = ctk.CTkFrame(self)
        self.input_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        self.create_input_widgets(self.input_frame)

        # 2. Frame de Opções (Options Frame)
        self.options_frame = ctk.CTkFrame(self)
        self.options_frame.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        self.create_options_widgets(self.options_frame)
        
        # 3. Frame de Log (Log Frame)
        self.log_frame = ctk.CTkFrame(self)
        self.log_frame.grid(row=2, column=0, padx=10, pady=10, sticky="nsew")
        self.create_log_widgets(self.log_frame)

    # --- Criação de Widgets ---
    
    def create_input_widgets(self, parent):
        parent.grid_columnconfigure(1, weight=1)
        
        ctk.CTkLabel(parent, text="Arquivo Multi-FASTA:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        
        # Campo de Texto para exibir o caminho
        ctk.CTkEntry(parent, textvariable=self.input_file_path, state="readonly").grid(row=0, column=1, padx=5, pady=5, sticky="ew")
        
        # Botão 'Browse'
        ctk.CTkButton(parent, text="Selecionar Arquivo", command=self.select_file).grid(row=0, column=2, padx=5, pady=5)
        
        # Botão 'Executar'
        self.start_button = ctk.CTkButton(parent, text="▶️ INICIAR PIPELINE", command=self.start_pipeline, fg_color="#1F538D") # Cor azul personalizada
        self.start_button.grid(row=1, column=0, columnspan=3, pady=10)

    def create_options_widgets(self, parent):
        parent.grid_columnconfigure((0, 1, 2, 3), weight=1)
        ctk.CTkLabel(parent, text="Método Filogenético:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        
        # Radio Buttons
        options = [
            ("1 - Neighbor-Joining (NJ)", '1'),
            ("2 - Máxima Verossimilhança (ML) - IQ-TREE", '2'),
            ("3 - Inferência Bayesiana (BI) - MrBayes (Interativo)", '3')
        ]
        
        for i, (text, value) in enumerate(options):
            ctk.CTkRadioButton(parent, text=text, variable=self.choice_var, value=value).grid(row=0, column=i+1, padx=10, pady=5, sticky="w")

    def create_log_widgets(self, parent):
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(1, weight=1)
        
        ctk.CTkLabel(parent, text="Log de Execução:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        
        # Campo de Texto para Log (ReadOnly)
        self.log_text = ctk.CTkTextbox(parent, wrap="word", height=250)
        self.log_text.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")
        self.log_text.configure(state="disabled")
        
        # Configurações de tag para cores no log
        self.log_text.tag_config('error', foreground='red')
        self.log_text.tag_config('success', foreground='green')
        self.log_text.tag_config('highlight', foreground='blue') # Para comandos

    # --- Funções de Lógica da GUI ---

    def select_file(self):
        filepath = filedialog.askopenfilename(
            title="Selecione o Arquivo Multi-FASTA",
            filetypes=(("FASTA files", "*.fasta *.fa"), ("All files", "*.*"))
        )
        if filepath:
            # Garante que o diretório de trabalho seja o do arquivo (simplifica a montagem Docker)
            os.chdir(Path(filepath).parent) 
            self.input_file_path.set(Path(filepath).name)
            self.log("Arquivo selecionado e diretório de trabalho definido para: " + str(Path(filepath).parent))

    def log(self, message: str, is_error: bool = False, is_success: bool = False, is_highlight: bool = False):
        """Atualiza a caixa de texto de log de forma segura."""
        self.log_text.configure(state="normal")
        
        tag = ''
        if is_error:
            tag = 'error'
        elif is_success:
            tag = 'success'
        elif is_highlight:
            tag = 'highlight'
            
        self.log_text.insert("end", message + "\n", tag)
        self.log_text.see("end") 
        self.log_text.configure(state="disabled")
        self.update() 

    def toggle_interface(self, state):
        """Habilita ou desabilita botões e opções."""
        self.start_button.configure(state=state)
        # Habilita/desabilita as opções (RadioButtons)
        for widget in self.options_frame.winfo_children():
            if isinstance(widget, ctk.CTkRadioButton):
                widget.configure(state=state)

    def start_pipeline(self):
        # Checa e prepara o log
        input_file = self.input_file_path.get()
        if not input_file:
            messagebox.showerror("Erro de Input", "Por favor, selecione um arquivo FASTA de entrada.")
            return

        # Limpar Log e Desabilitar Interface antes de começar
        self.toggle_interface(state="disabled")
        self.log_text.configure(state="normal")
        self.log_text.delete(1.0, "end")
        self.log_text.configure(state="disabled")
        
        self.log("Iniciando o Pipeline de Filogenia (usando Thread para não travar a GUI)...", is_highlight=True)
        
        # Cria e inicia a thread para rodar a lógica principal
        pipeline_thread = threading.Thread(target=self.run_pipeline_logic, daemon=True)
        pipeline_thread.start()

    # --- Função Principal do Pipeline (Executada em Thread) ---

    def run_pipeline_logic(self):
        # AQUI usamos o Path().name pois o diretório já foi definido e o Docker 
        # montará o diretório de trabalho atual (cwd)
        input_file_name = self.input_file_path.get()
        
        aligned_file = Path(input_file_name).with_suffix('').name + ".aligned.fasta"
        trimmed_file = Path(input_file_name).with_suffix('').name + ".trimmed.fasta"
        choice = self.choice_var.get()
        
        try:
            # --- PASSO 1: ALINHAMENTO (MAFFT) ---
            self.log(f"\n[PASSO 1/3] Alinhando '{input_file_name}' com MAFFT...")
            mafft_cmd = ['sh', '-c', f"mafft --auto /data/{input_file_name} > /data/{aligned_file}"]
            run_docker_command("mafft", mafft_cmd, self.log)

            # --- PASSO 2: TRIMAGEM (trimAl) ---
            self.log(f"\n[PASSO 2/3] Filtrando alinhamento com trimAl...")
            trimal_cmd = ['trimal', '-in', f'/data/{aligned_file}', '-out', f'/data/{trimmed_file}', '-automated1']
            run_docker_command("trimal", trimal_cmd, self.log)

            # --- PASSO 3: CONSTRUÇÃO DA ÁRVORE ---
            self.log(f"\n[PASSO 3/3] Executando o método de escolha (Opção {choice})...")

            if choice == '1':
                # --- Rota 1: Neighbor-Joining (Biopython) ---
                self.log("Executando Neighbor-Joining (NJ)...")
                nj_script_name = "_temp_nj_script.py"
                nj_output_file = Path(input_file_name).with_suffix('').name + ".nj_tree.newick"
                
                # Conteúdo do script Biopython
                nj_script_content = f"""
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import sys
sys.dont_write_bytecode = True
print('Lendo alinhamento...')
aln = AlignIO.read('/data/{trimmed_file}', 'fasta')
print('Calculando matriz de distância...')
calculator = DistanceCalculator('identity')
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
                
                biopython_cmd = ['python', f'/data/{nj_script_name}']
                run_docker_command("biopython", biopython_cmd, self.log)
                
                os.remove(nj_script_name)
                self.log(f"Árvore NJ salva em: {nj_output_file}", is_success=True)

            elif choice == '2':
                # --- Rota 2: Máxima Verossimilhança (IQ-TREE) ---
                self.log("Executando Máxima Verossimilhança (ML) - IQ-TREE (ModelFinder e 1000 Boostraps)...")
                iqtree_cmd = ['iqtree', '-s', f'/data/{trimmed_file}', '-m', 'MFP', '-B', '1000']
                run_docker_command("iqtree", iqtree_cmd, self.log)
                self.log(f"Árvore ML salva em: {trimmed_file}.treefile", is_success=True)

            elif choice == '3':
                # --- Rota 3: Inferência Bayesiana (MrBayes) ---
                self.log("Iniciando setup para Inferência Bayesiana (BI)...")
                nexus_file = Path(input_file_name).with_suffix('').name + ".trimmed.nex"
                
                # Conversão FASTA para NEXUS
                self.log("Convertendo FASTA para NEXUS (usando Biopython)...")
                biopython_convert_cmd = [
                    'python', '-c', 
                    f"from Bio import AlignIO; AlignIO.convert('/data/{trimmed_file}', 'fasta', '/data/{nexus_file}', 'nexus', molecule_type='DNA')"
                ]
                run_docker_command("biopython", biopython_convert_cmd, self.log)
                
                self.log("\n*** MODO INTERATIVO DO MRBAYES INICIADO ***", is_error=True)
                self.log("Atenção: MrBayes requer interação manual no terminal para rodar o MCMC e SUMT.")
                self.log("Comandos recomendados:", is_highlight=True)
                self.log(f"1. execute /data/{nexus_file}")
                self.log(f"2. lset nst=6 rates=gamma") # Exemplo de modelo
                self.log(f"3. mcmc ngen=100000 samplefreq=100")
                self.log(f"4. sumt")
                self.log(f"5. quit")

                mrbayes_cmd = ['mb']
                # Esta chamada VAI travar o terminal onde a GUI foi chamada!
                run_docker_command("mrbayes", mrbayes_cmd, self.log, interactive=True)
                
            else:
                self.log(f"Escolha '{choice}' inválida. Saindo.", is_error=True)

            self.log("\n*** PROJETO DE FILOGENIA CONCLUÍDO COM SUCESSO! ***", is_success=True)
            messagebox.showinfo("Sucesso", "Pipeline de Filogenia concluído!")

        except Exception:
            # O erro já foi reportado pelo log na função run_docker_command
            messagebox.showerror("Erro no Pipeline", "O pipeline foi interrompido devido a um erro. Verifique o log para detalhes.")

        finally:
            # Re-habilita a interface, independente de sucesso ou falha
            self.toggle_interface(state="normal")
            self.log("Interface re-habilitada para nova execução.")


# ==============================================================================
# --- 3. EXECUÇÃO ---
# ==============================================================================

if __name__ == "__main__":
    # É necessário que o Docker esteja rodando e o CustomTkinter esteja instalado!
    
    # Define o diretório inicial para a GUI (o diretório onde este script está)
    os.chdir(Path(__file__).parent) 
    
    app = PhylogenyApp()
    app.mainloop()