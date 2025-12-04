# gui_app.py (Usando NiceGUI)

from nicegui import ui
import subprocess
import os

# --- Ação que será chamada pelo botão ---
def run_pipeline_with_nicegui(fasta_file: str, method_choice: str):
    
    # 1. Defina o comando para chamar o seu script Python
    # O seu script DEVE estar adaptado para receber argumentos!
    command = [
        "python", "run_phylogeny.py",
        "--input", fasta_file,
        "--method", method_choice
    ]
    
    # 2. Mostre o comando e o status no log da NiceGUI
    ui.notify(f'Executando: {" ".join(command)}')
    
    # 3. Execute o script e capture a saída
    try:
        # Usa subprocess.run (bloqueante)
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        
        # Exibe a saída
        ui.label('--- Output do Pipeline ---')
        ui.code(result.stdout)
        ui.notify('Pipeline concluído com sucesso!', type='positive')

    except subprocess.CalledProcessError as e:
        ui.notify(f'ERRO no Pipeline! Código de saída: {e.returncode}', type='negative')
        ui.label('--- Erro do Pipeline ---')
        ui.code(e.stderr).classes('text-red')
        

# --- Layout da Interface (NiceGUI) ---

with ui.card():
    ui.label('Pipeline de Filogenia (Docker)')
    
    # Variáveis de entrada
    fasta_input = ui.input(label='Caminho do Arquivo FASTA', placeholder='ex: test.fasta')
    
    method_options = {
        '1': 'Neighbor-Joining (NJ)',
        '2': 'Máxima Verossimilhança (ML)',
        '3': 'Inferência Bayesiana (BI)'
    }
    
    method_select = ui.select(
        # Passa o dicionário
        options=method_options,
        # O valor inicial DEVE ser a chave do dicionário (a string '2')
        value='2',
        label='Método Filogenético'
    )
    
    # Botão de Execução
    ui.button('INICIAR PIPELINE', on_click=lambda: run_pipeline_with_nicegui(
        fasta_input.value, method_select.value
    ))

ui.run() # Inicia a aplicação web local