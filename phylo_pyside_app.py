import os
import subprocess
import sys
from pathlib import Path

# PySide6 Imports
from PySide6.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, 
                               QLabel, QLineEdit, QPushButton, QRadioButton, 
                               QGroupBox, QTextEdit, QFileDialog, QMessageBox, QButtonGroup)
from PySide6.QtCore import QThread, Signal, QObject, Qt # Qt is required for constants
from PySide6.QtGui import QFont, QColor, QTextCharFormat, QTextCursor # QTextCursor is now used directly

# ==============================================================================
# --- 1. CONFIGURATION AND HELPER FUNCTIONS ---
# ==============================================================================

# --- Container Names Configuration ---
CONTAINERS = {
    "mafft": "quay.io/biocontainers/mafft:7.525--h031d066_0",
    "trimal": "quay.io/biocontainers/trimal:1.5",
    "iqtree": "quay.io/biocontainers/iqtree:3.0.1--h503566f_0",
    "biopython": "quay.io/biocontainers/biopython:1.79",
    "mrbayes": "quay.io/biocontainers/mrbayes:3.2.7--hd0d793b_7"
}

# --- Docker Command Runner (Using Signal for Logging) ---
def run_docker_command(container_name: str, command: list, log_callback, interactive: bool = False):
    """
    Executes a Docker command. log_callback (a PySide Signal) is used to 
    send messages back to the GUI.
    """
    work_dir = Path.cwd().resolve()
    docker_base = [
        "docker", "run", "--rm",
        "-v", f"{work_dir}:/data",
        "-w", "/data"
    ]
    if interactive:
        docker_base.extend(["-it"]) 
        
    full_command = docker_base + [CONTAINERS[container_name]] + command

    log_callback.emit(f"\n--- Executing: {' '.join(full_command)} ---", 'highlight')
    
    try:
        process = subprocess.run(
            full_command, 
            text=True, 
            check=True, 
            capture_output=True,
            timeout=None 
        )
        
        if process.stdout:
             log_callback.emit("--- CONTAINER OUTPUT (STDOUT) ---", 'default')
             log_callback.emit(process.stdout, 'default')
        
        log_callback.emit("--- Command completed successfully! ---", 'success')
        return True
        
    except subprocess.CalledProcessError as e:
        log_callback.emit(f"EXECUTION ERROR: {container_name} command failed.", 'error')
        log_callback.emit(f"STATUS CODE: {e.returncode}", 'error')
        if e.stderr: log_callback.emit(f"STDERR (ERROR):\n{e.stderr}", 'error')
        raise 
        
    except Exception as e:
        log_callback.emit(f"UNKNOWN DOCKER ERROR: {e}", 'error')
        raise

# ==============================================================================
# --- 2. WORKER (Pipeline Logic in Separate Thread) ---
# ==============================================================================

class PhylogenyWorker(QObject):
    # Signals definition
    log_message = Signal(str, str) # Message and Type (default, success, error, highlight)
    finished = Signal(bool)        # True for success, False for error/canceled

    def __init__(self, input_file_name, choice):
        super().__init__()
        self.input_file_name = input_file_name
        self.choice = choice

    def run(self):
        """Contains the complete pipeline logic, including MrBayes automation."""
        
        input_file_name = self.input_file_name
        
        base_name = Path(input_file_name).with_suffix('').name
        aligned_file = base_name + ".aligned.fasta"
        trimmed_file = base_name + ".trimmed.fasta"
        
        try:
            # --- STEP 1: ALIGNMENT (MAFFT) ---
            self.log_message.emit(f"\n[STEP 1/3] Aligning '{input_file_name}' with MAFFT...", 'default')
            mafft_cmd = ['sh', '-c', f"mafft --auto /data/{input_file_name} > /data/{aligned_file}"]
            run_docker_command("mafft", mafft_cmd, self.log_message)

            # --- STEP 2: TRIMMING (trimAl) ---
            self.log_message.emit(f"\n[STEP 2/3] Filtering alignment with trimAl...", 'default')
            trimal_cmd = ['trimal', '-in', f'/data/{aligned_file}', '-out', f'/data/{trimmed_file}', '-automated1']
            run_docker_command("trimal", trimal_cmd, self.log_message)

            # --- STEP 3: TREE CONSTRUCTION ---
            self.log_message.emit(f"\n[STEP 3/3] Executing selected method (Option {self.choice})...", 'default')

            if self.choice == '1':
                # --- Route 1: Neighbor-Joining (Biopython) ---
                self.log_message.emit("Executing Neighbor-Joining (NJ)...", 'default')
                nj_script_name = "_temp_nj_script.py"
                nj_output_file = base_name + ".nj_tree.newick"
                
                # Biopython script content
                nj_script_content = f"""
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import sys
sys.dont_write_bytecode = True
print('Reading alignment...')
aln = AlignIO.read('/data/{trimmed_file}', 'fasta')
print('Calculating distance matrix...')
calculator = DistanceCalculator('identity')
constructor = DistanceTreeConstructor(calculator, 'nj')
tree = constructor.build_tree(aln)
tree.ladderize()
Phylo.write(tree, '/data/{nj_output_file}', 'newick')
print(f'NJ Tree saved to: {nj_output_file}')
print('\\n--- ASCII Tree ---')
Phylo.draw_ascii(tree)
"""
                with open(nj_script_name, "w") as f:
                    f.write(nj_script_content)
                
                biopython_cmd = ['python', f'/data/{nj_script_name}']
                run_docker_command("biopython", biopython_cmd, self.log_message)
                
                os.remove(nj_script_name)
                self.log_message.emit(f"NJ Tree saved to: {nj_output_file}", 'success')

            elif self.choice == '2':
                # --- Route 2: Maximum Likelihood (IQ-TREE) ---
                self.log_message.emit("Executing Maximum Likelihood (ML) - IQ-TREE (ModelFinder and 1000 Boostraps)...", 'default')
                iqtree_cmd = ['iqtree', '-s', f'/data/{trimmed_file}', '-m', 'MFP', '-B', '1000']
                run_docker_command("iqtree", iqtree_cmd, self.log_message)
                self.log_message.emit(f"ML Tree saved to: {trimmed_file}.treefile", 'success')

            elif self.choice == '3':
                # --- Route 3: Bayesian Inference (MrBayes) - AUTOMATION ---
                self.log_message.emit("Starting setup for Bayesian Inference (BI)...", 'default')
                nexus_file = base_name + ".trimmed.nex"
                
                # STEP 3.1: Convert FASTA to NEXUS
                self.log_message.emit("Converting FASTA to NEXUS (using Biopython)...", 'default')
                biopython_convert_cmd = [
                    'python', '-c', 
                    f"from Bio import AlignIO; AlignIO.convert('/data/{trimmed_file}', 'fasta', '/data/{nexus_file}', 'nexus', molecule_type='DNA')"
                ]
                run_docker_command("biopython", biopython_convert_cmd, self.log_message)
                
                # STEP 3.2: Create MrBayes Command Block (Automation)
                mrbayes_commands = f"""
[ MrBayes command block for automation ]

begin mrbayes;
    set autoclose=yes nowarn=yes; 
    
    # Substitution model (GTR+G, common for DNA)
    lset nst=6 rates=gamma;
    
    # MCMC Configuration 
    mcmc ngen=100000 samplefreq=100 nchains=4 nruns=2;
    
    # Summarize (discards the first 25% as burnin)
    sump burnin=250;
    sumt burnin=250;
    
    # Exit automatically
    quit;
end;
"""
                # Append commands to the NEXUS file
                try:
                    with open(nexus_file, 'a') as f:
                        f.write(mrbayes_commands)
                    self.log_message.emit(f"Automated MrBayes commands added to {nexus_file}", 'success')
                except Exception as e:
                    self.log_message.emit(f"ERROR: Failed to write commands to NEXUS: {e}", 'error')
                    raise

                # STEP 3.3: Automated MrBayes Execution
                self.log_message.emit("\nStarting MrBayes in NON-interactive and automated mode...", 'default')
                
                mrbayes_cmd = ['mb', f'/data/{nexus_file}']
                
                run_docker_command("mrbayes", mrbayes_cmd, self.log_message)
                
                self.log_message.emit(f"Bayesian Inference finished. Consensus tree saved to '{nexus_file}.con.tre'", 'success')

            self.finished.emit(True)

        except Exception:
            self.finished.emit(False)


# ==============================================================================
# --- 3. INTERFACE CLASS (PySide6) ---
# ==============================================================================

class PhylogenyApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Phylo-Pipe")
        self.setGeometry(100, 100, 900, 700)
        
        self.thread = None
        self.worker = None

        # CREATION OF THE LOG WIDGET
        self.log_text = QTextEdit(self) 
        
        self.set_colors()
        self.setup_ui()
        self.apply_stylesheet() 

    def apply_stylesheet(self):
        """Applies a Qt stylesheet (CSS-like) for the White/Dark Blue theme."""
        
        # Colors for the theme
        WHITE = "#FFFFFF"
        LIGHT_GRAY = "#F0F0F0"
        DARK_BLUE = "#1F456E" # Main dark blue
        DEEP_BLUE = "#1A3653"
        BUTTON_HOVER = "#2A68A6"

        stylesheet = f"""
        /* Main Window Style */
        QWidget {{
            background-color: {WHITE};
            color: #333333; 
            font-family: Arial, sans-serif;
        }}

        /* Groups (QGroupBox) */
        QGroupBox {{
            background-color: {LIGHT_GRAY};
            border: 1px solid {DARK_BLUE};
            border-radius: 5px;
            margin-top: 10px;
            padding: 10px;
            font-weight: bold;
            color: {DEEP_BLUE};
        }}
        QGroupBox::title {{
            subcontrol-origin: margin;
            subcontrol-position: top left;
            padding: 0 3px;
        }}

        /* START Button */
        #startButton {{ 
            background-color: {DARK_BLUE};
            color: {WHITE};
            border: 1px solid {DEEP_BLUE};
            border-radius: 5px;
            padding: 15px 20px;
            font-size: 16pt;
            font-weight: bold;
        }}
        #startButton:hover {{
            background-color: {BUTTON_HOVER};
        }}

        /* Normal Button (Select File) */
        #selectFileButton {{ /* Use the new objectName */
            background-color: {DARK_BLUE};
            color: {WHITE};
            border: 1px solid {DEEP_BLUE};
            border-radius: 3px;
            padding: 5px 10px;
        }}
        #selectFileButton:hover {{
            background-color: {BUTTON_HOVER};
        }}

        /* Line Edit (QLineEdit) */
        QLineEdit {{
            background-color: {WHITE};
            border: 1px solid #CCCCCC;
            border-radius: 3px;
            padding: 5px;
        }}
        
        /* Log Area (QTextEdit) */
        QTextEdit {{
            background-color: #F8F8F8;
            border: 1px solid #CCCCCC;
            padding: 5px;
        }}
        """
        self.setStyleSheet(stylesheet)

    def set_colors(self):
        """Defines colors and tags for the log."""
        self.log_text.setFont(QFont("Monospace", 9))
        
        self.fmt_error = QTextCharFormat()
        self.fmt_error.setForeground(QColor("red"))
        
        self.fmt_success = QTextCharFormat()
        self.fmt_success.setForeground(QColor("green"))
        
        self.fmt_highlight = QTextCharFormat()
        self.fmt_highlight.setForeground(QColor("blue")) 

    def setup_ui(self):
        main_layout = QVBoxLayout(self)

        # 1. INPUT
        input_group = QGroupBox("1. Input File", self)
        input_layout = QHBoxLayout(input_group)
        self.input_line = QLineEdit(self)
        self.input_line.setReadOnly(True) 
        input_layout.addWidget(QLabel("FASTA Path:"))
        input_layout.addWidget(self.input_line)
        
        browse_button = QPushButton("Select File")
        browse_button.setObjectName("selectFileButton") # 🌟 CORRECTION: Object name added
        browse_button.clicked.connect(self.select_file)
        
        input_layout.addWidget(browse_button)
        main_layout.addWidget(input_group)

        # 2. OPTIONS
        options_group = QGroupBox("2. Phylogenetic Method", self)
        options_layout = QHBoxLayout(options_group)
        self.radio_group = QButtonGroup(self)
        
        self.radio_nj = QRadioButton("1 - Neighbor-Joining (NJ)")
        self.radio_ml = QRadioButton("2 - Maximum Likelihood (ML) - IQ-TREE")
        self.radio_bi = QRadioButton("3 - Bayesian Inference (BI) - MrBayes (Automated)")
        
        self.radio_ml.setChecked(True) 
        
        self.radio_group.addButton(self.radio_nj, 1)
        self.radio_group.addButton(self.radio_ml, 2)
        self.radio_group.addButton(self.radio_bi, 3)
        
        options_layout.addWidget(self.radio_nj)
        options_layout.addWidget(self.radio_ml)
        options_layout.addWidget(self.radio_bi)
        options_layout.addStretch(1) 
        main_layout.addWidget(options_group)

        # 3. CONTROL & LOG
        control_layout = QHBoxLayout()
        self.start_button = QPushButton("RUN")
        self.start_button.setObjectName("startButton") 
        self.start_button.clicked.connect(self.start_pipeline)
        control_layout.addWidget(self.start_button)
        main_layout.addLayout(control_layout)
        
        log_label = QLabel("Execution Log:")
        main_layout.addWidget(log_label)

        self.log_text.setReadOnly(True)
        main_layout.addWidget(self.log_text)

    # --- Event and Log Functions ---

    def select_file(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Select Multi-FASTA File", "", "FASTA Files (*.fasta *.fa);;All Files (*)")
        if filepath:
            os.chdir(Path(filepath).parent)
            self.input_line.setText(Path(filepath).name)
            self.log(f"File selected. Working directory: {Path(filepath).parent}", 'default')

    def log(self, message: str, message_type: str):
        """Updates the log box with color formatting. Fixes the cursor error."""
        cursor = self.log_text.textCursor()
        
        # 🌟 CORRECTION: Use QTextCursor.End directly, as the class is imported.
        cursor.movePosition(QTextCursor.End) 
        
        # Select color format
        if message_type == 'error':
            fmt = self.fmt_error
        elif message_type == 'success':
            fmt = self.fmt_success
        elif message_type == 'highlight':
            fmt = self.fmt_highlight
        else: # default
            fmt = QTextCharFormat()

        cursor.insertText(message + "\n", fmt)
        self.log_text.setTextCursor(cursor) 

    def toggle_interface(self, enable: bool):
        """Enables/disables controls."""
        self.start_button.setEnabled(enable)
        self.radio_nj.setEnabled(enable)
        self.radio_ml.setEnabled(enable)
        self.radio_bi.setEnabled(enable)
        self.input_line.setReadOnly(not enable)
        
        # 🌟 CORRECTION: Use the defined objectName and check for NoneType
        select_button = self.findChild(QPushButton, "selectFileButton") 
        if select_button: 
            select_button.setEnabled(enable)

    def pipeline_finished(self, success: bool):
        """Called when the Worker thread finishes."""
        self.toggle_interface(True)
        self.thread.quit()
        self.thread.wait() 
        
        if success:
            QMessageBox.information(self, "Success", "Phylogeny Pipeline completed successfully!")
        else:
            QMessageBox.critical(self, "Fatal Error", "The pipeline was interrupted due to an error. Check the log.")

    def start_pipeline(self):
        input_file_name = self.input_line.text()
        if not input_file_name:
            QMessageBox.warning(self, "Attention", "Please select a FASTA input file.")
            return

        choice = self.radio_group.checkedId()
        
        self.log_text.clear()
        self.toggle_interface(False)
        self.log("Initializing pipeline...", 'highlight')
        
        # Create Worker and Thread
        self.thread = QThread()
        self.worker = PhylogenyWorker(input_file_name, str(choice))
        
        # Connect and start
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.pipeline_finished)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.worker.log_message.connect(self.log) 
        
        self.thread.start()


# ==============================================================================
# --- 4. EXECUTION ---
# ==============================================================================

if __name__ == "__main__":
    
    # Set the working directory to the script location
    os.chdir(Path(__file__).parent) 
    
    app = QApplication(sys.argv)
    window = PhylogenyApp()
    window.show()
    sys.exit(app.exec())