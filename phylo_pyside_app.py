import os
import subprocess
import sys
import json
from pathlib import Path
from datetime import datetime
from collections import Counter

# PySide6 Imports
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, 
    QPushButton, QRadioButton, QGroupBox, QTextEdit, QFileDialog, 
    QMessageBox, QButtonGroup, QProgressBar, QTabWidget, QSpinBox,
    QDoubleSpinBox, QCheckBox, QTableWidget, QTableWidgetItem, QComboBox,
    QScrollArea, QFrame, QSplitter
)
from PySide6.QtCore import QThread, Signal, QObject, Qt, QTimer
from PySide6.QtGui import QFont, QColor, QTextCharFormat, QTextCursor

# ==============================================================================
# --- 1. CONFIGURATION AND HELPER FUNCTIONS ---
# ==============================================================================

CONTAINERS = {
    "mafft": "quay.io/biocontainers/mafft:7.525--h031d066_0",
    "trimal": "quay.io/biocontainers/trimal:1.5",
    "iqtree": "quay.io/biocontainers/iqtree:3.0.1--h503566f_0",
    "biopython": "quay.io/biocontainers/biopython:1.79",
    "mrbayes": "quay.io/biocontainers/mrbayes:3.2.7--hd0d793b_7",
    "cladegene_plot": "cladegene_plot:latest"
}

HISTORY_FILE = "phylo_history.json"

# --- Sequence Validation Functions ---
def validate_fasta(filepath):
    """Validates FASTA file and returns statistics."""
    issues = []
    sequences = {}
    
    try:
        with open(filepath, 'r') as f:
            current_id = None
            current_seq = []
            
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:].split()[0]
                    if current_id in sequences:
                        issues.append(f"Duplicate ID: {current_id}")
                    current_seq = []
                else:
                    # Check for invalid characters
                    invalid_chars = set(line.upper()) - set('ACGTRYSWKMBDHVN-')
                    if invalid_chars:
                        issues.append(f"Line {line_num}: Invalid characters {invalid_chars}")
                    current_seq.append(line)
            
            if current_id:
                sequences[current_id] = ''.join(current_seq)
        
        if not sequences:
            issues.append("No sequences found in file")
            return False, issues, {}
        
        # Calculate statistics
        seq_lengths = [len(seq) for seq in sequences.values()]
        stats = {
            'num_sequences': len(sequences),
            'min_length': min(seq_lengths),
            'max_length': max(seq_lengths),
            'avg_length': sum(seq_lengths) / len(seq_lengths),
            'total_length': sum(seq_lengths)
        }
        
        return len(issues) == 0, issues, stats
        
    except Exception as e:
        return False, [f"Error reading file: {str(e)}"], {}

def analyze_alignment(filepath):
    """Analyzes alignment quality."""
    try:
        sequences = {}
        with open(filepath, 'r') as f:
            current_id = None
            current_seq = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:]
                    current_seq = []
                elif line:
                    current_seq.append(line)
            if current_id:
                sequences[current_id] = ''.join(current_seq)
        
        if not sequences:
            return {}
        
        alignment_length = len(list(sequences.values())[0])
        total_gaps = sum(seq.count('-') for seq in sequences.values())
        total_positions = len(sequences) * alignment_length
        
        # Calculate identity
        identity_sum = 0
        comparisons = 0
        seq_list = list(sequences.values())
        for i in range(len(seq_list)):
            for j in range(i+1, len(seq_list)):
                matches = sum(1 for a, b in zip(seq_list[i], seq_list[j]) 
                            if a == b and a != '-')
                identity_sum += matches / alignment_length
                comparisons += 1
        
        avg_identity = (identity_sum / comparisons * 100) if comparisons > 0 else 0
        
        return {
            'alignment_length': alignment_length,
            'gap_percentage': (total_gaps / total_positions * 100),
            'avg_identity': avg_identity
        }
    except:
        return {}

# --- Docker Command Runner ---
def run_docker_command(container_name, command, log_callback, progress_callback=None):
    """Executes Docker command with progress updates."""
    work_dir = Path.cwd().resolve()
    docker_base = [
        "docker", "run", "--rm",
        "-v", f"{work_dir}:/data",
        "-w", "/data"
    ]
    
    full_command = docker_base + [CONTAINERS[container_name]] + command
    log_callback.emit(f"\n🔧 Executing: {' '.join(full_command)}", 'highlight')
    
    try:
        process = subprocess.Popen(
            full_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        stdout_lines = []
        stderr_lines = []
        
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                stdout_lines.append(output.strip())
                log_callback.emit(output.strip(), 'default')
                if progress_callback:
                    progress_callback.emit(5)  # Increment progress
        
        stderr_output = process.stderr.read()
        if stderr_output:
            stderr_lines.append(stderr_output.strip())
        
        returncode = process.poll()
        
        if returncode != 0:
            log_callback.emit(f" ERROR: Command failed with code {returncode}", 'error')
            if stderr_output:
                log_callback.emit(f"STDERR:\n{stderr_output}", 'error')
            raise subprocess.CalledProcessError(returncode, full_command, stderr=stderr_output)
        
        log_callback.emit(" Command completed successfully!", 'success')
        return True
        
    except Exception as e:
        log_callback.emit(f" DOCKER ERROR: {e}", 'error')
        raise

# ==============================================================================
# --- 2. WORKER (Pipeline Logic) ---
# ==============================================================================

class PhylogenyWorker(QObject):
    log_message = Signal(str, str)
    progress_update = Signal(int)
    progress_increment = Signal(int)
    stage_update = Signal(str)
    finished = Signal(bool, dict)
    
    def __init__(self, input_file_name, choice, params):
        super().__init__()
        self.input_file_name = input_file_name
        self.choice = choice
        self.params = params
        self.start_time = None
        
    def run(self):
        self.start_time = datetime.now()
        input_file_name = self.input_file_name
        base_name = Path(input_file_name).with_suffix('').name
        aligned_file = base_name + ".aligned.fasta"
        trimmed_file = base_name + ".trimmed.fasta"
        
        results = {
            'success': False,
            'output_files': [],
            'execution_time': 0,
            'alignment_stats': {},
            'tree_file': None
        }
        
        try:
            # ==================================================================
            # STEP 1: ALIGNMENT
            # ==================================================================
            self.stage_update.emit("Alignment (MAFFT)")
            self.progress_update.emit(10)
            self.log_message.emit(f"\n [STEP 1/4] Aligning '{input_file_name}' with MAFFT...", 'default')
            
            mafft_strategy = self.params.get('mafft_strategy', 'auto')
            mafft_cmd = ['sh', '-c', f"mafft --{mafft_strategy} /data/{input_file_name} > /data/{aligned_file}"]
            run_docker_command("mafft", mafft_cmd, self.log_message, self.progress_increment)
            self.progress_update.emit(35)
            
            # Analyze alignment
            alignment_stats = analyze_alignment(aligned_file)
            if alignment_stats:
                self.log_message.emit(f"\n Alignment Quality:", 'highlight')
                self.log_message.emit(f"  • Length: {alignment_stats['alignment_length']} bp", 'default')
                self.log_message.emit(f"  • Gap %: {alignment_stats['gap_percentage']:.2f}%", 'default')
                self.log_message.emit(f"  • Avg Identity: {alignment_stats['avg_identity']:.2f}%", 'default')
                results['alignment_stats'] = alignment_stats
            
            # ==================================================================
            # STEP 2: TRIMMING
            # ==================================================================
            self.stage_update.emit("Trimming (trimAl)")
            self.progress_update.emit(40)
            self.log_message.emit(f"\n [STEP 2/4] Filtering alignment with trimAl...", 'default')
            
            trimal_method = self.params.get('trimal_method', 'automated1')
            trimal_cmd = ['trimal', '-in', f'/data/{aligned_file}', '-out', f'/data/{trimmed_file}', f'-{trimal_method}']
            run_docker_command("trimal", trimal_cmd, self.log_message, self.progress_increment)
            self.progress_update.emit(65)
            
            # ==================================================================
            # STEP 3: TREE CONSTRUCTION
            # ==================================================================
            self.stage_update.emit("Tree Construction")
            self.progress_update.emit(70)
            self.log_message.emit(f"\n [STEP 3/4] Building phylogenetic tree (Method {self.choice})...", 'default')
            
            if self.choice == '1':
                # Neighbor-Joining
                nj_output_file = base_name + ".nj_tree.newick"
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
print('\\n--- Tree Structure ---')
Phylo.draw_ascii(tree)
"""
                nj_script_name = "_temp_nj_script.py"
                with open(nj_script_name, "w") as f:
                    f.write(nj_script_content)
                
                biopython_cmd = ['python', f'/data/{nj_script_name}']
                run_docker_command("biopython", biopython_cmd, self.log_message, self.progress_increment)
                os.remove(nj_script_name)
                results['output_files'].append(nj_output_file)
                results['tree_file'] = nj_output_file
                
            elif self.choice == '2':
                # Maximum Likelihood
                bootstrap = self.params.get('iqtree_bootstrap', 1000)
                model = self.params.get('iqtree_model', 'MFP')
                
                iqtree_cmd = ['iqtree', '-s', f'/data/{trimmed_file}', '-m', model, '-B', str(bootstrap)]
                run_docker_command("iqtree", iqtree_cmd, self.log_message, self.progress_increment)
                tree_file = f"{trimmed_file}.treefile"
                results['output_files'].append(tree_file)
                results['tree_file'] = tree_file
                
            elif self.choice == '3':
                # Bayesian Inference
                nexus_file = base_name + ".trimmed.nex"
                
                biopython_convert_cmd = [
                    'python', '-c',
                    f"from Bio import AlignIO; AlignIO.convert('/data/{trimmed_file}', 'fasta', '/data/{nexus_file}', 'nexus', molecule_type='DNA')"
                ]
                run_docker_command("biopython", biopython_convert_cmd, self.log_message)
                
                ngen = self.params.get('mrbayes_ngen', 100000)
                nchains = self.params.get('mrbayes_nchains', 4)
                
                mrbayes_commands = f"""
begin mrbayes;
    set autoclose=yes nowarn=yes;
    lset nst=6 rates=gamma;
    mcmc ngen={ngen} samplefreq=100 nchains={nchains} nruns=2;
    sump burnin=250;
    sumt burnin=250;
    quit;
end;
"""
                with open(nexus_file, 'a') as f:
                    f.write(mrbayes_commands)
                
                mrbayes_cmd = ['mb', f'/data/{nexus_file}']
                run_docker_command("mrbayes", mrbayes_cmd, self.log_message, self.progress_increment)
                tree_file = f"{nexus_file}.con.tre"
                results['output_files'].append(tree_file)
                results['tree_file'] = tree_file

            # ==================================================================
            # STEP 4: PLOTTING (NOVO)
            # ==================================================================
            if results['tree_file']:
                self.stage_update.emit("Generating Plots (R)")
                self.progress_update.emit(90)
                self.log_message.emit(f"\n [STEP 4/4] Generating tree visualizations with R...", 'highlight')
                
                # Executa o container R passando o arquivo da árvore
                plot_cmd = [results['tree_file']]
                run_docker_command("cladegene_plot", plot_cmd, self.log_message)
                
                # Verifica visualmente se os plots foram criados para informar no log
                plots_dir = Path("plots")
                if plots_dir.exists():
                    plots_count = len(list(plots_dir.glob("tree_*.png")))
                    if plots_count > 0:
                        self.log_message.emit(f" Successfully generated {plots_count} plots in ./plots/", 'success')
            
            # Finalização
            self.progress_update.emit(100)
            
            # Calculate execution time
            execution_time = (datetime.now() - self.start_time).total_seconds()
            results['execution_time'] = execution_time
            results['success'] = True
            
            self.log_message.emit(f"\n  Total execution time: {execution_time:.2f} seconds", 'success')
            self.finished.emit(True, results)
            
        except Exception as e:
            self.log_message.emit(f"\n Pipeline failed: {str(e)}", 'error')
            self.finished.emit(False, results)

# ==============================================================================
# --- 3. MAIN APPLICATION ---
# ==============================================================================

class PhylogenyApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CladoGene 🧬")
        self.setGeometry(100, 100, 1200, 800)
        
        self.thread = None
        self.worker = None
        self.history = self.load_history()
        
        self.setup_ui()
        self.apply_stylesheet()
        
    def load_history(self):
        """Load execution history from file."""
        if Path(HISTORY_FILE).exists():
            try:
                with open(HISTORY_FILE, 'r') as f:
                    return json.load(f)
            except:
                return []
        return []
    
    def save_history(self, entry):
        """Save execution to history."""
        self.history.append(entry)
        with open(HISTORY_FILE, 'w') as f:
            json.dump(self.history[-50:], f, indent=2)  # Keep last 50
    
    def apply_stylesheet(self):
        """Enhanced dark blue/white theme."""
        stylesheet = """
        QWidget {
            background-color: #FFFFFF;
            color: #2C3E50;
            font-family: 'Segoe UI', Arial, sans-serif;
            font-size: 10pt;
        }
        QGroupBox {
            background-color: #F8F9FA;
            border: 2px solid #1F456E;
            border-radius: 8px;
            margin-top: 12px;
            padding: 15px;
            font-weight: bold;
            color: #1F456E;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            subcontrol-position: top left;
            padding: 0 5px;
            font-size: 11pt;
        }
        QPushButton {
            background-color: #1F456E;
            color: white;
            border: none;
            border-radius: 5px;
            padding: 10px 20px;
            font-weight: bold;
        }
        QPushButton:hover {
            background-color: #2A68A6;
        }
        QPushButton:disabled {
            background-color: #CCCCCC;
            color: #666666;
        }
        QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {
            background-color: white;
            border: 1px solid #BDC3C7;
            border-radius: 4px;
            padding: 6px;
        }
        QTextEdit {
            background-color: #F8F9FA;
            border: 1px solid #BDC3C7;
            border-radius: 4px;
            font-family: 'Courier New', monospace;
            font-size: 9pt;
        }
        QProgressBar {
            border: 1px solid #BDC3C7;
            border-radius: 5px;
            text-align: center;
            background-color: #ECF0F1;
        }
        QProgressBar::chunk {
            background-color: #27AE60;
            border-radius: 4px;
        }
        QTabWidget::pane {
            border: 1px solid #BDC3C7;
            border-radius: 4px;
            background-color: white;
        }
        QTabBar::tab {
            background-color: #ECF0F1;
            color: #2C3E50;
            padding: 8px 16px;
            margin-right: 2px;
            border-top-left-radius: 4px;
            border-top-right-radius: 4px;
        }
        QTabBar::tab:selected {
            background-color: #1F456E;
            color: white;
        }
        """
        self.setStyleSheet(stylesheet)
    
    def setup_ui(self):
        main_layout = QVBoxLayout(self)
        
        # Title
        title = QLabel("🧬 CladeGene - Phylogenetic Analysis")
        title.setFont(QFont('Arial', 16, QFont.Bold))
        title.setStyleSheet("color: #1F456E; padding: 10px;")
        main_layout.addWidget(title)
        
        # Tabs
        tabs = QTabWidget()
        tabs.addTab(self.create_pipeline_tab(), " Pipeline")
        tabs.addTab(self.create_parameters_tab(), "⚙️ Parameters")
        tabs.addTab(self.create_history_tab(), " History")
        main_layout.addWidget(tabs)
        
    def create_pipeline_tab(self):
        """Main pipeline tab with split layout."""
        tab = QWidget()
        main_layout = QHBoxLayout(tab)
        
        # Create splitter for left/right division
        splitter = QSplitter(Qt.Horizontal)
        
        # ========== LEFT SIDE: Controls ==========
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 5, 0)
        
        # Input Section
        input_group = QGroupBox(" Input File")
        input_layout = QVBoxLayout(input_group)
        
        file_layout = QHBoxLayout()
        self.input_line = QLineEdit()
        self.input_line.setReadOnly(True)
        self.input_line.setPlaceholderText("Select a FASTA file...")
        file_layout.addWidget(QLabel("FASTA:"))
        file_layout.addWidget(self.input_line)
        
        browse_btn = QPushButton(" Browse")
        browse_btn.clicked.connect(self.select_file)
        file_layout.addWidget(browse_btn)
        input_layout.addLayout(file_layout)
        
        # Preview area
        self.preview_text = QTextEdit()
        self.preview_text.setReadOnly(True)
        self.preview_text.setMaximumHeight(100)
        self.preview_text.setPlaceholderText("File statistics will appear here...")
        input_layout.addWidget(QLabel(" File Preview:"))
        input_layout.addWidget(self.preview_text)
        
        left_layout.addWidget(input_group)
        
        # Method Selection
        method_group = QGroupBox(" Phylogenetic Method")
        method_layout = QVBoxLayout(method_group)
        
        self.radio_group = QButtonGroup()
        self.radio_nj = QRadioButton(" Neighbor-Joining (Fast, distance-based)")
        self.radio_ml = QRadioButton(" Maximum Likelihood (Accurate, model-based)")
        self.radio_bi = QRadioButton(" Bayesian Inference (Probabilistic, posterior)")
        
        self.radio_ml.setChecked(True)
        
        self.radio_group.addButton(self.radio_nj, 1)
        self.radio_group.addButton(self.radio_ml, 2)
        self.radio_group.addButton(self.radio_bi, 3)
        
        method_layout.addWidget(self.radio_nj)
        method_layout.addWidget(self.radio_ml)
        method_layout.addWidget(self.radio_bi)
        
        left_layout.addWidget(method_group)
        
        # Control Section
        control_layout = QHBoxLayout()
        
        self.start_button = QPushButton(" RUN ANALYSIS")
        self.start_button.setMinimumHeight(50)
        self.start_button.setFont(QFont('Arial', 12, QFont.Bold))
        self.start_button.clicked.connect(self.start_pipeline)
        control_layout.addWidget(self.start_button)
        
        self.stop_button = QPushButton(" STOP")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self.stop_pipeline)
        control_layout.addWidget(self.stop_button)
        
        left_layout.addLayout(control_layout)
        
        # Progress Section
        progress_group = QGroupBox(" Progress")
        progress_layout = QVBoxLayout(progress_group)
        
        self.stage_label = QLabel("Ready to start...")
        self.stage_label.setFont(QFont('Arial', 10, QFont.Bold))
        progress_layout.addWidget(self.stage_label)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        progress_layout.addWidget(self.progress_bar)
        
        left_layout.addWidget(progress_group)
        left_layout.addStretch()  # Push logo to bottom
        
        # Logo Section at bottom - CladeGene logo
        logo_frame = QFrame()
        logo_frame.setStyleSheet("""
            QFrame {
                background-color: white;
                border: 1px solid #E0E0E0;
                border-radius: 5px;
                padding: 15px;
                min-height: 120px;
            }
        """)
        logo_layout = QVBoxLayout(logo_frame)
        logo_layout.setAlignment(Qt.AlignCenter)
        logo_layout.setContentsMargins(10, 10, 10, 10)
        
        # Load CladeGene logo from figures folder
        logo_path = Path(__file__).parent / "figures" / "CladeGene.png"
        
        from PySide6.QtGui import QPixmap
        if logo_path.exists():
            logo_label = QLabel()
            pixmap = QPixmap(str(logo_path))
            # Scale to fit the frame width (approximately 80% of left panel)
            target_width = 280
            scaled_pixmap = pixmap.scaledToWidth(target_width, Qt.SmoothTransformation)
            logo_label.setPixmap(scaled_pixmap)
            logo_label.setAlignment(Qt.AlignCenter)
            logo_layout.addWidget(logo_label)
        else:
            # Fallback if logo not found
            fallback = QLabel("🧬 CladeGene")
            fallback.setAlignment(Qt.AlignCenter)
            fallback.setStyleSheet("color: #1F456E; font-size: 16pt; font-weight: bold;")
            logo_layout.addWidget(fallback)
        
        left_layout.addWidget(logo_frame)
        
        # ========== RIGHT SIDE: Log ==========
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(5, 0, 0, 0)
        
        log_group = QGroupBox(" Execution Log")
        log_layout = QVBoxLayout(log_group)
        
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.set_log_colors()
        log_layout.addWidget(self.log_text)
        
        # Removido: Botão de plotagem manual
        
        right_layout.addWidget(log_group)
        
        # Add both sides to splitter
        splitter.addWidget(left_widget)
        splitter.addWidget(right_widget)
        
        # Set initial sizes (40% left, 60% right)
        splitter.setSizes([400, 600])
        splitter.setStretchFactor(0, 2)  # Left side less stretchable
        splitter.setStretchFactor(1, 3)  # Right side more stretchable
        
        main_layout.addWidget(splitter)
        
        return tab
    
    def create_parameters_tab(self):
        """Parameters configuration tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # MAFFT Parameters
        mafft_group = QGroupBox("MAFFT (Alignment)")
        mafft_layout = QVBoxLayout(mafft_group)
        
        mafft_h = QHBoxLayout()
        mafft_h.addWidget(QLabel("Strategy:"))
        self.mafft_strategy = QComboBox()
        self.mafft_strategy.addItems(['auto', 'localpair', 'globalpair', 'retree 2', 'fftns'])
        mafft_h.addWidget(self.mafft_strategy)
        mafft_layout.addLayout(mafft_h)
        
        scroll_layout.addWidget(mafft_group)
        
        # trimAl Parameters
        trimal_group = QGroupBox("trimAl (Filtering)")
        trimal_layout = QVBoxLayout(trimal_group)
        
        trimal_h = QHBoxLayout()
        trimal_h.addWidget(QLabel("Method:"))
        self.trimal_method = QComboBox()
        self.trimal_method.addItems(['automated1', 'strict', 'gappyout'])
        trimal_h.addWidget(self.trimal_method)
        trimal_layout.addLayout(trimal_h)
        
        scroll_layout.addWidget(trimal_group)
        
        # IQ-TREE Parameters
        iqtree_group = QGroupBox("IQ-TREE (Maximum Likelihood)")
        iqtree_layout = QVBoxLayout(iqtree_group)
        
        iqtree_model_h = QHBoxLayout()
        iqtree_model_h.addWidget(QLabel("Model:"))
        self.iqtree_model = QComboBox()
        self.iqtree_model.addItems(['MFP', 'GTR+G', 'HKY+G', 'JC'])
        iqtree_model_h.addWidget(self.iqtree_model)
        iqtree_layout.addLayout(iqtree_model_h)
        
        iqtree_boot_h = QHBoxLayout()
        iqtree_boot_h.addWidget(QLabel("Bootstrap:"))
        self.iqtree_bootstrap = QSpinBox()
        self.iqtree_bootstrap.setRange(100, 10000)
        self.iqtree_bootstrap.setValue(1000)
        self.iqtree_bootstrap.setSingleStep(100)
        iqtree_boot_h.addWidget(self.iqtree_bootstrap)
        iqtree_layout.addLayout(iqtree_boot_h)
        
        scroll_layout.addWidget(iqtree_group)
        
        # MrBayes Parameters
        mrbayes_group = QGroupBox("MrBayes (Bayesian Inference)")
        mrbayes_layout = QVBoxLayout(mrbayes_group)
        
        mb_ngen_h = QHBoxLayout()
        mb_ngen_h.addWidget(QLabel("Generations:"))
        self.mrbayes_ngen = QSpinBox()
        self.mrbayes_ngen.setRange(10000, 10000000)
        self.mrbayes_ngen.setValue(100000)
        self.mrbayes_ngen.setSingleStep(10000)
        mb_ngen_h.addWidget(self.mrbayes_ngen)
        mrbayes_layout.addLayout(mb_ngen_h)
        
        mb_chains_h = QHBoxLayout()
        mb_chains_h.addWidget(QLabel("Chains:"))
        self.mrbayes_nchains = QSpinBox()
        self.mrbayes_nchains.setRange(2, 20)
        self.mrbayes_nchains.setValue(4)
        mb_chains_h.addWidget(self.mrbayes_nchains)
        mrbayes_layout.addLayout(mb_chains_h)
        
        scroll_layout.addWidget(mrbayes_group)
        
        scroll_layout.addStretch()
        scroll.setWidget(scroll_widget)
        layout.addWidget(scroll)
        
        # Reset button
        reset_btn = QPushButton(" Reset to Defaults")
        reset_btn.clicked.connect(self.reset_parameters)
        layout.addWidget(reset_btn)
        
        return tab
    
    def create_history_tab(self):
        """Execution history tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        layout.addWidget(QLabel(" Recent Executions:"))
        
        self.history_table = QTableWidget()
        self.history_table.setColumnCount(5)
        self.history_table.setHorizontalHeaderLabels(['Date', 'File', 'Method', 'Time (s)', 'Status'])
        self.history_table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.history_table)
        
        self.update_history_table()
        
        clear_btn = QPushButton(" Clear History")
        clear_btn.clicked.connect(self.clear_history)
        layout.addWidget(clear_btn)
        
        return tab
    
    def update_history_table(self):
        """Update history table with data."""
        self.history_table.setRowCount(len(self.history))
        for i, entry in enumerate(reversed(self.history)):
            self.history_table.setItem(i, 0, QTableWidgetItem(entry.get('date', 'N/A')))
            self.history_table.setItem(i, 1, QTableWidgetItem(entry.get('file', 'N/A')))
            self.history_table.setItem(i, 2, QTableWidgetItem(entry.get('method', 'N/A')))
            self.history_table.setItem(i, 3, QTableWidgetItem(f"{entry.get('time', 0):.2f}"))
            self.history_table.setItem(i, 4, QTableWidgetItem('correct' if entry.get('success') else 'error'))
    
    def clear_history(self):
        """Clear execution history."""
        reply = QMessageBox.question(self, 'Confirm', 'Clear all history?',
                                      QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.history = []
            if Path(HISTORY_FILE).exists():
                os.remove(HISTORY_FILE)
            self.update_history_table()
    
    def set_log_colors(self):
        """Define log text formats."""
        self.log_text.setFont(QFont("Courier New", 9))
        
        self.fmt_error = QTextCharFormat()
        self.fmt_error.setForeground(QColor("#E74C3C"))
        self.fmt_error.setFontWeight(QFont.Bold)
        
        self.fmt_success = QTextCharFormat()
        self.fmt_success.setForeground(QColor("#27AE60"))
        self.fmt_success.setFontWeight(QFont.Bold)
        
        self.fmt_highlight = QTextCharFormat()
        self.fmt_highlight.setForeground(QColor("#3498DB"))
        self.fmt_highlight.setFontWeight(QFont.Bold)
    
    def reset_parameters(self):
        """Reset all parameters to defaults."""
        self.mafft_strategy.setCurrentText('auto')
        self.trimal_method.setCurrentText('automated1')
        self.iqtree_model.setCurrentText('MFP')
        self.iqtree_bootstrap.setValue(1000)
        self.mrbayes_ngen.setValue(100000)
        self.mrbayes_nchains.setValue(4)
        QMessageBox.information(self, "Reset", "Parameters reset to defaults!")
    
    def select_file(self):
        """File selection with validation."""
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Select Multi-FASTA File", "", 
            "FASTA Files (*.fasta *.fa *.fna);;All Files (*)"
        )
        
        if filepath:
            os.chdir(Path(filepath).parent)
            self.input_line.setText(Path(filepath).name)
            self.log(f" File selected: {Path(filepath).name}", 'highlight')
            
            # Validate and show preview
            valid, issues, stats = validate_fasta(filepath)
            
            preview_text = ""
            if valid:
                preview_text = f" Valid FASTA file\n"
                preview_text += f" Sequences: {stats['num_sequences']}\n"
                preview_text += f" Length range: {stats['min_length']}-{stats['max_length']} bp\n"
                preview_text += f" Average length: {stats['avg_length']:.1f} bp\n"
                preview_text += f" Total size: {stats['total_length']:,} bp"
            else:
                preview_text = " VALIDATION ISSUES:\n"
                for issue in issues:
                    preview_text += f"  • {issue}\n"
            
            self.preview_text.setPlainText(preview_text)
            
            if not valid:
                QMessageBox.warning(self, "Validation Warning", 
                                   f"File has {len(issues)} issue(s). Check preview.")
    
    def log(self, message, message_type='default'):
        """Add message to log with formatting."""
        cursor = self.log_text.textCursor()
        cursor.movePosition(QTextCursor.End)
        
        if message_type == 'error':
            fmt = self.fmt_error
        elif message_type == 'success':
            fmt = self.fmt_success
        elif message_type == 'highlight':
            fmt = self.fmt_highlight
        else:
            fmt = QTextCharFormat()
        
        cursor.insertText(message + "\n", fmt)
        self.log_text.setTextCursor(cursor)
        self.log_text.ensureCursorVisible()
    
    def get_parameters(self):
        """Get current parameter values."""
        return {
            'mafft_strategy': self.mafft_strategy.currentText(),
            'trimal_method': self.trimal_method.currentText(),
            'iqtree_model': self.iqtree_model.currentText(),
            'iqtree_bootstrap': self.iqtree_bootstrap.value(),
            'mrbayes_ngen': self.mrbayes_ngen.value(),
            'mrbayes_nchains': self.mrbayes_nchains.value()
        }
    
    def toggle_interface(self, enable):
        """Enable/disable interface elements."""
        self.start_button.setEnabled(enable)
        self.stop_button.setEnabled(not enable)
        self.radio_nj.setEnabled(enable)
        self.radio_ml.setEnabled(enable)
        self.radio_bi.setEnabled(enable)
        self.input_line.setReadOnly(not enable)
        
        browse_btn = self.findChild(QPushButton, "📁 Browse")
        if browse_btn:
            browse_btn.setEnabled(enable)
    
    def stop_pipeline(self):
        """Stop running pipeline."""
        if self.thread and self.thread.isRunning():
            reply = QMessageBox.question(
                self, 'Confirm Stop',
                'Are you sure you want to stop the pipeline?',
                QMessageBox.Yes | QMessageBox.No
            )
            
            if reply == QMessageBox.Yes:
                self.log(" Stopping pipeline...", 'error')
                self.thread.terminate()
                self.thread.wait()
                self.toggle_interface(True)
                self.progress_bar.setValue(0)
                self.stage_label.setText("Stopped by user")
    
    def start_pipeline(self):
        #Start the phylogenetic pipeline.
        input_file_name = self.input_line.text()
        
        if not input_file_name:
            QMessageBox.warning(self, "Input Required", 
                              "Please select a FASTA input file.")
            return
        
        if not Path(input_file_name).exists():
            QMessageBox.critical(self, "File Not Found",
                               f"File '{input_file_name}' not found!")
            return
        
        choice = str(self.radio_group.checkedId())
        params = self.get_parameters()
        
        method_names = {
            '1': 'Neighbor-Joining',
            '2': 'Maximum Likelihood',
            '3': 'Bayesian Inference'
        }
        
        self.log_text.clear()
        self.log("="*60, 'default')
        self.log("🧬 CladoGene - Analysis Started", 'highlight')
        self.log("="*60, 'default')
        self.log(f"📁 Input: {input_file_name}", 'default')
        self.log(f" Method: {method_names[choice]}", 'default')
        self.log(f"  Parameters: {params}", 'default')
        self.log(f" Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", 'default')
        self.log("="*60, 'default')
        
        self.toggle_interface(False)
        self.progress_bar.setValue(0)
        self.stage_label.setText("Initializing...")
        
        # Create worker thread
        self.thread = QThread()
        self.worker = PhylogenyWorker(input_file_name, choice, params)
        self.worker.moveToThread(self.thread)
        
        # Connect signals
        self.thread.started.connect(self.worker.run)
        self.worker.log_message.connect(self.log)
        self.worker.progress_update.connect(self.progress_bar.setValue)
        self.worker.stage_update.connect(self.stage_label.setText)
        self.worker.finished.connect(self.pipeline_finished)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        
        self.thread.start()
    
    def pipeline_finished(self, success, results):
        """Handle pipeline completion."""
        self.toggle_interface(True)
        self.thread.quit()
        self.thread.wait()
        
        if success:
            self.log("="*60, 'success')
            self.log("✅ ANALYSIS COMPLETED SUCCESSFULLY!", 'success')
            self.log("="*60, 'success')
            
            if results.get('alignment_stats'):
                stats = results['alignment_stats']
                self.log(f"\n Final Alignment Statistics:", 'highlight')
                self.log(f"  • Length: {stats.get('alignment_length', 'N/A')} bp", 'default')
                self.log(f"  • Gap percentage: {stats.get('gap_percentage', 0):.2f}%", 'default')
                self.log(f"  • Average identity: {stats.get('avg_identity', 0):.2f}%", 'default')
            
            if results.get('output_files'):
                self.log(f"\n Output files:", 'highlight')
                for f in results['output_files']:
                    self.log(f"  • {f}", 'default')
            
            self.log(f"\n Total time: {results.get('execution_time', 0):.2f} seconds", 'success')
            
            # Save to history
            history_entry = {
                'date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'file': self.input_line.text(),
                'method': ['NJ', 'ML', 'BI'][int(str(self.radio_group.checkedId()))-1],
                'time': results.get('execution_time', 0),
                'success': True
            }
            self.save_history(history_entry)
            self.update_history_table()
            
            QMessageBox.information(self, "Success ", 
                                  f"Phylogenetic analysis completed!\n\n"
                                  f"Time: {results.get('execution_time', 0):.2f}s\n"
                                  f"Output: {len(results.get('output_files', []))} file(s)")
        else:
            self.log("="*60, 'error')
            self.log(" ANALYSIS FAILED", 'error')
            self.log("="*60, 'error')
            
            history_entry = {
                'date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'file': self.input_line.text(),
                'method': ['NJ', 'ML', 'BI'][int(str(self.radio_group.checkedId()))-1],
                'time': results.get('execution_time', 0),
                'success': False
            }
            self.save_history(history_entry)
            self.update_history_table()
            
            QMessageBox.critical(self, "Error ",
                               "Pipeline failed. Check the log for details.")
        
        self.progress_bar.setValue(100 if success else 0)
        self.stage_label.setText(" Completed" if success else " Failed")

# ==============================================================================
# --- 4. MAIN EXECUTION ---
# ==============================================================================

if __name__ == "__main__":
    os.chdir(Path(__file__).parent)
    
    app = QApplication(sys.argv)
    app.setStyle('Fusion')  # Modern cross-platform style
    
    window = PhylogenyApp()
    window.show()
    
    sys.exit(app.exec())