#!/usr/bin/env python3
"""
Phylo-Pipe Pro - Command Line Interface
Advanced phylogenetic analysis pipeline with Docker containers
"""

import os
import sys
import argparse
import subprocess
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple

# ==============================================================================
# CONFIGURATION
# ==============================================================================

CONTAINERS = {
    "mafft": "quay.io/biocontainers/mafft:7.525--h031d066_0",
    "trimal": "quay.io/biocontainers/trimal:1.5",
    "iqtree": "quay.io/biocontainers/iqtree:3.0.1--h503566f_0",
    "biopython": "quay.io/biocontainers/biopython:1.79",
    "mrbayes": "quay.io/biocontainers/mrbayes:3.2.7--hd0d793b_7"
}

# Color codes for terminal output
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

# ==============================================================================
# VALIDATION FUNCTIONS
# ==============================================================================

def validate_fasta(filepath: Path) -> Tuple[bool, List[str], Dict]:
    """
    Validate FASTA file format and content.
    Returns: (is_valid, issues_list, statistics_dict)
    """
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
                    if not current_id:
                        issues.append(f"Line {line_num}: Empty sequence ID")
                    elif current_id in sequences:
                        issues.append(f"Line {line_num}: Duplicate ID '{current_id}'")
                    current_seq = []
                else:
                    # Validate characters
                    invalid_chars = set(line.upper()) - set('ACGTRYSWKMBDHVN-')
                    if invalid_chars:
                        issues.append(
                            f"Line {line_num}: Invalid characters {invalid_chars} "
                            f"(only IUPAC DNA codes allowed)"
                        )
                    current_seq.append(line.upper())
            
            # Add last sequence
            if current_id:
                sequences[current_id] = ''.join(current_seq)
        
        if not sequences:
            issues.append("No sequences found in file")
            return False, issues, {}
        
        if len(sequences) < 3:
            issues.append(f"Only {len(sequences)} sequences found (minimum 3 recommended)")
        
        # Calculate statistics
        seq_lengths = [len(seq.replace('-', '')) for seq in sequences.values()]
        
        stats = {
            'num_sequences': len(sequences),
            'min_length': min(seq_lengths),
            'max_length': max(seq_lengths),
            'avg_length': sum(seq_lengths) / len(seq_lengths),
            'total_length': sum(seq_lengths),
            'sequences': sequences
        }
        
        # Check for very different lengths
        if stats['max_length'] > stats['min_length'] * 3:
            issues.append(
                f"WARNING: Large length variation "
                f"({stats['min_length']}-{stats['max_length']} bp)"
            )
        
        return len([i for i in issues if not i.startswith('WARNING')]) == 0, issues, stats
        
    except Exception as e:
        return False, [f"Error reading file: {str(e)}"], {}

def analyze_alignment(filepath: Path) -> Dict:
    """Analyze alignment quality metrics."""
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
                    current_seq.append(line.upper())
            if current_id:
                sequences[current_id] = ''.join(current_seq)
        
        if not sequences:
            return {}
        
        alignment_length = len(list(sequences.values())[0])
        total_gaps = sum(seq.count('-') for seq in sequences.values())
        total_positions = len(sequences) * alignment_length
        
        # Pairwise identity
        identity_sum = 0
        comparisons = 0
        seq_list = list(sequences.values())
        
        for i in range(len(seq_list)):
            for j in range(i+1, len(seq_list)):
                matches = sum(
                    1 for a, b in zip(seq_list[i], seq_list[j])
                    if a == b and a != '-'
                )
                identity_sum += matches / alignment_length
                comparisons += 1
        
        avg_identity = (identity_sum / comparisons * 100) if comparisons > 0 else 0
        
        # Calculate conserved positions
        conserved = 0
        for pos in range(alignment_length):
            column = [seq[pos] for seq in seq_list if seq[pos] != '-']
            if column and len(set(column)) == 1:
                conserved += 1
        
        return {
            'alignment_length': alignment_length,
            'gap_percentage': (total_gaps / total_positions * 100),
            'avg_identity': avg_identity,
            'conserved_positions': conserved,
            'conserved_percentage': (conserved / alignment_length * 100)
        }
    except Exception as e:
        print(f"{Colors.YELLOW}Warning: Could not analyze alignment: {e}{Colors.END}")
        return {}

# ==============================================================================
# DOCKER EXECUTION
# ==============================================================================

def run_docker_command(container_name: str, command: List[str], verbose: bool = True) -> bool:
    """Execute Docker command with real-time output."""
    work_dir = Path.cwd().resolve()
    
    docker_base = [
        "docker", "run", "--rm",
        "-v", f"{work_dir}:/data",
        "-w", "/data"
    ]
    
    full_command = docker_base + [CONTAINERS[container_name]] + command
    
    if verbose:
        print(f"\n{Colors.CYAN}🔧 Executing:{Colors.END} {' '.join(full_command)}")
    
    try:
        process = subprocess.Popen(
            full_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1
        )
        
        # Stream output in real-time
        stdout_lines = []
        stderr_lines = []
        
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                stdout_lines.append(output.strip())
                if verbose:
                    print(output.strip())
        
        stderr_output = process.stderr.read()
        if stderr_output:
            stderr_lines.append(stderr_output.strip())
        
        returncode = process.poll()
        
        if returncode != 0:
            print(f"\n{Colors.RED} ERROR: Command failed with exit code {returncode}{Colors.END}")
            if stderr_output and verbose:
                print(f"{Colors.RED}STDERR:\n{stderr_output}{Colors.END}")
            return False
        
        if verbose:
            print(f"{Colors.GREEN} Command completed successfully!{Colors.END}")
        return True
        
    except Exception as e:
        print(f"{Colors.RED} Docker execution error: {e}{Colors.END}")
        return False

# ==============================================================================
# PIPELINE STAGES
# ==============================================================================

def stage_alignment(input_file: Path, strategy: str = 'auto', verbose: bool = True) -> Path:
    """Stage 1: Multiple sequence alignment with MAFFT."""
    aligned_file = input_file.with_suffix('').with_suffix('.aligned.fasta')
    
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.HEADER}{Colors.BOLD}STAGE 1: MULTIPLE SEQUENCE ALIGNMENT{Colors.END}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.BLUE}Tool:{Colors.END} MAFFT")
    print(f"{Colors.BLUE}Strategy:{Colors.END} {strategy}")
    print(f"{Colors.BLUE}Input:{Colors.END} {input_file}")
    print(f"{Colors.BLUE}Output:{Colors.END} {aligned_file}")
    
    cmd = ['sh', '-c', f"mafft --{strategy} /data/{input_file.name} > /data/{aligned_file.name}"]
    
    if not run_docker_command("mafft", cmd, verbose):
        raise RuntimeError("MAFFT alignment failed")
    
    # Analyze alignment
    stats = analyze_alignment(aligned_file)
    if stats:
        print(f"\n{Colors.CYAN} Alignment Quality:{Colors.END}")
        print(f"  • Length: {stats['alignment_length']} bp")
        print(f"  • Gap percentage: {stats['gap_percentage']:.2f}%")
        print(f"  • Average identity: {stats['avg_identity']:.2f}%")
        print(f"  • Conserved positions: {stats['conserved_positions']} ({stats['conserved_percentage']:.2f}%)")
    
    return aligned_file

def stage_trimming(aligned_file: Path, method: str = 'automated1', verbose: bool = True) -> Path:
    """Stage 2: Alignment trimming with trimAl."""
    trimmed_file = aligned_file.with_suffix('').with_suffix('.trimmed.fasta')
    
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.HEADER}{Colors.BOLD}STAGE 2: ALIGNMENT TRIMMING{Colors.END}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.BLUE}Tool:{Colors.END} trimAl")
    print(f"{Colors.BLUE}Method:{Colors.END} {method}")
    print(f"{Colors.BLUE}Input:{Colors.END} {aligned_file}")
    print(f"{Colors.BLUE}Output:{Colors.END} {trimmed_file}")
    
    cmd = [
        'trimal',
        '-in', f'/data/{aligned_file.name}',
        '-out', f'/data/{trimmed_file.name}',
        f'-{method}'
    ]
    
    if not run_docker_command("trimal", cmd, verbose):
        raise RuntimeError("trimAl filtering failed")
    
    # Analyze trimmed alignment
    stats = analyze_alignment(trimmed_file)
    if stats:
        print(f"\n{Colors.CYAN} Trimmed Alignment Quality:{Colors.END}")
        print(f"  • Length: {stats['alignment_length']} bp")
        print(f"  • Gap percentage: {stats['gap_percentage']:.2f}%")
        print(f"  • Average identity: {stats['avg_identity']:.2f}%")
    
    return trimmed_file

def stage_tree_nj(trimmed_file: Path, verbose: bool = True) -> Path:
    """Stage 3A: Neighbor-Joining tree construction."""
    output_file = trimmed_file.with_suffix('').with_suffix('.nj_tree.newick')
    
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.HEADER}{Colors.BOLD}STAGE 3: TREE CONSTRUCTION - NEIGHBOR-JOINING{Colors.END}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.BLUE}Method:{Colors.END} Distance-based (Neighbor-Joining)")
    print(f"{Colors.BLUE}Input:{Colors.END} {trimmed_file}")
    print(f"{Colors.BLUE}Output:{Colors.END} {output_file}")
    
    script_content = f"""
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import sys
sys.dont_write_bytecode = True

print('Reading alignment...')
aln = AlignIO.read('/data/{trimmed_file.name}', 'fasta')

print('Calculating distance matrix...')
calculator = DistanceCalculator('identity')
constructor = DistanceTreeConstructor(calculator, 'nj')

print('Building NJ tree...')
tree = constructor.build_tree(aln)
tree.ladderize()

print('Saving tree...')
Phylo.write(tree, '/data/{output_file.name}', 'newick')

print('\\n--- Tree Visualization ---')
Phylo.draw_ascii(tree)
"""
    
    script_file = Path("_temp_nj_script.py")
    script_file.write_text(script_content)
    
    try:
        cmd = ['python', f'/data/{script_file.name}']
        if not run_docker_command("biopython", cmd, verbose):
            raise RuntimeError("NJ tree construction failed")
    finally:
        script_file.unlink(missing_ok=True)
    
    print(f"\n{Colors.GREEN} NJ tree saved to: {output_file}{Colors.END}")
    return output_file

def stage_tree_ml(trimmed_file: Path, model: str = 'MFP', bootstrap: int = 1000, 
                  verbose: bool = True) -> Path:
    """Stage 3B: Maximum Likelihood tree with IQ-TREE."""
    output_file = trimmed_file.with_suffix('.treefile')
    
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.HEADER}{Colors.BOLD}STAGE 3: TREE CONSTRUCTION - MAXIMUM LIKELIHOOD{Colors.END}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.BLUE}Tool:{Colors.END} IQ-TREE")
    print(f"{Colors.BLUE}Model:{Colors.END} {model} (ModelFinder if MFP)")
    print(f"{Colors.BLUE}Bootstrap:{Colors.END} {bootstrap} replicates")
    print(f"{Colors.BLUE}Input:{Colors.END} {trimmed_file}")
    print(f"{Colors.BLUE}Output:{Colors.END} {output_file}")
    
    cmd = [
        'iqtree',
        '-s', f'/data/{trimmed_file.name}',
        '-m', model,
        '-B', str(bootstrap)
    ]
    
    if not run_docker_command("iqtree", cmd, verbose):
        raise RuntimeError("IQ-TREE analysis failed")
    
    print(f"\n{Colors.GREEN} ML tree saved to: {output_file}{Colors.END}")
    return output_file

def stage_tree_bayesian(trimmed_file: Path, ngen: int = 100000, nchains: int = 4,
                       verbose: bool = True) -> Path:
    """Stage 3C: Bayesian Inference with MrBayes."""
    nexus_file = trimmed_file.with_suffix('.nex')
    output_file = Path(str(nexus_file) + '.con.tre')
    
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.HEADER}{Colors.BOLD}STAGE 3: TREE CONSTRUCTION - BAYESIAN INFERENCE{Colors.END}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.BLUE}Tool:{Colors.END} MrBayes")
    print(f"{Colors.BLUE}Model:{Colors.END} GTR+G")
    print(f"{Colors.BLUE}Generations:{Colors.END} {ngen:,}")
    print(f"{Colors.BLUE}Chains:{Colors.END} {nchains}")
    print(f"{Colors.BLUE}Input:{Colors.END} {trimmed_file}")
    print(f"{Colors.BLUE}Output:{Colors.END} {output_file}")
    
    # Convert FASTA to NEXUS
    print(f"\n{Colors.CYAN}Converting FASTA to NEXUS format...{Colors.END}")
    convert_cmd = [
        'python', '-c',
        f"from Bio import AlignIO; AlignIO.convert('/data/{trimmed_file.name}', 'fasta', "
        f"'/data/{nexus_file.name}', 'nexus', molecule_type='DNA')"
    ]
    
    if not run_docker_command("biopython", convert_cmd, verbose):
        raise RuntimeError("FASTA to NEXUS conversion failed")
    
    # Add MrBayes commands
    burnin = int(ngen * 0.25 / 100)  # 25% burnin
    
    mrbayes_block = f"""
begin mrbayes;
    set autoclose=yes nowarn=yes;
    
    [ Substitution model: GTR+Gamma ]
    lset nst=6 rates=gamma;
    
    [ MCMC parameters ]
    mcmc ngen={ngen} samplefreq=100 nchains={nchains} nruns=2;
    
    [ Summarize results (25% burnin) ]
    sump burnin={burnin};
    sumt burnin={burnin};
    
    quit;
end;
"""
    
    with open(nexus_file, 'a') as f:
        f.write(mrbayes_block)
    
    # Run MrBayes
    print(f"\n{Colors.CYAN}Running MrBayes (this may take a while)...{Colors.END}")
    cmd = ['mb', f'/data/{nexus_file.name}']
    
    if not run_docker_command("mrbayes", cmd, verbose):
        raise RuntimeError("MrBayes analysis failed")
    
    print(f"\n{Colors.GREEN} Bayesian consensus tree saved to: {output_file}{Colors.END}")
    return output_file

# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

def run_pipeline(input_file: Path, method: str, params: Dict, verbose: bool = True):
    """Run complete phylogenetic analysis pipeline."""
    
    start_time = datetime.now()
    
    print(f"\n{Colors.BOLD}{Colors.HEADER}{'='*70}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.HEADER} PHYLO-PIPE - PHYLOGENETIC ANALYSIS{Colors.END}")
    print(f"{Colors.BOLD}{Colors.HEADER}{'='*70}{Colors.END}")
    print(f"{Colors.BLUE}Input file:{Colors.END} {input_file}")
    print(f"{Colors.BLUE}Method:{Colors.END} {method}")
    print(f"{Colors.BLUE}Started:{Colors.END} {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{Colors.BOLD}{Colors.HEADER}{'='*70}{Colors.END}")
    
    # Validate input
    print(f"\n{Colors.CYAN} Validating input file...{Colors.END}")
    valid, issues, stats = validate_fasta(input_file)
    
    if issues:
        print(f"\n{Colors.YELLOW}  Validation issues found:{Colors.END}")
        for issue in issues:
            color = Colors.YELLOW if issue.startswith('WARNING') else Colors.RED
            print(f"  {color}• {issue}{Colors.END}")
    
    if not valid:
        print(f"\n{Colors.RED} Input file validation failed. Cannot proceed.{Colors.END}")
        sys.exit(1)
    
    print(f"\n{Colors.GREEN} Input validation passed{Colors.END}")
    print(f"{Colors.CYAN} Input Statistics:{Colors.END}")
    print(f"  • Sequences: {stats['num_sequences']}")
    print(f"  • Length range: {stats['min_length']}-{stats['max_length']} bp")
    print(f"  • Average length: {stats['avg_length']:.1f} bp")
    
    try:
        # Stage 1: Alignment
        aligned_file = stage_alignment(
            input_file,
            strategy=params.get('mafft_strategy', 'auto'),
            verbose=verbose
        )
        
        # Stage 2: Trimming
        trimmed_file = stage_trimming(
            aligned_file,
            method=params.get('trimal_method', 'automated1'),
            verbose=verbose
        )
        
        # Stage 3: Tree construction
        if method == 'nj':
            output_file = stage_tree_nj(trimmed_file, verbose)
        elif method == 'ml':
            output_file = stage_tree_ml(
                trimmed_file,
                model=params.get('iqtree_model', 'MFP'),
                bootstrap=params.get('iqtree_bootstrap', 1000),
                verbose=verbose
            )
        elif method == 'bayesian':
            output_file = stage_tree_bayesian(
                trimmed_file,
                ngen=params.get('mrbayes_ngen', 100000),
                nchains=params.get('mrbayes_nchains', 4),
                verbose=verbose
            )
        else:
            raise ValueError(f"Unknown method: {method}")
        
        # Summary
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        print(f"\n{Colors.BOLD}{Colors.GREEN}{'='*70}{Colors.END}")
        print(f"{Colors.BOLD}{Colors.GREEN} ANALYSIS COMPLETED SUCCESSFULLY!{Colors.END}")
        print(f"{Colors.BOLD}{Colors.GREEN}{'='*70}{Colors.END}")
        print(f"{Colors.BLUE}Final tree:{Colors.END} {output_file}")
        print(f"{Colors.BLUE}Execution time:{Colors.END} {duration:.2f} seconds ({duration/60:.1f} minutes)")
        print(f"{Colors.BLUE}Finished:{Colors.END} {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{Colors.BOLD}{Colors.GREEN}{'='*70}{Colors.END}\n")
        
        return True
        
    except Exception as e:
        print(f"\n{Colors.BOLD}{Colors.RED}{'='*70}{Colors.END}")
        print(f"{Colors.BOLD}{Colors.RED} PIPELINE FAILED{Colors.END}")
        print(f"{Colors.BOLD}{Colors.RED}{'='*70}{Colors.END}")
        print(f"{Colors.RED}Error: {str(e)}{Colors.END}\n")
        return False

# ==============================================================================
# CLI INTERFACE
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Phylo-Pipe',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Neighbor-Joining (fast)
  python run_phylogeny.py -i sequences.fasta -m nj
  
  # Maximum Likelihood with 2000 bootstraps
  python run_phylogeny.py -i sequences.fasta -m ml --bootstrap 2000
  
  # Bayesian Inference with 500k generations
  python run_phylogeny.py -i sequences.fasta -m bayesian --ngen 500000
  
  # Custom MAFFT strategy
  python run_phylogeny.py -i sequences.fasta -m ml --mafft-strategy localpair
"""
    )
    
    # Required arguments
    parser.add_argument(
        '-i', '--input',
        type=Path,
        required=True,
        help='Input FASTA file with sequences'
    )
    
    parser.add_argument(
        '-m', '--method',
        choices=['nj', 'ml', 'bayesian'],
        required=True,
        help='Phylogenetic method (nj=Neighbor-Joining, ml=Maximum Likelihood, bayesian=Bayesian Inference)'
    )
    
    # MAFFT parameters
    parser.add_argument(
        '--mafft-strategy',
        default='auto',
        choices=['auto', 'localpair', 'globalpair', 'retree', 'fftns'],
        help='MAFFT alignment strategy (default: auto)'
    )
    
    # trimAl parameters
    parser.add_argument(
        '--trimal-method',
        default='automated1',
        choices=['automated1', 'strict', 'gappyout'],
        help='trimAl filtering method (default: automated1)'
    )
    
    # IQ-TREE parameters
    parser.add_argument(
        '--iqtree-model',
        default='MFP',
        help='IQ-TREE substitution model (default: MFP for ModelFinder)'
    )
    
    parser.add_argument(
        '--bootstrap',
        type=int,
        default=1000,
        help='Number of bootstrap replicates for ML (default: 1000)'
    )
    
    # MrBayes parameters
    parser.add_argument(
        '--ngen',
        type=int,
        default=100000,
        help='Number of MCMC generations for Bayesian (default: 100000)'
    )
    
    parser.add_argument(
        '--nchains',
        type=int,
        default=4,
        help='Number of MCMC chains for Bayesian (default: 4)'
    )
    
    # General options
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Suppress detailed output'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='Phylo-Pipe'
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not args.input.exists():
        print(f"{Colors.RED}Error: Input file '{args.input}' not found{Colors.END}")
        sys.exit(1)
    
    # Build parameters dictionary
    params = {
        'mafft_strategy': args.mafft_strategy,
        'trimal_method': args.trimal_method,
        'iqtree_model': args.iqtree_model,
        'iqtree_bootstrap': args.bootstrap,
        'mrbayes_ngen': args.ngen,
        'mrbayes_nchains': args.nchains
    }
    
    # Run pipeline
    success = run_pipeline(
        args.input,
        args.method,
        params,
        verbose=not args.quiet
    )
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()