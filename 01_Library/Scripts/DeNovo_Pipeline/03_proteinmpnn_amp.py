import os, sys, argparse, yaml, sqlite3, subprocess, glob, random
from datetime import datetime
from Bio import SeqIO
import pandas as pd

def trigger_fatal_error(step, reason, action, log_path=None):
    err_snippet = ""
    log_info = f"\n 로그 경로 : {log_path}" if log_path else ""
    if log_path and os.path.exists(log_path):
        try:
            with open(log_path, 'r') as f:
                lines = f.readlines()
                if lines:
                    err_snippet = "\n\n [원본 에러 메시지 (최근 10줄)]:\n" + "".join(lines[-10:])
        except Exception:
            pass
    print(f"\n==================================================\n [FATAL ERROR] {step} 가동 중단!\n 원인 : {reason}\n 조치 : {action}{log_info}{err_snippet}\n==================================================\n")
    sys.exit(1)

def load_config(config_path):
    with open(config_path, 'r') as f: return yaml.safe_load(f)

def run_proteinmpnn(config_path, mode):
    start_time = datetime.now()
    run_id = start_time.strftime("RUN_%Y%m%d_%H%M%S")

    config = load_config(config_path)
    target_id = config.get('project', {}).get('target_pdb_id', '6M0J')
    project_name = config.get('project', {}).get('project_name', target_id)
    global_seed = config.get('project', {}).get('global_seed', 42)
    mpnn_seed = config.get('advanced_overrides', {}).get('mpnn_seed', global_seed)
    random.seed(mpnn_seed)
    
    workspace_dir = config.get('paths', {}).get('base_dir', os.path.join("03_Workspace", project_name))
    db_path = config.get('paths', {}).get('db_path', os.path.join(workspace_dir, "00_Master_Control", "insilico_master.db"))
    current_run_dir = config.get('paths', {}).get('current_run_dir', workspace_dir)
    
    mpnn_out_dir = os.path.join(current_run_dir, "03_ProteinMPNN")
    log_dir = os.path.join(mpnn_out_dir, "00_logs"); os.makedirs(log_dir, exist_ok=True)
    symlink_dir = os.path.join(mpnn_out_dir, "01_scaffold_links"); os.makedirs(symlink_dir, exist_ok=True)
    fasta_dir = os.path.join(mpnn_out_dir, "02_fasta_outputs"); os.makedirs(fasta_dir, exist_ok=True)
    report_dir = os.path.join(mpnn_out_dir, "03_seq_reports"); os.makedirs(report_dir, exist_ok=True)

    engine_log = os.path.join(log_dir, f"proteinmpnn_engine_{run_id}_SEED{mpnn_seed}.log")
    with open(engine_log, 'w') as f:
        f.write(f"=== SCRIPT START: {start_time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        f.write(f"=== SCRIPT END: [RUNNING...] ===\n")
        f.write(f"=== RUN ID: {run_id} | GLOBAL SEED: {global_seed} ===\n")
        f.write(f"=================================\n\n")

    conn = sqlite3.connect(db_path); cur = conn.cursor()

    if mode == 'resume':

        cur.execute("SELECT COUNT(*) FROM Binders WHERE target_id = ? AND generation = 'DeNovo_MPNN'", (target_id,))
        if cur.fetchone()[0] > 0: 
            print("    [Mode: RESUME] Found existing sequences. Bypassing...")
            with open(engine_log, 'a') as lf: lf.write("Bypassed via RESUME mode.\n")
            conn.close(); finalize_log(engine_log, start_time); return

    if mode == 'overwrite':
        print("    [Mode: OVERWRITE] Soft Overwrite Active - 기존 서열 보존.")

    cur.execute("SELECT binder_id FROM Binders WHERE status IN ('Passed_RF_RMSD', 'Passed_MPNN', 'Passed_FastQC', 'Passed_FoldX', 'Passed_Radar', 'Passed_Radar_Champion') AND target_id = ? AND generation = 'DeNovo_RF'", (target_id,))
    scaffolds = cur.fetchall()
    if not scaffolds: trigger_fatal_error("[Phase 2-3]", "No scaffolds found.", "2호기 확인")

    print(f"    [Amplification]  IGNITING AI: Generating sequences (Run ID: {run_id})...")
    mpnn_script = os.environ.get("PROTEINMPNN_PATH", "/opt/ProteinMPNN/protein_mpnn_run.py")

    ai_python = os.environ.get("AI_PYTHON", "python3")

    success_count = 0
    total_scaffolds = len(scaffolds)
    
    for idx, row in enumerate(scaffolds):
        print(f"\r      -> [ProteinMPNN] Generating sequences: {idx+1} / {total_scaffolds}...", end="", flush=True)

        scaffold_id = row[0]
        

        if mode == 'resume':
            cur.execute("SELECT COUNT(*) FROM Binders WHERE parent_id = ? AND generation = 'DeNovo_MPNN'", (scaffold_id,))
            if cur.fetchone()[0] > 0:
                continue

        src_pdb = os.path.join(current_run_dir, "02_RF_Scaffolds", "01_generated_pdbs", f"{scaffold_id}.pdb")
        link_pdb = os.path.join(symlink_dir, f"{scaffold_id}.pdb")
        
        if os.path.exists(src_pdb):
            if not os.path.exists(link_pdb): os.symlink(os.path.abspath(src_pdb), link_pdb)
        else: continue
            
        fasta_file = os.path.join(fasta_dir, f"{scaffold_id}.fa")
        
        if not os.path.exists(fasta_file):
            try:
                with open(engine_log, 'a') as lf:
                    subprocess.run([ai_python, mpnn_script, "--pdb_path", link_pdb, "--out_folder", fasta_dir, "--num_seq_per_target", "10", "--sampling_temp", "0.1", "--batch_size", "1", "--seed", str(mpnn_seed)], check=True, stdout=lf, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                pass

        actual_fasta = os.path.join(fasta_dir, "seqs", f"{scaffold_id}.fa")
        if os.path.exists(actual_fasta):
            records = list(SeqIO.parse(actual_fasta, "fasta"))
            for i, rec in enumerate(records):
                if i == 0: continue
                new_binder_id, designed_seq = f"{scaffold_id}_seq{i}", str(rec.seq)
                cur.execute("INSERT OR IGNORE INTO Binders (binder_id, target_id, sequence, status, parent_id, generation) VALUES (?, ?, ?, 'Passed_MPNN', ?, 'DeNovo_MPNN')", (new_binder_id, target_id, designed_seq, scaffold_id))
                for table in ['QC_Physiochemical', 'QC_Interaction', 'QC_Structural', 'QC_Production']: cur.execute(f"INSERT OR IGNORE INTO {table} (binder_id) VALUES (?)", (new_binder_id,))
                success_count += 1
                
    if total_scaffolds > 0: print()
    
    conn.commit()
    cur.execute("SELECT binder_id, parent_id, sequence FROM Binders WHERE status != 'Failed_RF_RMSD' AND target_id = ? AND generation = 'DeNovo_MPNN'", (target_id,))
    pd.DataFrame(cur.fetchall(), columns=['Binder_ID', 'Parent_Scaffold', 'Sequence']).to_csv(os.path.join(report_dir, f"mpnn_output_sequences_{run_id}.csv"), index=False)
    
    conn.close()
    if success_count == 0: trigger_fatal_error("[Phase 3]", "No sequences generated.", "로그 확인")
    print(f"    [SUCCESS] {success_count} AI sequences extracted.")
    
    finalize_log(engine_log, start_time)

def finalize_log(log_path, start_time):
    end_time = datetime.now()
    duration = end_time - start_time
    with open(log_path, 'r') as f: lines = f.readlines()
    if len(lines) > 1: lines[1] = f"=== SCRIPT END: {end_time.strftime('%Y-%m-%d %H:%M:%S')} (Duration: {duration}) ===\n"
    with open(log_path, 'w') as f: f.writelines(lines)

def main(config_path):
    config = load_config(config_path)
    mode = config.get('pipeline_control', {}).get('step03_proteinmpnn', 'overwrite').lower()
    if mode == 'skip': sys.exit(0)
    print(f">>> [Phase 3] ProteinMPNN Sequence Amplification Started...")
    run_proteinmpnn(config_path, mode)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); main(parser.parse_args().config)
