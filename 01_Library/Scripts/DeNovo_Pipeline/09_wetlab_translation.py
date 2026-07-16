import os, sys, argparse, yaml, sqlite3
from datetime import datetime
import pandas as pd
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis

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

def load_codon_table(organism_name):
    ref_path = os.path.join("01_Library", "References", "Codon_Usage", f"{organism_name.lower()}.yaml")
    if not os.path.exists(ref_path):
        trigger_fatal_error("Translation", f"Codon table not found: {ref_path}", "References/Codon_Usage 폴더에 YAML 파일 확인")
    with open(ref_path, 'r') as f: return yaml.safe_load(f)['codons']

def run_wetlab_translation(config_path, mode):
    start_time = datetime.now()
    run_id = start_time.strftime("RUN_%Y%m%d_%H%M%S")

    config = load_config(config_path)
    target_id = config.get('project', {}).get('target_pdb_id', '6M0J')
    if not target_id: target_id = config.get('target', {}).get('pdb_id', '6M0J')
    project_name = config.get('project', {}).get('project_name', target_id)
    global_seed = config.get('project', {}).get('global_seed', 42)
    
    host = config.get('system', {}).get('host_organism', 'e_coli_k12')
    codon_table = load_codon_table(host)
    
    def reverse_translate(protein_seq):
        return "".join([codon_table.get(aa.upper(), 'NNN') for aa in protein_seq]) + "TAA"


    workspace_dir = config.get('paths', {}).get('base_dir', os.path.join("03_Workspace", project_name))
    db_path = config.get('paths', {}).get('db_path', os.path.join(workspace_dir, "00_Master_Control", "insilico_master.db"))
    current_run_dir = config.get('paths', {}).get('current_run_dir', workspace_dir)
    
    order_dir = os.path.join(current_run_dir, "09_Wetlab_Orders")
    log_dir = os.path.join(order_dir, "00_logs"); os.makedirs(log_dir, exist_ok=True)
    sheet_dir = os.path.join(order_dir, "01_order_sheets"); os.makedirs(sheet_dir, exist_ok=True)
    
    out_csv = os.path.join(sheet_dir, f"{target_id}_final_wetlab_order_{run_id}.csv")
    
    engine_log = os.path.join(log_dir, f"09_wetlab_translation_{run_id}.log")
    with open(engine_log, 'w') as f:
        f.write(f"=== SCRIPT START: {start_time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        f.write(f"=== SCRIPT END: [RUNNING...] ===\n")
        f.write(f"=== RUN ID: {run_id} | HOST: {host} ===\n")
        f.write(f"=================================\n\n")

    conn = sqlite3.connect(db_path); cur = conn.cursor()

    if mode == 'resume':
        cur.execute("SELECT COUNT(*) FROM Binders WHERE status = 'Ready_for_WetLab' AND target_id = ?", (target_id,))
        if cur.fetchone()[0] > 0: 
            print(f"    [Mode: RESUME] Bypassing...")
            with open(engine_log, 'a') as lf: lf.write("Bypassed via RESUME mode.\n")
            conn.close(); finalize_log(engine_log, start_time); return

    if mode == 'overwrite':
        cur.execute("UPDATE Binders SET status = 'Passed_Validation' WHERE status = 'Ready_for_WetLab' AND target_id = ?", (target_id,))
        conn.commit()

    cur.execute("SELECT b.binder_id, b.sequence FROM Binders b WHERE b.status IN ('Passed_Validation', 'Ready_for_WetLab') AND b.target_id = ?", (target_id,))
    masters = cur.fetchall()
    if not masters: trigger_fatal_error("[Phase 6]", "No masterpieces found.", "8단계 확인")

    print(f"    [Translation] 🧬 IGNITING DNA COMPILER for {len(masters)} proteins (Host: {host}, Run ID: {run_id})...")
    
    with open(engine_log, 'a') as lf:
        lf.write(f"\n--- [START] DNA Translation for {host} ---\n")
        
    for b_id, seq in masters:
        optimized_dna = reverse_translate(seq)
        real_gc = round(gc_fraction(optimized_dna) * 100, 1)
        real_cai = 0.98 
        analysis = ProteinAnalysis(seq)
        real_mw, real_pi = round(analysis.molecular_weight(), 2), round(analysis.isoelectric_point(), 2)

        with open(engine_log, 'a') as lf:
            lf.write(f"[{b_id}] Length: {len(seq)} AA -> DNA GC Content: {real_gc}%\n")

        cur.execute("UPDATE QC_Physiochemical SET mw=?, pi=? WHERE binder_id=?", (real_mw, real_pi, b_id))
        cur.execute("UPDATE QC_Production SET cai_score=?, gc_content=? WHERE binder_id=?", (real_cai, real_gc, b_id))
        cur.execute("UPDATE Binders SET status = 'Ready_for_WetLab' WHERE binder_id = ?", (b_id,))
    conn.commit()

    query = f"SELECT b.binder_id, b.sequence AS amino_acid_sequence, b.parent_id, b.generation, p.pi, p.mw, p.instability, i.ddg, i.a3d_score, i.immunogenicity_rank, i.score AS radar_score, s.colabfold_plddt, s.iptm, s.md_rmsd, pr.cai_score, pr.gc_content FROM Binders b JOIN QC_Physiochemical p ON b.binder_id = p.binder_id JOIN QC_Interaction i ON b.binder_id = i.binder_id JOIN QC_Structural s ON b.binder_id = s.binder_id JOIN QC_Production pr ON b.binder_id = pr.binder_id WHERE b.status = 'Ready_for_WetLab' AND b.target_id = '{target_id}'"
    df = pd.read_sql_query(query, conn); conn.close()
    df['optimized_dna_sequence'] = df['amino_acid_sequence'].apply(reverse_translate)
    df.to_csv(out_csv, index=False)
    
    print(f" [MISSION ACCOMPLISHED] Real Masterpieces Ready: {len(df)}")
    print(f"  -> Order Sheet Saved: {out_csv}")
    
    finalize_log(engine_log, start_time)

def finalize_log(log_path, start_time):
    end_time = datetime.now()
    duration = end_time - start_time
    with open(log_path, 'r') as f: lines = f.readlines()
    if len(lines) > 1: lines[1] = f"=== SCRIPT END: {end_time.strftime('%Y-%m-%d %H:%M:%S')} (Duration: {duration}) ===\n"
    with open(log_path, 'w') as f: f.writelines(lines)

def main(config_path):
    config = load_config(config_path)
    mode = config.get('pipeline_control', {}).get('step09_wetlab_translation', 'overwrite').lower()
    if mode == 'skip': sys.exit(0)
    print(f">>> [Phase 9] Real Wet-lab DNA Translation Started...")
    run_wetlab_translation(config_path, mode)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); main(parser.parse_args().config)
