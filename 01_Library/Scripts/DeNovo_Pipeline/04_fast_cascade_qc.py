import os, sys, argparse, yaml, sqlite3, subprocess, tempfile, random
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def trigger_fatal_error(step, reason, action, log_path=None):
    err_snippet = ""
    log_info = f"\n📄 로그 경로 : {log_path}" if log_path else ""
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

def run_fast_cascade_qc(config_path, mode):
    start_time = datetime.now()
    run_id = start_time.strftime("RUN_%Y%m%d_%H%M%S")

    config = load_config(config_path)
    target_id = config.get('project', {}).get('target_pdb_id', '6M0J')
    if not target_id: target_id = config.get('target', {}).get('pdb_id', '6M0J')
    project_name = config.get('project', {}).get('project_name', target_id)
    global_seed = config.get('project', {}).get('global_seed', 42)
    random.seed(global_seed)

    workspace_dir = config.get('paths', {}).get('base_dir', os.path.join("03_Workspace", project_name))
    db_path = config.get('paths', {}).get('db_path', os.path.join(workspace_dir, "00_Master_Control", "insilico_master.db"))
    current_run_dir = config.get('paths', {}).get('current_run_dir', workspace_dir)
    
    qc_dir = os.path.join(current_run_dir, "04_Fast_QC")
    log_dir = os.path.join(qc_dir, "00_logs"); os.makedirs(log_dir, exist_ok=True)
    link_dir = os.path.join(qc_dir, "01_target_links"); os.makedirs(link_dir, exist_ok=True)
    report_dir = os.path.join(qc_dir, "02_qc_reports"); os.makedirs(report_dir, exist_ok=True)
    
    c = config.get('filtering_strictness', {})
    max_instability, max_a3d, phys_top_percent = c.get('max_instability', 40.0), c.get('max_a3d_score', 1.0), c.get('physicochemical_top_percent', 100.0)

    overrides = config.get('advanced_overrides', {})
    if 'netmhcpan_alleles' in overrides and overrides['netmhcpan_alleles']:
        hla_alleles = overrides['netmhcpan_alleles']
    else:
        ref_path = os.path.join("01_Library", "References", "HLA_Alleles", "global_standard.yaml")
        if os.path.exists(ref_path):
            with open(ref_path, 'r') as f: hla_alleles = yaml.safe_load(f).get('alleles', ["HLA-A*02:01"])
        else:
            hla_alleles = ["HLA-A*02:01"]
        hla_param_str = ",".join(hla_alleles).replace("*", "")

    engine_log = os.path.join(log_dir, f"04_cascade_qc_engine_{run_id}.log") 
    with open(engine_log, 'w') as f:
        f.write(f"=== SCRIPT START: {start_time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        f.write(f"=== SCRIPT END: [RUNNING...] ===\n")
        f.write(f"=== RUN ID: {run_id} | HLA_TARGETS: {hla_param_str} ===\n")
        f.write(f"===============================\n\n")

    conn = sqlite3.connect(db_path); cur = conn.cursor()

    if mode == 'resume':
        cur.execute("SELECT COUNT(*) FROM Binders WHERE status IN ('Passed_FastQC', 'Passed_FoldX', 'Failed_FoldX_Abs', 'Failed_FoldX_Rel', 'Passed_Radar', 'Passed_Radar_Champion') AND target_id = ?", (target_id,))
        if cur.fetchone()[0] > 0: 
            print(f"    [Mode: RESUME] FastQC already completed. Bypassing...")
            conn.close(); finalize_log(engine_log, start_time); return


    if mode == 'overwrite':

        cur.execute("UPDATE Binders SET status = 'Passed_MPNN' WHERE status LIKE 'Failed%' AND generation = 'DeNovo_MPNN' AND target_id = ?", (target_id,))
        conn.commit()

    cur.execute("SELECT b.binder_id, b.sequence, p.instability, b.parent_id FROM Binders b JOIN QC_Physiochemical p ON b.binder_id = p.binder_id WHERE b.status = 'Passed_MPNN' AND b.target_id = ?", (target_id,))
    candidates = cur.fetchall()
    if not candidates: trigger_fatal_error("[Phase 3]", "No 'Passed_MPNN' sequences found.", "3호기 확인")

    print(f"    [FastQC] 🔬 IGNITING QC CASCADE for {len(candidates)} sequences (Run ID: {run_id})...")
    a3d_script = os.environ.get("A3D_PATH", "aggrescan3d")

    netmhcpan_bin = os.environ.get("NETMHCPAN_PATH", "/usr/local/bin/netMHCpan")


    cur.execute("SELECT COUNT(*) FROM Binders WHERE generation = 'DeNovo_MPNN' AND target_id = ?", (target_id,))
    total_mpnn_generated = cur.fetchone()[0]
    already_done_qc = total_mpnn_generated - len(candidates)


    max_immuno = 2.0  
    total_len = len(candidates)
    
    for idx, (b_id, raw_seq, existing_instability, p_id) in enumerate(candidates): 
        print(f"\r      -> Calculating Metrics: {idx+1} / {total_len}...", end="", flush=True)
        

        if '/' in raw_seq:
            target_seq, binder_seq = raw_seq.rsplit('/', 1)
        else:
            binder_seq = raw_seq
            
        clean_seq = binder_seq.replace("X", "").replace("-", "").strip()
        

        sim_instab = round(ProteinAnalysis(clean_seq).instability_index(), 2)
        cur.execute("UPDATE QC_Physiochemical SET instability=? WHERE binder_id=?", (sim_instab, b_id))
        

        src_pdb = os.path.join(current_run_dir, "02_RF_Scaffolds", "01_generated_pdbs", f"{p_id}.pdb")
        link_pdb = os.path.join(link_dir, f"{p_id}.pdb")
        if os.path.exists(src_pdb) and not os.path.exists(link_pdb): os.symlink(os.path.abspath(src_pdb), link_pdb)
        try:
            res_a3d = subprocess.run(["python3", a3d_script, "-i", link_pdb], capture_output=True, text=True)
            sim_a3d = float(res_a3d.stdout.split("A3D Score:")[1].split()[0]) if "A3D Score:" in res_a3d.stdout else round(ProteinAnalysis(clean_seq).gravy(), 3)
        except Exception: sim_a3d = 99.9
        cur.execute("UPDATE QC_Interaction SET a3d_score=? WHERE binder_id=?", (sim_a3d, b_id))
        

        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
                tmp.write(f">{b_id}\n{clean_seq}\n"); tmp_fasta = tmp.name

            res_im = subprocess.run([netmhcpan_bin, "-f", tmp_fasta, "-a", hla_param_str], capture_output=True, text=True)
            os.remove(tmp_fasta)
            if "Rank" not in res_im.stdout: 
                sim_immuno = 99.9
                with open(engine_log, 'a') as lf:
                    lf.write(f"\n[NetMHCpan ERROR - {b_id}]\nSTDOUT: {res_im.stdout}\nSTDERR: {res_im.stderr}\n")
            else:
                epitope_count = res_im.stdout.count("<= SB") + res_im.stdout.count("<= WB")
                sim_immuno = float(epitope_count)
        except Exception as e: 
            sim_immuno = 99.9
            with open(engine_log, 'a') as lf: lf.write(f"\n[Python ERROR - {b_id}] {str(e)}\n")
        cur.execute("UPDATE QC_Interaction SET immunogenicity_rank=? WHERE binder_id=?", (sim_immuno, b_id))
        
    print("\n    [FastQC] 📊 Applying Cutoffs and Generating Reports...")
    
   
    report_query = f"SELECT b.binder_id, p.instability, i.a3d_score, i.immunogenicity_rank FROM Binders b JOIN QC_Physiochemical p ON b.binder_id = p.binder_id JOIN QC_Interaction i ON b.binder_id = i.binder_id WHERE b.target_id = '{target_id}' AND b.generation = 'DeNovo_MPNN'"
    df = pd.read_sql_query(report_query, conn)
    
    # 99.9(에러)는 무조건 탈락 처리
    df['Status'] = 'Failed'
    pass_condition = (df['instability'] <= max_instability) & (df['a3d_score'] <= max_a3d) & (df['immunogenicity_rank'] <= max_immuno) & (df['immunogenicity_rank'] != 99.9) 
    df.loc[pass_condition, 'Status'] = 'Passed_FastQC'

    for _, row in df.iterrows():
        cur.execute("UPDATE Binders SET status = ? WHERE binder_id = ?", (row['Status'], row['binder_id']))
    

    df.to_csv(os.path.join(report_dir, f"fast_cascade_ALL_report_{run_id}.csv"), index=False)
    passed_df = df[df['Status'] == 'Passed_FastQC'].copy()
    passed_df = passed_df.sort_values(by=['immunogenicity_rank', 'a3d_score', 'instability']) 
    passed_df.to_csv(os.path.join(report_dir, f"PASSED_fastqc_list_{run_id}.csv"), index=False)
    

    plt.figure(figsize=(18, 5))
    colors = {'Passed_FastQC': 'green', 'Failed': 'red'}
    
    plt.subplot(1, 3, 1)
    for status in colors:
        subset = df[df['Status'] == status]
        plt.scatter(subset.index, subset['instability'], color=colors[status], label=status, alpha=0.6)
    plt.axhline(y=max_instability, color='blue', linestyle='--', label=f'Cutoff ({max_instability})')
    plt.title('Instability'); plt.legend()

    plt.subplot(1, 3, 2)
    for status in colors:
        subset = df[(df['Status'] == status) & (df['a3d_score'] != 99.9)]
        plt.scatter(subset.index, subset['a3d_score'], color=colors[status], label=status, alpha=0.6)
    plt.axhline(y=max_a3d, color='blue', linestyle='--', label=f'Cutoff ({max_a3d})')
    plt.title('Aggregation (A3D)'); plt.legend()

    plt.subplot(1, 3, 3)
    for status in colors:
        subset = df[(df['Status'] == status) & (df['immunogenicity_rank'] != 99.9)]
        plt.scatter(subset.index, subset['immunogenicity_rank'], color=colors[status], label=status, alpha=0.6)
    plt.axhline(y=max_immuno, color='blue', linestyle='--', label='Cutoff (0.0)')
    plt.title('Immunogenicity (Epitopes)'); plt.legend()

    plt.tight_layout()
    graph_path = os.path.join(report_dir, f"fastqc_metrics_graph_{run_id}.png")
    plt.savefig(graph_path, dpi=300)
    plt.close()

    conn.commit(); conn.close()
    print(f"    [SUCCESS] Analysis Complete. {len(passed_df)} survivors advance.")
    print(f"     Saved ALL Report: fast_cascade_ALL_report_{run_id}.csv")
    print(f"     Saved PASSED List: PASSED_fastqc_list_{run_id}.csv")
    print(f"     Saved Metrics Graph: {graph_path}")
    finalize_log(engine_log, start_time)

def finalize_log(log_path, start_time):
    end_time = datetime.now()
    duration = end_time - start_time
    with open(log_path, 'r') as f: lines = f.readlines()
    if len(lines) > 1: lines[1] = f"=== SCRIPT END: {end_time.strftime('%Y-%m-%d %H:%M:%S')} (Duration: {duration}) ===\n"
    with open(log_path, 'w') as f: f.writelines(lines)

def main(config_path):
    config = load_config(config_path)
    mode = config.get('pipeline_control', {}).get('step04_fast_cascade', 'overwrite').lower()
    if mode == 'skip': sys.exit(0)
    print(f">>> [Phase 4] Real Fast QC Cascade Started...")
    run_fast_cascade_qc(config_path, mode)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); main(parser.parse_args().config)
