#!/usr/bin/env python3
import os, sys, argparse, yaml, sqlite3, glob, subprocess, random
from datetime import datetime
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

def run_rfdiffusion_triage(config_path, mode):
    start_time = datetime.now()
    run_id = start_time.strftime("RUN_%Y%m%d_%H%M%S")

    config = load_config(config_path)
    sub1 = config.get('micro_control', {}).get('step02', {}).get('sub_1', mode)
    sub2 = config.get('micro_control', {}).get('step02', {}).get('sub_2', mode)
    target_id = config.get('project', {}).get('target_pdb_id', '6M0J')
    project_name = config.get('project', {}).get('project_name', target_id)
    global_seed = config.get('project', {}).get('global_seed', 42)
    rfd_seed = config.get('rfdiffusion', {}).get('seed', global_seed)
    random.seed(rfd_seed)
    
    num_designs = config.get('rfdiffusion', {}).get('num_designs', 10000)
    contig_str = config.get('rfdiffusion', {}).get('contig', '15-65')    
    abs_cutoff = config.get('filtering_strictness', {}).get('rmsd_absolute_cutoff', 2.5)
    top_percent = config.get('filtering_strictness', {}).get('rmsd_relative_top_percent', 1.0)
    
    workspace_dir = config.get('paths', {}).get('base_dir', os.path.join("03_Workspace", project_name))
    db_path = config.get('paths', {}).get('db_path', os.path.join(workspace_dir, "00_Master_Control", "insilico_master.db"))
    current_run_dir = config.get('paths', {}).get('current_run_dir', workspace_dir)
    
    scaffold_dir = os.path.join(current_run_dir, "02_RF_Scaffolds")
    log_dir = os.path.join(scaffold_dir, "00_logs"); os.makedirs(log_dir, exist_ok=True)
    pdb_dir = os.path.join(scaffold_dir, "01_generated_pdbs"); os.makedirs(pdb_dir, exist_ok=True)
    plot_dir = os.path.join(scaffold_dir, "02_rmsd_plots"); os.makedirs(plot_dir, exist_ok=True)
    report_dir = os.path.join(scaffold_dir, "03_survivors"); os.makedirs(report_dir, exist_ok=True)
    
    engine_log = os.path.join(log_dir, f"rfdiffusion_engine_{run_id}_SEED{rfd_seed}.log")
    with open(engine_log, 'w') as f:
        f.write(f"=== SCRIPT START: {start_time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        f.write(f"=== SCRIPT END: [RUNNING...] ===\n")
        f.write(f"=== RUN ID: {run_id} | GLOBAL SEED: {global_seed} ===\n")
        f.write(f"=== CONFIG INFO (RFdiffusion) ===\n{yaml.dump(config.get('rfdiffusion', {}))}\n=================================\n\n")

    target_pdb = os.path.join(current_run_dir, "01_Target_Profile", "01_raw_pdb", f"{target_id}_raw.pdb")
    
    if not os.path.exists(target_pdb):
        trigger_fatal_error("[Phase 2] PDB File Missing", f"Target PDB가 존재하지 않습니다: {target_pdb}", "1호기(Target Profiling)를 'OVERWRITE' 모드로 실행하여 PDB를 먼저 다운로드하십시오.", engine_log)
    
    conn = sqlite3.connect(db_path); cur = conn.cursor()


    cur.execute("SELECT COUNT(*) FROM Binders WHERE target_id = ? AND generation = 'DeNovo_RF'", (target_id,))
    already_done = cur.fetchone()[0]
    out_prefix = os.path.join(pdb_dir, f"design_AF_{target_id}")

    if sub1 == 'overwrite':
        print("    [Phase 2-1: OVERWRITE] 뼈대 생성 초기화 (기존 기록 무시하고 덮어씁니다).")
        already_done = 0
    elif sub1 == 'resume':
        print(f"    [Phase 2-1: RESUME] 기존 생성된 뼈대 파일을 유지합니다. (현재 {already_done}개 완료됨)")

    if sub1 != 'skip':
        if already_done < num_designs:
            designs_to_add = num_designs - already_done
            print(f"    [Generation]  IGNITING GPU: Designing {designs_to_add} additional scaffolds (Total Target: {num_designs})...")
            rfdiffusion_script = os.environ.get("RFDIFFUSION_PATH", "/opt/RFdiffusion/run_inference.py")
            ai_python = os.environ.get("AI_PYTHON", "python3")
            
            try:
                with open(engine_log, 'a') as lf:
                    process = subprocess.Popen([
                        ai_python, rfdiffusion_script, 
                        f"inference.output_prefix={out_prefix}", 
                        f"inference.num_designs={num_designs}", # 총 목표 개수 전달
                        f"inference.input_pdb={target_pdb}", 
                        f"contigmap.contigs=[{contig_str}]", 
                        f"inference.deterministic=False",
                        f"+inference.seed={rfd_seed}" 
                    ], stdout=lf, stderr=subprocess.STDOUT)
                    
                    import time
                    while process.poll() is None:
                        count = len(glob.glob(f"{out_prefix}_*.pdb"))
                        print(f"\r      -> [Generation Progress] {count} / {num_designs} scaffolds complete...", end="", flush=True)
                        time.sleep(2)
                    
                    final_count = len(glob.glob(f"{out_prefix}_*.pdb"))
                    print(f"\r      -> [Generation Progress] {final_count} / {num_designs} scaffolds complete... Done!\n", flush=True)
                    
            except Exception as e:
                trigger_fatal_error("[Phase 2] GPU Engine", f"RFdiffusion crashed: {e}", f"로그 파일 점검", engine_log)
        else:
            print(f"    [Generation Progress] 목표 스캐폴드({num_designs}개)가 이미 모두 생성되어 있습니다.")
    else:
        print("    [Phase 2-1: SKIP] 뼈대 생성 단계를 건너뜁니다.")

    generated_pdbs = glob.glob(f"{out_prefix}_*.pdb")
    for pdb_path in generated_pdbs:
        b_id = os.path.basename(pdb_path).replace('.pdb', '')
        cur.execute("INSERT OR IGNORE INTO Binders (binder_id, target_id, sequence, status, parent_id, generation) VALUES (?, ?, 'AWAITING_SEQ', 'Failed_RF_RMSD', 'NONE', 'DeNovo_RF')", (b_id, target_id))
        cur.execute("INSERT OR IGNORE INTO QC_Structural (binder_id) VALUES (?)", (b_id,))
    conn.commit()

    if sub2 == 'skip':
        print("    [Phase 2-2: SKIP] 필터링 단계를 건너뜁니다.")
        conn.close(); finalize_log(engine_log, start_time); return

    if sub2 == 'overwrite':
        print("    [Phase 2-2: OVERWRITE] 필터링 기준 재평가를 위해 기존 DB 상태를 초기화합니다.")
        cur.execute("UPDATE Binders SET status = 'Generated_RF' WHERE target_id = ? AND generation = 'DeNovo_RF'", (target_id,))
        cur.execute("UPDATE QC_Structural SET rfd_rmsd = 99.9 WHERE binder_id IN (SELECT binder_id FROM Binders WHERE target_id = ? AND generation = 'DeNovo_RF')", (target_id,))
        conn.commit()
    elif sub2 == 'resume':
        cur.execute("SELECT COUNT(*) FROM Binders WHERE status != 'Failed_RF_RMSD' AND target_id = ? AND generation = 'DeNovo_RF'", (target_id,))
        if cur.fetchone()[0] > 0: 
            print(f"    [Phase 2-2: RESUME] 필터링이 이미 완료되었습니다. (Bypassing...)")
            conn.close(); finalize_log(engine_log, start_time); return

    print("    [Phase 2-2] Triage Started (Extracting metrics & pLDDT...)")

    cur.execute("SELECT b.binder_id, s.rfd_rmsd FROM Binders b JOIN QC_Structural s ON b.binder_id = s.binder_id WHERE b.target_id = ? AND b.generation = 'DeNovo_RF'", (target_id,))
    scaffolds = cur.fetchall()
    processed_scaffolds = []
    
    for b_id, rmsd in scaffolds:
        if rmsd is None or rmsd == 99.9:
            trb_file = os.path.join(pdb_dir, f"{b_id}.trb")
            real_rmsd = 99.9
            if os.path.exists(trb_file):

                if os.path.getsize(trb_file) == 0:
                    with open(engine_log, 'a') as lf: 
                        lf.write(f"TRB Error for {b_id}: File is 0 bytes (Disk full or interrupted).\n")
                    real_rmsd = 99.9
                else:
                    import torch, pickle, numpy as np
                    try:

                        data = torch.load(trb_file, map_location='cpu', weights_only=False)
                    except Exception as e1:
                        try:

                            with open(trb_file, 'rb') as tf: data = pickle.load(tf)
                        except Exception as e2:
                            try:

                                data = np.load(trb_file, allow_pickle=True).item()
                            except Exception as e3:
                                with open(engine_log, 'a') as lf: 
                                    lf.write(f"TRB Load Error for {b_id}: Torch/Pickle/Numpy all failed.\n")
                                data = {} # 에러 방지용 빈 딕셔너리
                                

                    val = 99.9
                    for k in ['inpaint_seq_rmsd', 'rmsd', 'complex_rmsd', 'binder_rmsd']:
                        if k in data:
                            val = data[k]
                            break
                    
                    if val == 99.9:
                        rmsd_keys = [k for k in data.keys() if isinstance(k, str) and 'rmsd' in k.lower()]
                        if rmsd_keys: val = data[rmsd_keys[0]]
                    

                    if val == 99.9 and 'plddt' in data:
                        import numpy as np
                        plddt_array = np.array(data['plddt'])
                        mean_plddt = np.mean(plddt_array) # 0 ~ 100 점 (높을수록 우수)
                        

                        if mean_plddt <= 1.0:
                            mean_plddt = mean_plddt * 100.0


                        val = (100.0 - mean_plddt) / 10.0
                    

                    if hasattr(val, 'item'): val = val.item()
                    real_rmsd = round(float(val), 3)
            
            cur.execute("UPDATE QC_Structural SET rfd_rmsd = ? WHERE binder_id = ?", (real_rmsd, b_id))
            processed_scaffolds.append((b_id, real_rmsd))
        else: processed_scaffolds.append((b_id, rmsd))
    conn.commit()

    absolute_passed = [s for s in processed_scaffolds if s[1] <= abs_cutoff]
    absolute_failed = [s for s in processed_scaffolds if s[1] > abs_cutoff]
    
    if absolute_passed or absolute_failed:
        plt.figure(figsize=(10, 6))
        all_rmsd = [s[1] for s in processed_scaffolds if s[1] != 99.9]
        colors = ['green' if r <= abs_cutoff else 'red' for r in all_rmsd]
        plt.scatter(range(len(all_rmsd)), all_rmsd, c=colors, alpha=0.5, s=10)
        plt.axhline(y=abs_cutoff, color='blue', linestyle='--', label=f'Cutoff ({abs_cutoff} Å)')
        plt.title(f'RFdiffusion RMSD Distribution (Run ID: {run_id})')
        plt.xlabel('Scaffold Index'); plt.ylabel('Pseudo-RMSD (from pLDDT)')
        plt.legend()
        plt.savefig(os.path.join(plot_dir, f"rfdiffusion_selection_dot_graph_{run_id}.png"))
        plt.close()
        
        record_list = [(x[0], x[1], 'Passed') for x in absolute_passed] + [(x[0], x[1], 'Failed') for x in absolute_failed]
        pd.DataFrame(record_list, columns=['Binder_ID', 'RMSD', 'Status']).to_csv(os.path.join(report_dir, f"all_scaffolds_rmsd_report_{run_id}.csv"), index=False)

    if not absolute_passed: trigger_fatal_error("[Phase 2]", "0 scaffolds passed RMSD cutoff.", "컷오프 완화 필요", engine_log)
    
    final_survivors = sorted(absolute_passed, key=lambda x: x[1])[:max(1, int(len(scaffolds) * (top_percent / 100.0)))]
    
    cur.execute("UPDATE Binders SET status = 'Failed_RF_RMSD' WHERE target_id = ? AND generation = 'DeNovo_RF'", (target_id,))
    for b_id, _ in final_survivors: cur.execute("UPDATE Binders SET status = 'Passed_RF_RMSD' WHERE binder_id = ?", (b_id,))
    conn.commit(); conn.close()

    pd.DataFrame(final_survivors, columns=['Binder_ID', 'RMSD']).to_csv(os.path.join(report_dir, f"passed_scaffolds_list_{run_id}.csv"), index=False)
    print(f"    [SUCCESS] Triage Complete. Data extracted to {scaffold_dir}")
    
    finalize_log(engine_log, start_time)

def finalize_log(log_path, start_time):
    end_time = datetime.now()
    duration = end_time - start_time
    with open(log_path, 'r') as f: lines = f.readlines()
    if len(lines) > 1: lines[1] = f"=== SCRIPT END: {end_time.strftime('%Y-%m-%d %H:%M:%S')} (Duration: {duration}) ===\n"
    with open(log_path, 'w') as f: f.writelines(lines)

def main(config_path):
    config = load_config(config_path)
    mode = config.get('pipeline_control', {}).get('step02_rfdiffusion', 'overwrite').lower()
    if mode == 'skip': sys.exit(0)
    print(f">>> [Phase 2] RFdiffusion Generation & Triage Started...")
    run_rfdiffusion_triage(config_path, mode)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); main(parser.parse_args().config)
