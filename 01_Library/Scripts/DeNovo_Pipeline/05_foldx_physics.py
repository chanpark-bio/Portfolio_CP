import os, sys, argparse, yaml, sqlite3, json, subprocess, glob, re, random, math
from datetime import datetime
import pandas as pd

def fix_and_split_chains_3d(src_pdb, dst_pdb):
    with open(src_pdb, 'r') as f: lines = f.readlines()
    chains = {line[21] for line in lines if line.startswith("ATOM")}
    if len(chains) > 1:
        import shutil; shutil.copy(src_pdb, dst_pdb); return

    ca_coords = []
    for line in lines:
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_num = line[22:27]
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            ca_coords.append((res_num, x, y, z))
    
    split_res_str = None
    for i in range(1, len(ca_coords)):
        prev, curr = ca_coords[i-1], ca_coords[i]
        dist = math.sqrt((curr[1]-prev[1])**2 + (curr[2]-prev[2])**2 + (curr[3]-prev[3])**2)
        if dist > 8.0:  # 펩타이드 결합(3.8Å)을 아득히 넘는 물리적 단절 구간 발견!
            split_res_str = curr[0]
            break

    out_lines = []
    current_chain = 'A'
    passed_break = False
    
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            res_num = line[22:27]
            if split_res_str and not passed_break and res_num == split_res_str:
                out_lines.append("TER\n") 
                current_chain = 'B'
                passed_break = True
            elif split_res_str and passed_break:
                current_chain = 'B'
            new_line = line[:21] + current_chain + line[22:]
            out_lines.append(new_line)
        else:
            out_lines.append(line)
            
    with open(dst_pdb, 'w') as f: f.writelines(out_lines)

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

def run_foldx_physics(config_path, mode):
    start_time = datetime.now()
    run_id = start_time.strftime("RUN_%Y%m%d_%H%M%S")

    config = load_config(config_path)
    target_id = config.get('project', {}).get('target_pdb_id', '6M0J')
    project_name = config.get('project', {}).get('project_name', target_id)
    global_seed = config.get('project', {}).get('global_seed', 42)
    random.seed(global_seed)
    
    max_ddg_cutoff, top_percent = config.get('filtering_strictness', {}).get('max_ddg_cutoff', -3.0), config.get('filtering_strictness', {}).get('ddg_top_percent', 20.0)
    
    workspace_dir = config.get('paths', {}).get('base_dir', os.path.join("03_Workspace", project_name))
    db_path = config.get('paths', {}).get('db_path', os.path.join(workspace_dir, "00_Master_Control", "insilico_master.db"))
    current_run_dir = config.get('paths', {}).get('current_run_dir', workspace_dir)
    
    foldx_dir = os.path.join(current_run_dir, "05_FoldX_Physics")
    log_dir = os.path.join(foldx_dir, "00_logs"); os.makedirs(log_dir, exist_ok=True)
    link_dir = os.path.join(foldx_dir, "01_target_links"); os.makedirs(link_dir, exist_ok=True)
    report_dir = os.path.join(foldx_dir, "02_energy_reports"); os.makedirs(report_dir, exist_ok=True)
    
    engine_log = os.path.join(log_dir, f"foldx_engine_{run_id}.log")
    with open(engine_log, 'w') as f:
        f.write(f"=== SCRIPT START: {start_time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        f.write(f"=== SCRIPT END: [RUNNING...] ===\n")
        f.write(f"=== RUN ID: {run_id} | GLOBAL SEED: {global_seed} ===\n")
        f.write(f"=== CONFIG INFO (FoldX) ===\nCutoff: {max_ddg_cutoff}, Top: {top_percent}%\n===========================\n\n")

    conn = sqlite3.connect(db_path); cur = conn.cursor()

    if mode == 'resume':
        cur.execute("SELECT COUNT(*) FROM Binders WHERE status IN ('Passed_FoldX', 'Passed_Radar', 'Passed_Radar_Champion', 'Passed_Maturation', 'Mutant_Generated') AND target_id = ?", (target_id,))
        if cur.fetchone()[0] > 0: 
            print(f"    [Mode: RESUME] Bypassing...")
            with open(engine_log, 'a') as lf: lf.write("Bypassed via RESUME mode.\n")
            conn.close(); finalize_log(engine_log, start_time); return

    if mode == 'overwrite':
        cur.execute("UPDATE Binders SET status = 'Passed_FastQC' WHERE status LIKE 'Failed_FoldX%' AND target_id = ?", (target_id,))
        conn.commit()

    cur.execute("SELECT b.binder_id, b.parent_id, i.ddg FROM Binders b JOIN QC_Interaction i ON b.binder_id = i.binder_id WHERE b.status IN ('Passed_FastQC', 'Passed_FoldX', 'Passed_Radar', 'Passed_Radar_Champion') AND b.target_id = ?", (target_id,))
    candidates = cur.fetchall()
    if not candidates: trigger_fatal_error("[Phase 5]", "No candidates found.", "이전 단계 점검")

    print(f"    [FoldX]  EVALUATING PHYSICS for {len(candidates)} binders (Run ID: {run_id})...")
    foldx_bin = os.environ.get("FOLDX_PATH", "foldx")
    

    cur.execute("SELECT COUNT(*) FROM Binders WHERE status IN ('Passed_FastQC', 'Passed_FoldX', 'Failed_FoldX_Abs', 'Failed_FoldX_Rel', 'Passed_Radar', 'Passed_Radar_Champion', 'Passed_Maturation', 'Mutant_Generated') AND target_id = ?", (target_id,))
    total_foldx_pool = cur.fetchone()[0]
    already_done_foldx = total_foldx_pool - len(candidates)

    absolute_passed, calc_count, bypass_count, cache_hit_count = [], 0, 0, 0
    all_results = [] 
    scaffold_cache = {}
    
    for idx, (b_id, parent_id, existing_ddg) in enumerate(candidates):
        print(f"\r      -> [FoldX Engine] Calculating Thermodynamics: {already_done_foldx + idx + 1} / {total_foldx_pool}...", end="", flush=True)
        if existing_ddg is not None and existing_ddg != 99.9:
            real_ddg, bypass_count = existing_ddg, bypass_count + 1
            scaffold_cache[parent_id] = real_ddg 
        elif parent_id in scaffold_cache:
            real_ddg, cache_hit_count = scaffold_cache[parent_id], cache_hit_count + 1
        else:
            src_pdb = os.path.join(current_run_dir, "02_RF_Scaffolds", "01_generated_pdbs", f"{parent_id}.pdb")
            link_pdb = os.path.join(link_dir, f"{parent_id}.pdb")
            if os.path.exists(src_pdb) and not os.path.exists(link_pdb): 
                fix_and_split_chains_3d(src_pdb, link_pdb)

            pdb_filename = f"{parent_id}.pdb" 
            try:
                repair_cmd = [foldx_bin, "--command=RepairPDB", "--pdb", pdb_filename, "--pdb-dir", link_dir, "--output-dir", foldx_dir]
                subprocess.run(repair_cmd, capture_output=True, text=True)
                repaired_pdb = f"{parent_id}_Repair.pdb"
                
                analyse_cmd = [foldx_bin, "--command=AnalyseComplex", "--pdb", repaired_pdb, "--pdb-dir", foldx_dir, "--output-dir", foldx_dir]
                res = subprocess.run(analyse_cmd, capture_output=True, text=True)
                
                with open(engine_log, 'a') as lf: lf.write(f"[{b_id} Repair & Analyse]\n{res.stdout}\n")
                
                matches = re.findall(r'Total\s*=\s*([-+]?\d*\.\d+|\d+)', res.stdout)
                real_ddg = float(matches[-1]) if matches else 99.9
            except Exception as e: 
                with open(engine_log, 'a') as lf: lf.write(f"[{b_id} ERROR] 물리 엔진 가동 실패: {str(e)}\n")
                real_ddg = 99.9
            scaffold_cache[parent_id] = real_ddg
            calc_count += 1

        cur.execute("UPDATE QC_Interaction SET ddg=? WHERE binder_id=?", (real_ddg, b_id))
        all_results.append((b_id, parent_id, real_ddg)) 

    if len(candidates) > 0: print()
    if calc_count > 0 or cache_hit_count > 0: conn.commit()


    df = pd.DataFrame(all_results, columns=['Binder_ID', 'Parent_Scaffold', 'ddG_kcal_mol'])
    df = df.sort_values('ddG_kcal_mol')
    absolute_passed_df = df[df['ddG_kcal_mol'] <= max_ddg_cutoff]
    survivor_count = int(max(1, len(absolute_passed_df) * (top_percent / 100.0))) if not absolute_passed_df.empty else 0
    df['Status'] = ['Passed_FoldX' if i < survivor_count else 'Failed_FoldX_Rel' for i in range(len(df))]
    
    report_path = os.path.join(report_dir, f"foldx_binding_energy_results_{run_id}.csv")
    df.to_csv(report_path, index=False)
    print(f"    [FoldX] 📄 Energy Report Saved: {report_path}")


    if not absolute_passed_df.empty:
        for b_id in absolute_passed_df['Binder_ID'].iloc[:survivor_count]:
            cur.execute("UPDATE Binders SET status = 'Passed_FoldX' WHERE binder_id = ?", (b_id,))
        for b_id in absolute_passed_df['Binder_ID'].iloc[survivor_count:]:
            cur.execute("UPDATE Binders SET status = 'Failed_FoldX_Rel' WHERE binder_id = ?", (b_id,))
        failed_abs = df[df['ddG_kcal_mol'] > max_ddg_cutoff]
        for b_id in failed_abs['Binder_ID']:
            cur.execute("UPDATE Binders SET status = 'Failed_FoldX_Abs' WHERE binder_id = ?", (b_id,))
    else:
        for b_id in df['Binder_ID']:
            cur.execute("UPDATE Binders SET status = 'Failed_FoldX_Abs' WHERE binder_id = ?", (b_id,))
        conn.commit(); conn.close()
        trigger_fatal_error("[Phase 5]", "0 binders passed.", "컷오프 완화 필요", engine_log)

    conn.commit(); conn.close()
    print(f"    [SUCCESS] Physics Engine Complete. {survivor_count} survivors advance.")
    
    finalize_log(engine_log, start_time)

def finalize_log(log_path, start_time):
    end_time = datetime.now()
    duration = end_time - start_time
    with open(log_path, 'r') as f: lines = f.readlines()
    if len(lines) > 1: lines[1] = f"=== SCRIPT END: {end_time.strftime('%Y-%m-%d %H:%M:%S')} (Duration: {duration}) ===\n"
    with open(log_path, 'w') as f: f.writelines(lines)

def main(config_path):
    config = load_config(config_path)
    mode = config.get('pipeline_control', {}).get('step05_foldx_physics', 'overwrite').lower()
    if mode == 'skip': sys.exit(0)
    print(f">>> [Phase 5] FoldX Physics Engine Started...")
    run_foldx_physics(config_path, mode)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); main(parser.parse_args().config)
