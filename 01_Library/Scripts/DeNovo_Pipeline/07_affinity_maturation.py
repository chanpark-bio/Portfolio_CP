import os, sys, argparse, yaml, sqlite3, json, random, glob, subprocess, re
import pandas as pd
from datetime import datetime
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

def get_binder_chain_info(pdb_path, expected_seq_len):
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    chains = {}
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                chain_id, res_num, res_name = line[21], int(line[22:26].strip()), line[17:20].strip()
                if chain_id not in chains: chains[chain_id] = {}
                chains[chain_id][res_num] = d3to1.get(res_name, 'X')
    best_chain, min_diff = None, 9999
    for cid, seq_dict in chains.items():
        diff = abs(len(seq_dict) - expected_seq_len)
        if diff < min_diff: min_diff, best_chain = diff, cid
    if not best_chain: return "", 1, 'B'
    seq_dict = chains[best_chain]
    min_res, max_res = min(seq_dict.keys()), max(seq_dict.keys())
    return "".join([seq_dict.get(i, 'X') for i in range(min_res, max_res + 1)]), min_res, best_chain

def run_affinity_maturation(config_path, mode):
    start_time = datetime.now()
    run_id = start_time.strftime("RUN_%Y%m%d_%H%M%S")

    config = load_config(config_path)
    target_id = config.get('project', {}).get('target_pdb_id', '6M0J')
    project_name = config.get('project', {}).get('project_name', target_id)
    global_seed = config.get('project', {}).get('global_seed', 42)
    random.seed(global_seed)
    
    c = config.get('filtering_strictness', {})
    lib_size, ddg_improve, elite_ddg_threshold, top_count, max_instability = c.get('mutagenesis_library_size', 500), c.get('mutagenesis_ddg_improve', -1.5), c.get('absolute_elite_ddg', -3.0), c.get('final_top_count', 3), c.get('max_instability', 40.0)


    workspace_dir = config.get('paths', {}).get('base_dir', os.path.join("03_Workspace", project_name))
    db_path = config.get('paths', {}).get('db_path', os.path.join(workspace_dir, "00_Master_Control", "insilico_master.db"))
    current_run_dir = config.get('paths', {}).get('current_run_dir', workspace_dir)
    

    radar_dir = os.path.join(current_run_dir, "06_Champion_Radar")
    mut_dir = os.path.join(current_run_dir, "07_Maturation")
    log_dir = os.path.join(mut_dir, "00_logs"); os.makedirs(log_dir, exist_ok=True)
    pdb_dir = os.path.join(mut_dir, "01_mutant_pdbs"); os.makedirs(pdb_dir, exist_ok=True)
    report_dir = os.path.join(mut_dir, "02_evo_reports"); os.makedirs(report_dir, exist_ok=True)
    
    scaffold_src_dir = os.path.join(current_run_dir, "02_RF_Scaffolds", "01_generated_pdbs")

    engine_log = os.path.join(log_dir, f"maturation_engine_{run_id}.log")
    with open(engine_log, 'w') as f:
        f.write(f"=== SCRIPT START: {start_time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        f.write(f"=== SCRIPT END: [RUNNING...] ===\n")
        f.write(f"=== GLOBAL SEED: {global_seed} ===\n")
        f.write(f"=== CONFIG INFO (Filtering) ===\n{yaml.dump(c)}\n===============================\n\n")

    conn = sqlite3.connect(db_path); cur = conn.cursor()

    if mode == 'resume':
        cur.execute("SELECT COUNT(*) FROM Binders WHERE status = 'Passed_Maturation' AND target_id = ?", (target_id,))
        if cur.fetchone()[0] > 0: 
            print(f"    [Mode: RESUME] Bypassing...")
            with open(engine_log, 'a') as lf: lf.write("Bypassed via RESUME mode.\n")
            conn.close(); finalize_log(engine_log, start_time); return

    champ_json = os.path.join(radar_dir, "02_champion_data", "champion_identity.json")
    if not os.path.exists(champ_json): trigger_fatal_error("[Phase 4]", "champion JSON missing.", "6호기 확인")
    with open(champ_json, 'r') as f: champ_data = json.load(f)
    parent_id, parent_ddg = champ_data['champion_id'], champ_data['ddg_kcal_mol']

    cur.execute("SELECT b.sequence, p.instability FROM Binders b JOIN QC_Physiochemical p ON b.binder_id = p.binder_id WHERE b.binder_id = ?", (parent_id,))
    db_p_seq, p_instab = cur.fetchone()
    if p_instab is None: p_instab = 40.0

    if mode == 'overwrite':
        print("    [Mode: OVERWRITE] Soft Overwrite Active - 기존 돌연변이 보존.")
        cur.execute("UPDATE Binders SET status = 'Mutant_Generated' WHERE status IN ('Failed_Mut_ddG', 'Failed_Mut_Phys') AND parent_id = ? AND generation = 'Mutant_AlaScan'", (parent_id,))
        conn.commit()

    cur.execute("SELECT COUNT(*) FROM Binders WHERE parent_id = ? AND generation = 'Mutant_AlaScan'", (parent_id,))
    existing_mutants_count = cur.fetchone()[0]
    foldx_bin = os.environ.get("FOLDX_PATH", "foldx")
    
    if existing_mutants_count < lib_size:
        print(f" ->  IGNITING EVOLUTION: Generating {lib_size - existing_mutants_count} mutations...")
        parent_scaffold = parent_id.split('_seq')[0]
        parent_pdb_filename = f"{parent_scaffold}.pdb"
        pdb_seq, start_idx, chain_id = get_binder_chain_info(os.path.join(scaffold_src_dir, parent_pdb_filename), len(db_p_seq))

        with open(engine_log, 'a') as lf: lf.write(f"\n--- [START] FoldX BuildModel ---\n")
        for i in range(existing_mutants_count + 1, lib_size + 1):
            mut_id = f"{parent_id}_M{i:03d}"
            aas = "ACDEFGHIKLMNPQRSTVWY"
            seq_list = list(pdb_seq)
            valid_positions = [idx for idx, aa in enumerate(seq_list) if aa in aas]
            if not valid_positions: trigger_fatal_error("Mut", "No valid residues in Binder Chain", "PDB 확인")
            
            rel_pos = random.choice(valid_positions)
            orig_aa = seq_list[rel_pos]
            new_aa = random.choice(aas.replace(orig_aa, ""))
            seq_list[rel_pos] = new_aa
            mut_seq = "".join(seq_list)
            
            abs_pos = rel_pos + start_idx
            mut_filename = f"individual_list_{mut_id}.txt"
            with open(mut_filename, 'w') as mf: mf.write(f"{orig_aa}{chain_id}{abs_pos}{new_aa};\n")
            
            try:
                res = subprocess.run([foldx_bin, "--command=BuildModel", "--pdb", parent_pdb_filename, "--pdb-dir", scaffold_src_dir, "--mutant-file", mut_filename, "--output-dir", pdb_dir], capture_output=True, text=True)
                with open(engine_log, 'a') as lf: lf.write(f"[{mut_id} BuildModel]\n{res.stdout}\n{res.stderr}\n")
                foldx_out_pdb = os.path.join(pdb_dir, f"{parent_scaffold}_1.pdb")
                if os.path.exists(foldx_out_pdb): os.replace(foldx_out_pdb, os.path.join(pdb_dir, f"{mut_id}.pdb"))
            except Exception as e:
                with open(engine_log, 'a') as lf: lf.write(f"ERROR on {mut_id}: {e}\n")
            
            if os.path.exists(mut_filename): os.replace(mut_filename, os.path.join(pdb_dir, mut_filename))
            cur.execute("INSERT OR IGNORE INTO Binders (binder_id, target_id, sequence, status, parent_id, generation) VALUES (?, ?, ?, 'Mutant_Generated', ?, 'Mutant_AlaScan')", (mut_id, target_id, mut_seq, parent_id))
            for table in ['QC_Physiochemical', 'QC_Interaction', 'QC_Structural', 'QC_Production']: cur.execute(f"INSERT OR IGNORE INTO {table} (binder_id) VALUES (?)", (mut_id,))
        conn.commit()

    cur.execute("SELECT b.binder_id, b.sequence, p.instability, i.ddg FROM Binders b JOIN QC_Physiochemical p ON b.binder_id = p.binder_id JOIN QC_Interaction i ON b.binder_id = i.binder_id WHERE b.parent_id = ? AND b.generation = 'Mutant_AlaScan'", (parent_id,))
    mutants = cur.fetchall()
    survivors, calc_count = [], 0
    print(f" ->  EVALUATING MUTANTS: Running Cascade QC & Thermodynamic physics for {len(mutants)} variants...")
    
    with open(engine_log, 'a') as lf: lf.write(f"\n--- [START] FoldX AnalyseComplex (Mutants) ---\n")
    for m_id, m_seq, sim_instab, sim_ddg in mutants:
        if sim_ddg is not None and sim_ddg != 99.9:
            real_instab, real_ddg = sim_instab, sim_ddg
        else:
            real_instab = round(ProteinAnalysis(m_seq).instability_index(), 2) if len(m_seq) > 0 else 99.9
            mut_pdb_filename = f"{m_id}.pdb"
            try:
                if os.path.exists(os.path.join(pdb_dir, mut_pdb_filename)):
                    res = subprocess.run([foldx_bin, "--command=AnalyseComplex", "--pdb", mut_pdb_filename, "--pdb-dir", pdb_dir, "--output-dir", pdb_dir], capture_output=True, text=True)
                    matches = re.findall(r'Total\s*=\s*([-+]?\d*\.\d+|\d+)', res.stdout)
                    real_ddg = float(matches[-1]) if matches else 99.9
                else: real_ddg = 99.9
            except Exception: real_ddg = 99.9
            cur.execute("UPDATE QC_Physiochemical SET instability=? WHERE binder_id=?", (real_instab, m_id))
            cur.execute("UPDATE QC_Interaction SET ddg=? WHERE binder_id=?", (real_ddg, m_id))
            calc_count += 1
        
        relaxed_instab_limit = max(max_instability, p_instab + 5.0)
        if real_instab <= relaxed_instab_limit:
            if real_ddg <= (parent_ddg + ddg_improve) or real_ddg <= elite_ddg_threshold:
                survivors.append({"Mutant_ID": m_id, "ddG": real_ddg})
            else: cur.execute("UPDATE Binders SET status = 'Failed_Mut_ddG' WHERE binder_id = ?", (m_id,))
        else: cur.execute("UPDATE Binders SET status = 'Failed_Mut_Phys' WHERE binder_id = ?", (m_id,))

    if calc_count > 0: conn.commit()
    if not survivors: trigger_fatal_error("[Phase 4]", "No mutants passed.", "조건 완화 필요.")

    survivors.sort(key=lambda x: x['ddG']); top_mutants = survivors[:top_count]
    for mut in top_mutants: cur.execute("UPDATE Binders SET status = 'Passed_Maturation' WHERE binder_id = ?", (mut['Mutant_ID'],))
    conn.commit(); conn.close()
    pd.DataFrame(top_mutants).to_csv(os.path.join(report_dir, f"matured_top3_report_{run_id}.csv"), index=False)
    print(f"    [SUCCESS] Real Maturation complete. Top {len(top_mutants)} elite variants saved.")
    
    finalize_log(engine_log, start_time)

def finalize_log(log_path, start_time):
    end_time = datetime.now()
    duration = end_time - start_time
    with open(log_path, 'r') as f: lines = f.readlines()
    if len(lines) > 1:
        lines[1] = f"=== SCRIPT END: {end_time.strftime('%Y-%m-%d %H:%M:%S')} (Duration: {duration}) ===\n"
    with open(log_path, 'w') as f: f.writelines(lines)

def main(config_path):
    config = load_config(config_path)
    mode = config.get('pipeline_control', {}).get('step07_maturation', 'overwrite').lower()
    if mode == 'skip': sys.exit(0)
    print(f">>> [Phase 7] Real Affinity Maturation Started...")
    run_affinity_maturation(config_path, mode)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); main(parser.parse_args().config)
