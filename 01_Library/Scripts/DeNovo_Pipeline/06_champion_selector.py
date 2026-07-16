import os, sys, argparse, yaml, sqlite3, json
from datetime import datetime
import pandas as pd, numpy as np
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

def normalize_series(series, reverse=False):
    s_min, s_max = series.min(), series.max()
    if s_min == s_max: return pd.Series([100.0] * len(series), index=series.index)
    if reverse: return 100.0 * (s_max - series) / (s_max - s_min + 1e-9)
    else: return 100.0 * (series - s_min) / (s_max - s_min + 1e-9)

def finalize_log(log_path, start_time):
    end_time = datetime.now()
    duration = end_time - start_time
    with open(log_path, 'r') as f: lines = f.readlines()
    if len(lines) > 1: lines[1] = f"=== SCRIPT END: {end_time.strftime('%Y-%m-%d %H:%M:%S')} (Duration: {duration}) ===\n"
    with open(log_path, 'w') as f: f.writelines(lines)

def run_champion_selector(config_path, mode):
    start_time = datetime.now()
    run_id = start_time.strftime("RUN_%Y%m%d_%H%M%S")

    config = load_config(config_path)
    target_id = config.get('project', {}).get('target_pdb_id', '6M0J')
    if not target_id: target_id = config.get('target', {}).get('pdb_id', '6M0J')
    project_name = config.get('project', {}).get('project_name', target_id)
    global_seed = config.get('project', {}).get('global_seed', 42)

    weights = config.get('filtering_strictness', {}).get('weights', {'ddg': 0.4, 'a3d': 0.2, 'immuno': 0.2, 'instability': 0.2})
    

    workspace_dir = config.get('paths', {}).get('base_dir', os.path.join("03_Workspace", project_name))
    db_path = config.get('paths', {}).get('db_path', os.path.join(workspace_dir, "00_Master_Control", "insilico_master.db"))
    current_run_dir = config.get('paths', {}).get('current_run_dir', workspace_dir)
    
    radar_dir = os.path.join(current_run_dir, "06_Champion_Radar")
    log_dir = os.path.join(radar_dir, "00_logs"); os.makedirs(log_dir, exist_ok=True)
    chart_dir = os.path.join(radar_dir, "01_radar_charts"); os.makedirs(chart_dir, exist_ok=True)
    data_dir = os.path.join(radar_dir, "02_champion_data"); os.makedirs(data_dir, exist_ok=True)
    
    engine_log = os.path.join(log_dir, f"06_champion_selector_{run_id}.log")
    with open(engine_log, 'w') as f:
        f.write(f"=== SCRIPT START: {start_time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        f.write(f"=== SCRIPT END: [RUNNING...] ===\n")
        f.write(f"=== RUN ID: {run_id} | GLOBAL SEED: {global_seed} ===\n")
        f.write(f"=== CONFIG INFO (Weights) ===\n{yaml.dump(weights)}\n=================================\n\n")

    conn = sqlite3.connect(db_path); cur = conn.cursor()

    if mode == 'resume':
        cur.execute("SELECT COUNT(*) FROM Binders WHERE status LIKE 'Passed_Radar%' AND target_id = ?", (target_id,))
        if cur.fetchone()[0] > 0: 
            print(f"    [Mode: RESUME] Champion already selected. Bypassing...")
            with open(engine_log, 'a') as lf: lf.write("Bypassed via RESUME mode.\n")
            conn.close(); finalize_log(engine_log, start_time); return

    if mode == 'overwrite':
        cur.execute("UPDATE Binders SET status = 'Passed_FoldX' WHERE status LIKE 'Passed_Radar%' AND target_id = ?", (target_id,))
        cur.execute("UPDATE QC_Interaction SET score = NULL WHERE binder_id IN (SELECT binder_id FROM Binders WHERE target_id = ?)", (target_id,))
        conn.commit()

    query = f"SELECT b.binder_id, p.instability, i.a3d_score, i.immunogenicity_rank, i.ddg FROM Binders b JOIN QC_Physiochemical p ON b.binder_id = p.binder_id JOIN QC_Interaction i ON b.binder_id = i.binder_id WHERE b.status = 'Passed_FoldX' AND b.target_id = '{target_id}'"
    df = pd.read_sql_query(query, conn)
    
    if df.empty: trigger_fatal_error("[Phase 3-5]", "No candidates survived up to 'Passed_FoldX'.", "컷오프 완화 필요")

    print(f"    [Ranking]  IGNITING DATA ANALYTICS for {len(df)} elite survivors (Run ID: {run_id})...")
    
    with open(engine_log, 'a') as lf:
        lf.write(f"\n--- [START] Champion Analytics ---\n")
        lf.write(f"Loaded {len(df)} candidates. Applied Weights: {weights}\n")

    df['norm_instability'] = normalize_series(df['instability'], reverse=True)
    df['norm_a3d'] = normalize_series(df['a3d_score'], reverse=True)
    df['norm_ddg'] = normalize_series(df['ddg'], reverse=True)
    df['norm_immuno'] = normalize_series(df['immunogenicity_rank'], reverse=False)

    df['total_score'] = (df['norm_ddg'] * weights.get('ddg', 0.4) + df['norm_a3d'] * weights.get('a3d', 0.2) + df['norm_immuno'] * weights.get('immuno', 0.2) + df['norm_instability'] * weights.get('instability', 0.2))
    df = df.sort_values(by='total_score', ascending=False).reset_index(drop=True)
    champion = df.iloc[0]
    
    for _, row in df.iterrows():
        cur.execute("UPDATE QC_Interaction SET score=? WHERE binder_id=?", (float(row['total_score']), row['binder_id']))
        status = 'Passed_Radar_Champion' if row['binder_id'] == champion['binder_id'] else 'Passed_Radar'
        cur.execute("UPDATE Binders SET status = ? WHERE binder_id = ?", (status, row['binder_id']))
    conn.commit(); conn.close()

    df.head(10).to_csv(os.path.join(data_dir, f"top_champions_ranking_{run_id}.csv"), index=False)
    print("    [Visuals] Rendering high-definition radar chart...")
    
    try:
        labels = ['Instability (Physical)', 'A3D (Aggregation)', 'FoldX (Affinity)', 'Immunogenicity (Stealth)']
        angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist(); angles += angles[:1]
        fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True))
        colors = ['#FF4B4B', '#4B8BFF', '#22C55E']
        for i in range(min(3, len(df))):
            row = df.iloc[i]
            values = [row['norm_instability'], row['norm_a3d'], row['norm_ddg'], row['norm_immuno']]; values += values[:1]
            ax.plot(angles, values, color=colors[i], linewidth=2.5, linestyle='solid', label=f"Rank {i+1}: {row['binder_id']}")
            ax.fill(angles, values, color=colors[i], alpha=0.15)
        ax.set_yticklabels([]); ax.set_xticks(angles[:-1]); ax.set_xticklabels(labels, fontsize=12, fontweight='bold', color='#333333')
        plt.legend(loc='upper right', bbox_to_anchor=(1.35, 1.15), fontsize=10)
        plt.savefig(os.path.join(chart_dir, f"final_champion_radar_chart_{run_id}.png"), dpi=300, bbox_inches='tight'); plt.close()
        
        with open(engine_log, 'a') as lf: lf.write("Radar chart rendered successfully.\n")
    except Exception as e: 
        with open(engine_log, 'a') as lf: lf.write(f"ERROR rendering chart: {e}\n")

    champ_dict = {"champion_id": champion['binder_id'], "total_score": round(champion['total_score'], 2), "ddg_kcal_mol": round(champion['ddg'], 2)}
    with open(os.path.join(data_dir, "champion_identity.json"), 'w') as f: json.dump(champ_dict, f, indent=4)
    print(f"    [SUCCESS] Champion identified: {champion['binder_id']} (Score: {champion['total_score']:.1f})")


    print("    [Visuals] Rendering Champion 3D Structure & Movie...")
    pymol_bin = os.environ.get("PYMOL_PATH", "pymol")
    champion_pdb = os.path.join(current_run_dir, "02_RF_Scaffolds", "01_generated_pdbs", f"{champion['binder_id'].split('_seq')[0]}.pdb")
    target_pdb = os.path.join(current_run_dir, "01_Target_Profile", "01_raw_pdb", f"{target_id}_raw.pdb")
    
    if os.path.exists(champion_pdb) and os.path.exists(target_pdb):
        champ_pml = os.path.join(chart_dir, "render_champion.pml")
        with open(champ_pml, 'w') as f:
            f.write(f"load {target_pdb}, target\n")
            f.write(f"load {champion_pdb}, champion\n")
            f.write(f"hide all; show cartoon\n")
            f.write(f"bg_color white; set ray_opaque_background, off\n")
            f.write(f"color gray80, target\n")
            f.write(f"color brightorange, champion\n")
            f.write(f"viewport 1920, 1080\n")
            f.write(f"zoom all, buffer=5.0\n")
            f.write(f"mset 1 x120; util.mroll 1, 120, 1; set cache_frames, 1; set ray_trace_frames, 0; set antialias, 2\n")
            f.write(f"png {os.path.join(chart_dir, f'champion_{champion['binder_id']}.png')}, dpi=600, ray=1\n")
            f.write(f"movie.produce {os.path.join(chart_dir, f'champion_{champion['binder_id']}.mp4')}, quality=100\n")
            f.write(f"quit\n")
        try:
            with open(engine_log, 'a') as lf:
                subprocess.run([pymol_bin, "-cq", champ_pml], check=True, stdout=lf, stderr=subprocess.STDOUT)
        except Exception: pass
    # ========================================================
    
    finalize_log(engine_log, start_time)

def main(config_path):
    config = load_config(config_path)
    mode = config.get('pipeline_control', {}).get('step06_champion_selector', 'overwrite').lower()
    if mode == 'skip': sys.exit(0)
    print(f">>> [Phase 6] Radar Scoring Started...")
    run_champion_selector(config_path, mode)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); main(parser.parse_args().config)
