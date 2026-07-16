import os, sys, argparse, yaml, sqlite3, requests, subprocess
from datetime import datetime

def trigger_fatal_error(step, reason, action, log_path=None):
    err_snippet = ""
    log_info = f"\n 로그 경로 : {log_path}" if log_path else ""
    if log_path and os.path.exists(log_path):
        try:
            with open(log_path, 'r') as f:
                lines = f.readlines()
                if lines: err_snippet = "\n\n [원본 에러 메시지 (최근 10줄)]:\n" + "".join(lines[-10:])
        except Exception: pass
    print(f"\n==================================================\n [FATAL ERROR] {step} 가동 중단!\n 원인 : {reason}\n💡조치 : {action}{log_info}{err_snippet}\n==================================================\n")
    sys.exit(1)

def load_config(config_path):
    with open(config_path, 'r') as f: return yaml.safe_load(f)

def download_pdb(pdb_id, out_path):
    if os.path.exists(out_path):
        print(f"    [OK] Target PDB already exists. Skipping download.")
        return out_path
    print(f"    [RCSB] Downloading PDB {pdb_id}...")
    response_rcsb = requests.get(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb")
    if response_rcsb.status_code == 200:
        with open(out_path, 'w') as f: f.write(response_rcsb.text)
        return out_path
    else:
        response_af = requests.get(f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id.upper()}-F1-model_v4.pdb")
        if response_af.status_code == 200:
            with open(out_path, 'w') as f: f.write(response_af.text)
            return out_path
        else: trigger_fatal_error("[Phase 1]", f"Failed to fetch {pdb_id}.", "PDB ID 확인")

def run_structural_profiling(target_id, pdb_file, render_dir, config, engine_log):
    print(f"    [Profiling]  IGNITING PyMOL: Rendering  Hydrophobicity 3D maps & MP4s for {target_id}...")
    pymol_bin = os.environ.get("PYMOL_PATH", "pymol")
    hotspots = config.get('project', {}).get('hotspots', [])
    pml_path = os.path.join(render_dir, "render_target.pml")
    
    hotspot_str_parts = []
    for h in hotspots:
        if '-' in h: hotspot_str_parts.append(f"(chain {h[0]} and resi {h[1:]})")
        else: hotspot_str_parts.append(f"(chain {h[0]} and resi {h[1:]})")
    hotspot_str = " or ".join(hotspot_str_parts) if hotspot_str_parts else "none"

    ref_pml = "01_Library/References/eisenberg_scale.pml"
    if os.path.exists(ref_pml):
        with open(ref_pml, 'r') as rf: hydro_script = rf.read()
    else:
        hydro_script = "alter target, b=0.0\nsort\n"

    with open(pml_path, 'w') as f:
        f.write(f"load {pdb_file}, target\n")
        f.write(f"hide all; show cartoon, target\n")
        f.write(f"bg_color white; set ray_opaque_background, off\n")
        f.write(f"set cartoon_flat_sheets, 0\n")
        
        f.write(f"viewport 1920, 1080\n")
        f.write(f"zoom target, buffer=10.0\n")
        f.write(f"set antialias, 2\n")
        
        f.write(f"mset 1 x120\n")
        f.write(f"util.mroll 1, 120, 1\n")
        f.write(f"set cache_frames, 1\n") 
        f.write(f"set ray_trace_frames, 0\n") 
        
        f.write(f"color palecyan, chain A\n")
        f.write(f"color paleyellow, chain B\n")
        f.write(f"color lightpink, chain C\n")
        f.write(f"color palegreen, chain D\n")
        f.write(f"color gray70, (not chain A+B+C+D)\n")
        f.write(f"png {os.path.join(render_dir, '01_chains_colored.png')}, dpi=600, ray=1\n")
        f.write(f"movie.produce {os.path.join(render_dir, '01_chains_colored.mp4')}, quality=100\n")
        
        f.write(f"color gray80, target\n")
        if hotspots:
            f.write(f"select hotspots, {hotspot_str}\n")
            f.write(f"color tv_red, hotspots\n")
        f.write(f"png {os.path.join(render_dir, '02_hotspot_highlighted.png')}, dpi=600, ray=1\n")
        f.write(f"movie.produce {os.path.join(render_dir, '02_hotspot_highlighted.mp4')}, quality=100\n")
        
        f.write(hydro_script)
        f.write(f"show surface, target\n")
        f.write(f"spectrum b, blue_white_red, target, minimum=-4.5, maximum=4.5\n")
        f.write(f"set transparency, 0.2\n")
        f.write(f"png {os.path.join(render_dir, '03_hydrophobicity_surface.png')}, dpi=600, ray=1\n")
        f.write(f"movie.produce {os.path.join(render_dir, '03_hydrophobicity_surface.mp4')}, quality=100\n")
        f.write(f"quit\n")
        
    try:
        with open(engine_log, 'a') as lf:
            subprocess.run([pymol_bin, "-cq", pml_path], check=True, stdout=lf, stderr=subprocess.STDOUT)
            print("    [OK] High-Def PNGs and MP4 Movies rendered.")
    except Exception as e:
        with open(engine_log, 'a') as lf: lf.write(f"PyMOL Error: {e}\n")
        print("    [WARNING] PyMOL MP4 rendering failed (Check if FFmpeg is installed).")

def finalize_log(log_path, start_time):
    end_time = datetime.now()
    duration = end_time - start_time
    with open(log_path, 'r') as f: lines = f.readlines()
    if len(lines) > 1: lines[1] = f"=== SCRIPT END: {end_time.strftime('%Y-%m-%d %H:%M:%S')} (Duration: {duration}) ===\n"
    with open(log_path, 'w') as f: f.writelines(lines)

def main(config_path):
    start_time = datetime.now()
    run_id = start_time.strftime("RUN_%Y%m%d_%H%M%S")
    config = load_config(config_path)
    mode = config.get('pipeline_control', {}).get('step01_target_profiling', 'overwrite').lower()
    
    if mode == 'skip': 
        print(f">>> [Phase 1] Target Profiling Bypassed (Mode: SKIP)")
        sys.exit(0)
        
    print(f">>> [Phase 1] Target Profiling Started (Run ID: {run_id})...")
    target_id = config.get('project', {}).get('target_pdb_id', '6M0J')
    project_name = config.get('project', {}).get('project_name', target_id)
    global_seed = config.get('project', {}).get('global_seed', 42)
    
    workspace_dir = config.get('paths', {}).get('base_dir', os.path.join("03_Workspace", project_name))
    db_path = config.get('paths', {}).get('db_path', os.path.join(workspace_dir, "00_Master_Control", "insilico_master.db"))
    current_run_dir = config.get('paths', {}).get('current_run_dir', workspace_dir)
    
    target_dir = os.path.join(current_run_dir, "01_Target_Profile")
    log_dir = os.path.join(target_dir, "00_logs"); os.makedirs(log_dir, exist_ok=True)
    raw_dir = os.path.join(target_dir, "01_raw_pdb"); os.makedirs(raw_dir, exist_ok=True)
    render_dir = os.path.join(target_dir, "02_renders"); os.makedirs(render_dir, exist_ok=True)
    
    engine_log = os.path.join(log_dir, f"01_profiling_engine_{run_id}.log")
    with open(engine_log, 'w') as f:
        f.write(f"=== SCRIPT START: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"=== SCRIPT END: [RUNNING...]\n")
        f.write(f"=== RUN ID: {run_id} | GLOBAL SEED: {global_seed}\n")
        f.write(f"=== CONFIG INFO ===\n{yaml.dump(config.get('project', {}))}\n=================================\n\n")
    
    conn = sqlite3.connect(db_path)
    conn.execute("INSERT OR IGNORE INTO Targets (target_id) VALUES (?)", (target_id,))
    conn.commit(); conn.close()
    
    pdb_out_path = os.path.join(raw_dir, f"{target_id}_raw.pdb")
    download_pdb(target_id, pdb_out_path)
    run_structural_profiling(target_id, pdb_out_path, render_dir, config, engine_log)
    print(f"<<< [Phase 1] Completed.")
    finalize_log(engine_log, start_time)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); main(parser.parse_args().config)
