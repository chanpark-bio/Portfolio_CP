import os, argparse, yaml, sqlite3

def init_database(db_path):
    conn = sqlite3.connect(db_path); cur = conn.cursor()
    cur.execute('''CREATE TABLE IF NOT EXISTS Targets (target_id TEXT PRIMARY KEY)''')
    cur.execute('''CREATE TABLE IF NOT EXISTS Binders (binder_id TEXT PRIMARY KEY, target_id TEXT, sequence TEXT, status TEXT, parent_id TEXT, generation TEXT)''')
    cur.execute('''CREATE TABLE IF NOT EXISTS QC_Structural (binder_id TEXT PRIMARY KEY, rfd_rmsd REAL, colabfold_plddt REAL, iptm REAL, md_rmsd REAL)''')
    cur.execute('''CREATE TABLE IF NOT EXISTS QC_Physiochemical (binder_id TEXT PRIMARY KEY, instability REAL, mw REAL, pi REAL)''')
    cur.execute('''CREATE TABLE IF NOT EXISTS QC_Interaction (binder_id TEXT PRIMARY KEY, a3d_score REAL, immunogenicity_rank REAL, ddg REAL, score REAL)''')
    cur.execute('''CREATE TABLE IF NOT EXISTS QC_Production (binder_id TEXT PRIMARY KEY, cai_score REAL, gc_content REAL)''')
    conn.commit(); conn.close()

def create_project(project_name, target_name, pdb_id):
    workspace_name = f"{project_name}_{target_name}"
    base_dir = os.path.join("03_Workspace", workspace_name)
    
    if os.path.exists(base_dir):
        print(f" 이미 존재하는 프로젝트입니다: {base_dir}")
        return
        
    folders = [
        "00_Master_Control"
    ]
    for folder in folders: os.makedirs(os.path.join(base_dir, folder), exist_ok=True)

    config_path = os.path.join(base_dir, "00_Master_Control", "config.yaml")
    db_path = os.path.join(base_dir, "00_Master_Control", "insilico_master.db")
    
    config_data = {
        "project": {"project_name": workspace_name, "target_pdb_id": pdb_id, "global_seed": 42, "hotspots": ["A56", "A58", "A113", "A115"]},
        "paths": {"base_dir": base_dir, "control_dir": "00_Master_Control", "db_path": db_path},
        "pipeline_control": {
            "step01_target_profiling": "overwrite", "step02_rfdiffusion": "overwrite", "step03_proteinmpnn": "overwrite",
            "step04_fast_cascade": "overwrite", "step05_foldx_physics": "overwrite", "step06_champion_selector": "overwrite",
            "step07_maturation": "overwrite", "step08_validation": "overwrite", "step09_wetlab_translation": "overwrite"
        },
        "rfdiffusion": {"contig": "15-65", "num_designs": 10000},
        "filtering_strictness": {
            "rmsd_absolute_cutoff": 2.5, "rmsd_relative_top_percent": 1.0,
            "physicochemical_top_percent": 50.0, "max_instability": 40.0, "max_a3d_score": 1.0, "max_ddg_cutoff": -3.0, "ddg_top_percent": 20.0,
            "weights": {"ddg": 0.4, "a3d": 0.2, "immuno": 0.2, "instability": 0.2},
            "mutagenesis_library_size": 500, "mutagenesis_ddg_improve": -1.5, "absolute_elite_ddg": -3.0, "final_top_count": 3
        },
        "system": {"max_threads": 16, "host_organism": "e_coli_k12"},
        "validation_params": {"min_plddt": 80.0, "min_iptm": 0.6, "max_rmsd": 2.0, "md_topology_profile": "standard_soluble_protein", "md_params": {"temperature": 300, "pressure": 1.0, "em_steps": 500, "nvt_steps": 500, "npt_steps": 500, "md_steps": 500000, "dt": 0.002}}
    }
    
    with open(config_path, 'w') as f: yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)
    
    overrides_text = """
# =======================================================
# [Tier 2] 전문가용 강제 덮어쓰기 구역 (Advanced Overrides)
# 주의: 아래 항목들의 주석(#)을 해제하면 01_Library/References/ 의 기본 설정을 무시합니다.
# =======================================================
advanced_overrides: {}
  # --- 1. AI/딥러닝 엔진 설정 ---
  # mpnn_sampling_temp: 0.1
  
  # --- 2. 물리 엔진 설정 ---
  # foldx_temperature: "298K"
  
  # --- 3. MD 시뮬레이션 설정 (GROMACS) ---
  # md_force_field: "charmm36"
  # md_water_model: "tip3p"
  # md_box_type: "dodecahedron"
  # md_rms_reference_group: 4
  # md_rms_target_group: 4
  # md_rmsf_target_group: 4
  
  # --- 4. 면역원성 타겟 설정 (NetMHCpan) ---
  # netmhcpan_alleles: ["HLA-A*02:01", "HLA-A*24:02"]
"""
    with open(config_path, 'a') as f: f.write(overrides_text)
    
    init_database(db_path)
    print(f"    [CFG & DB]  프로젝트 환경(Tier 1 & 2) 및 장부 생성 완료")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True); parser.add_argument("--target", required=True); parser.add_argument("--pdb_id", required=True)
    args = parser.parse_args()
    create_project(args.project, args.target, args.pdb_id)
