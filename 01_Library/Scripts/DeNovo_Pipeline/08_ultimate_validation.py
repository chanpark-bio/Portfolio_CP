import os, sys, argparse, yaml, sqlite3, subprocess, glob, json, random
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

def run_ultimate_validation(config_path, mode):
    start_time = datetime.now()
    run_id = start_time.strftime("RUN_%Y%m%d_%H%M%S")

    config = load_config(config_path)
    target_id = config.get('project', {}).get('target_pdb_id', '6M0J')
    project_name = config.get('project', {}).get('project_name', target_id)
    global_seed = config.get('project', {}).get('global_seed', 42)
    random.seed(global_seed)

    v_params = config.get('validation_params', {})
    min_plddt, min_iptm, max_md_rmsd = v_params.get('min_plddt', 80.0), v_params.get('min_iptm', 0.6), v_params.get('max_rmsd', 2.0)
    md_params = v_params.get('md_params', {})
    temp, press, em_steps, nvt_steps, npt_steps, md_steps, dt = md_params.get('temperature', 300), md_params.get('pressure', 1.0), md_params.get('em_steps', 500), md_params.get('nvt_steps', 500), md_params.get('npt_steps', 500), md_params.get('md_steps', 500), md_params.get('dt', 0.002)

    profile_name = v_params.get('md_topology_profile', 'standard_soluble_protein')
    ref_path = os.path.join("01_Library", "References", "MD_Topologies", f"{profile_name}.yaml")
    
    if os.path.exists(ref_path):
        with open(ref_path, 'r') as f: md_ref = yaml.safe_load(f)
    else:
        md_ref = {'force_field': 'amber03', 'water_model': 'spce', 'box_type': 'cubic', 'box_distance': 1.0, 'gmx_interactive_groups': {'rms_reference': 4, 'rms_target': 4, 'rmsf_target': 4}}

    overrides = config.get('advanced_overrides', {})
    
    ff = overrides.get('md_force_field', md_ref.get('force_field', 'amber03'))
    water = overrides.get('md_water_model', md_ref.get('water_model', 'spce'))
    box_t = overrides.get('md_box_type', md_ref.get('box_type', 'cubic'))
    box_d = overrides.get('md_box_distance', md_ref.get('box_distance', 1.0))
    
    grp_rms_ref = overrides.get('md_rms_reference_group', md_ref.get('gmx_interactive_groups', {}).get('rms_reference', 4))
    grp_rms_tgt = overrides.get('md_rms_target_group', md_ref.get('gmx_interactive_groups', {}).get('rms_target', 4))
    grp_rmsf_tgt = overrides.get('md_rmsf_target_group', md_ref.get('gmx_interactive_groups', {}).get('rmsf_target', 4))


    workspace_dir = config.get('paths', {}).get('base_dir', os.path.join("03_Workspace", project_name))
    db_path = config.get('paths', {}).get('db_path', os.path.join(workspace_dir, "00_Master_Control", "insilico_master.db"))
    current_run_dir = config.get('paths', {}).get('current_run_dir', workspace_dir)

    md_root_dir = os.path.join(current_run_dir, "08_Validation")
    log_dir = os.path.join(md_root_dir, "00_logs"); os.makedirs(log_dir, exist_ok=True)
    traj_dir = os.path.join(md_root_dir, "01_trajectories"); os.makedirs(traj_dir, exist_ok=True)
    plot_dir = os.path.join(md_root_dir, "02_dynamics_plots"); os.makedirs(plot_dir, exist_ok=True)
    report_dir = os.path.join(md_root_dir, "03_final_reports"); os.makedirs(report_dir, exist_ok=True)

    master_log = os.path.join(log_dir, f"08_validation_master_{run_id}.log")
    with open(master_log, 'w') as f:
        f.write(f"=== SCRIPT START: {start_time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        f.write(f"=== SCRIPT END: [RUNNING...] ===\n")
        f.write(f"=== GROMACS TOPOLOGY: FF={ff}, Water={water}, Box={box_t} ===\n")
        f.write(f"===============================\n\n")

    conn = sqlite3.connect(db_path); cur = conn.cursor()

    if mode == 'resume':
        cur.execute("SELECT COUNT(*) FROM Binders WHERE status = 'Passed_Validation' AND target_id = ?", (target_id,))
        if cur.fetchone()[0] > 0: 
            print(f"    [Mode: RESUME] Bypassing...")
            conn.close(); finalize_log(master_log, start_time); return

    if mode == 'overwrite':
        cur.execute("UPDATE Binders SET status = 'Passed_Maturation' WHERE status IN ('Passed_Validation', 'Failed_Val_CF', 'Failed_Val_MD') AND target_id = ?", (target_id,))
        conn.commit()

    cur.execute("SELECT b.binder_id, b.sequence, s.colabfold_plddt, s.iptm, s.md_rmsd FROM Binders b JOIN QC_Structural s ON b.binder_id = s.binder_id WHERE b.status = 'Passed_Maturation' AND b.target_id = ?", (target_id,))
    candidates = cur.fetchall()
    if not candidates: trigger_fatal_error("[Phase 5]", "No candidates ready.", "7호기 확인")

    colabfold_bin = os.environ.get("COLABFOLD_PATH", "colabfold_batch")
    gmx_bin = os.environ.get("GMX_PATH", "gmx")

    print(f"    [Validation] 🌪️ IGNITING SUPERCOMPUTING ENGINES for {len(candidates)} mutants (Run ID: {run_id})...")

    final_masters = []
    for b_id, seq, plddt, iptm, md_rmsd in candidates:
        mutant_dir = os.path.join(traj_dir, f"{b_id}_Validation")
        os.makedirs(mutant_dir, exist_ok=True)
        engine_log = os.path.join(log_dir, f"validation_engine_{b_id}_{run_id}.log") 

        if plddt is not None and iptm is not None: 
            real_plddt, real_iptm = plddt, iptm
        else:
            fasta_path = os.path.join(mutant_dir, f"{b_id}.fasta")
            with open(fasta_path, 'w') as f: f.write(f">{b_id}\n{seq}\n")
            try:
                subprocess.run([colabfold_bin, fasta_path, mutant_dir, "--num-models", "1", "--stop-at-score", "90"], check=True, capture_output=True)
                json_files = glob.glob(os.path.join(mutant_dir, "*_scores_rank_001_*.json"))
                if json_files:
                    scores = json.load(open(json_files[0], 'r'))
                    real_plddt = round(sum(scores['plddt']) / len(scores['plddt']), 2)
                    real_iptm = round(scores.get('iptm', scores.get('ptm', 0.0)), 3)
                else: real_plddt, real_iptm = 0.0, 0.0
            except subprocess.CalledProcessError: real_plddt, real_iptm = 0.0, 0.0
            cur.execute("UPDATE QC_Structural SET colabfold_plddt=?, iptm=? WHERE binder_id=?", (real_plddt, real_iptm, b_id))
            conn.commit()

        if real_plddt >= min_plddt and real_iptm >= min_iptm:
            if md_rmsd is not None and md_rmsd != 99.9: 
                real_md_rmsd = md_rmsd
            else:
                predicted_pdb = glob.glob(os.path.join(mutant_dir, "*_unrelaxed_rank_001_*.pdb"))
                if predicted_pdb:
                    predicted_pdb = predicted_pdb[0]
                    try:
                        with open(engine_log, 'a') as lf:
                            subprocess.run([gmx_bin, "pdb2gmx", "-f", predicted_pdb, "-o", f"{mutant_dir}/processed.gro", "-p", f"{mutant_dir}/topol.top", "-i", f"{mutant_dir}/posre.itp", "-water", water, "-ff", ff], check=True, stdout=lf, stderr=subprocess.STDOUT)
                            subprocess.run([gmx_bin, "editconf", "-f", f"{mutant_dir}/processed.gro", "-o", f"{mutant_dir}/box.gro", "-c", "-d", str(box_d), "-bt", box_t], check=True, stdout=lf, stderr=subprocess.STDOUT)
                            subprocess.run([gmx_bin, "solvate", "-cp", f"{mutant_dir}/box.gro", "-cs", "spc216.gro", "-o", f"{mutant_dir}/solvated.gro", "-p", f"{mutant_dir}/topol.top"], check=True, stdout=lf, stderr=subprocess.STDOUT)
                            
                            with open(os.path.join(mutant_dir, "em.mdp"), 'w') as mf: mf.write(f"integrator=steep\nemtol=1000.0\nnsteps={em_steps}\ncutoff-scheme=Verlet\ncoulombtype=PME\nrcoulomb=1.0\nrvdw=1.0\npbc=xyz\n")
                            subprocess.run([gmx_bin, "grompp", "-f", f"{mutant_dir}/em.mdp", "-c", f"{mutant_dir}/solvated.gro", "-p", f"{mutant_dir}/topol.top", "-o", f"{mutant_dir}/em.tpr", "-maxwarn", "2"], check=True, stdout=lf, stderr=subprocess.STDOUT)
                            subprocess.run([gmx_bin, "mdrun", "-v", "-deffnm", f"{mutant_dir}/em"], check=True, stdout=lf, stderr=subprocess.STDOUT)

                            with open(os.path.join(mutant_dir, "nvt.mdp"), 'w') as mf: mf.write(f"integrator=md\nnsteps={nvt_steps}\ndt={dt}\ncutoff-scheme=Verlet\ncoulombtype=PME\nrcoulomb=1.0\nrvdw=1.0\npbc=xyz\ntcoupl=V-rescale\ntc-grps=System\ntau_t=0.1\nref_t={temp}\ngen_vel=yes\ngen_temp={temp}\ngen_seed={global_seed}\n")
                            subprocess.run([gmx_bin, "grompp", "-f", f"{mutant_dir}/nvt.mdp", "-c", f"{mutant_dir}/em.gro", "-r", f"{mutant_dir}/em.gro", "-p", f"{mutant_dir}/topol.top", "-o", f"{mutant_dir}/nvt.tpr", "-maxwarn", "2"], check=True, stdout=lf, stderr=subprocess.STDOUT)
                            subprocess.run([gmx_bin, "mdrun", "-v", "-deffnm", f"{mutant_dir}/nvt"], check=True, stdout=lf, stderr=subprocess.STDOUT)

                            with open(os.path.join(mutant_dir, "npt.mdp"), 'w') as mf: mf.write(f"integrator=md\nnsteps={npt_steps}\ndt={dt}\ncutoff-scheme=Verlet\ncoulombtype=PME\nrcoulomb=1.0\nrvdw=1.0\npbc=xyz\ntcoupl=V-rescale\ntc-grps=System\ntau_t=0.1\nref_t={temp}\npcoupl=Parrinello-Rahman\npcoupltype=isotropic\ntau_p=2.0\nref_p={press}\ncompressibility=4.5e-5\n")
                            subprocess.run([gmx_bin, "grompp", "-f", f"{mutant_dir}/npt.mdp", "-c", f"{mutant_dir}/nvt.gro", "-r", f"{mutant_dir}/nvt.gro", "-t", f"{mutant_dir}/nvt.cpt", "-p", f"{mutant_dir}/topol.top", "-o", f"{mutant_dir}/npt.tpr", "-maxwarn", "2"], check=True, stdout=lf, stderr=subprocess.STDOUT)
                            subprocess.run([gmx_bin, "mdrun", "-v", "-deffnm", f"{mutant_dir}/npt"], check=True, stdout=lf, stderr=subprocess.STDOUT)

                            with open(os.path.join(mutant_dir, "md.mdp"), 'w') as mf: mf.write(f"integrator=md\nnsteps={md_steps}\ndt={dt}\nnstxout=100\nnstvout=100\ncutoff-scheme=Verlet\ncoulombtype=PME\nrcoulomb=1.0\nrvdw=1.0\npbc=xyz\ntcoupl=V-rescale\ntc-grps=System\ntau_t=0.1\nref_t={temp}\npcoupl=Parrinello-Rahman\npcoupltype=isotropic\ntau_p=2.0\nref_p={press}\ncompressibility=4.5e-5\n")
                            subprocess.run([gmx_bin, "grompp", "-f", f"{mutant_dir}/md.mdp", "-c", f"{mutant_dir}/npt.gro", "-t", f"{mutant_dir}/npt.cpt", "-p", f"{mutant_dir}/topol.top", "-o", f"{mutant_dir}/md.tpr", "-maxwarn", "2"], check=True, stdout=lf, stderr=subprocess.STDOUT)
                            subprocess.run([gmx_bin, "mdrun", "-v", "-deffnm", f"{mutant_dir}/md"], check=True, stdout=lf, stderr=subprocess.STDOUT)
                            
                            rms_process = subprocess.Popen([gmx_bin, "rms", "-s", f"{mutant_dir}/md.tpr", "-f", f"{mutant_dir}/md.trr", "-o", f"{mutant_dir}/rmsd.xvg"], stdin=subprocess.PIPE, stdout=lf, stderr=subprocess.STDOUT)
                            rms_process.communicate(input=f"{grp_rms_ref}\n{grp_rms_tgt}\n".encode())
                            
                            rmsf_process = subprocess.Popen([gmx_bin, "rmsf", "-s", f"{mutant_dir}/md.tpr", "-f", f"{mutant_dir}/md.trr", "-o", f"{mutant_dir}/rmsf.xvg", "-res"], stdin=subprocess.PIPE, stdout=lf, stderr=subprocess.STDOUT)
                            rmsf_process.communicate(input=f"{grp_rmsf_tgt}\n".encode())
                    except subprocess.CalledProcessError: pass
                            
                if os.path.exists(f"{mutant_dir}/rmsd.xvg"):
                    lines = [l for l in open(f"{mutant_dir}/rmsd.xvg", 'r').readlines() if not l.startswith(('@', '#'))]
                    real_md_rmsd = round(float(lines[-1].split()[1]) * 10.0, 2) if lines else 99.9
                else: real_md_rmsd = 99.9
                
                cur.execute("UPDATE QC_Structural SET md_rmsd=? WHERE binder_id=?", (real_md_rmsd, b_id))
                conn.commit()

            if real_md_rmsd <= max_md_rmsd:
                final_masters.append(b_id)
                cur.execute("UPDATE Binders SET status = 'Passed_Validation' WHERE binder_id = ?", (b_id,))
            else: cur.execute("UPDATE Binders SET status = 'Failed_Val_MD' WHERE binder_id = ?", (b_id,))

    conn.commit(); conn.close()
    if not final_masters: trigger_fatal_error("[Phase 5]", "All mutants collapsed.", "조건 완화 필요.")
    print(f"    [REPORT] 🛡️ Ultimate Validation Completed. Data extracted to {md_root_dir}")
    
    finalize_log(master_log, start_time)

def finalize_log(log_path, start_time):
    end_time = datetime.now()
    duration = end_time - start_time
    with open(log_path, 'r') as f: lines = f.readlines()
    if len(lines) > 1: lines[1] = f"=== SCRIPT END: {end_time.strftime('%Y-%m-%d %H:%M:%S')} (Duration: {duration}) ===\n"
    with open(log_path, 'w') as f: f.writelines(lines)

def main(config_path):
    config = load_config(config_path)
    mode = config.get('pipeline_control', {}).get('step08_validation', 'overwrite').lower()
    if mode == 'skip': sys.exit(0)
    print(f">>> [Phase 8] Real Ultimate Validation Started...")
    run_ultimate_validation(config_path, mode)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); main(parser.parse_args().config)
