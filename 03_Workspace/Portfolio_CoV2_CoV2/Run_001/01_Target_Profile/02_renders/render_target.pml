load 03_Workspace/Portfolio_CoV2_CoV2/Run_002/01_Target_Profile/01_raw_pdb/6M0J_raw.pdb, target
hide all; show cartoon, target
bg_color white; set ray_opaque_background, off
set cartoon_flat_sheets, 0
viewport 1920, 1080
zoom target, buffer=10.0
set antialias, 2
mset 1 x120
util.mroll 1, 120, 1
set cache_frames, 1
set ray_trace_frames, 0
color palecyan, chain A
color paleyellow, chain B
color lightpink, chain C
color palegreen, chain D
color gray70, (not chain A+B+C+D)
png 03_Workspace/Portfolio_CoV2_CoV2/Run_002/01_Target_Profile/02_renders/01_chains_colored.png, dpi=600, ray=1
movie.produce 03_Workspace/Portfolio_CoV2_CoV2/Run_002/01_Target_Profile/02_renders/01_chains_colored.mp4, quality=100
color gray80, target
select hotspots, (chain E and resi 484) or (chain E and resi 493) or (chain E and resi 501)
color tv_red, hotspots
png 03_Workspace/Portfolio_CoV2_CoV2/Run_002/01_Target_Profile/02_renders/02_hotspot_highlighted.png, dpi=600, ray=1
movie.produce 03_Workspace/Portfolio_CoV2_CoV2/Run_002/01_Target_Profile/02_renders/02_hotspot_highlighted.mp4, quality=100
alter target, b=0.0
alter target and resn ILE, b=4.5
alter target and resn VAL, b=4.2
alter target and resn LEU, b=3.8
alter target and resn PHE, b=2.8
alter target and resn CYS, b=2.5
alter target and resn MET, b=1.9
alter target and resn ALA, b=1.8
alter target and resn GLY, b=-0.4
alter target and resn THR, b=-0.7
alter target and resn SER, b=-0.8
alter target and resn TRP, b=-0.9
alter target and resn TYR, b=-1.3
alter target and resn PRO, b=-1.6
alter target and resn HIS, b=-3.2
alter target and resn GLU, b=-3.5
alter target and resn GLN, b=-3.5
alter target and resn ASP, b=-3.5
alter target and resn ASN, b=-3.5
alter target and resn LYS, b=-3.9
alter target and resn ARG, b=-4.5
sort
show surface, target
spectrum b, blue_white_red, target, minimum=-4.5, maximum=4.5
set transparency, 0.2
png 03_Workspace/Portfolio_CoV2_CoV2/Run_002/01_Target_Profile/02_renders/03_hydrophobicity_surface.png, dpi=600, ray=1
movie.produce 03_Workspace/Portfolio_CoV2_CoV2/Run_002/01_Target_Profile/02_renders/03_hydrophobicity_surface.mp4, quality=100
quit
