# PyMOL visualization script
# Generated from vis_pairs equivalent
# Includes Goodsell-style matte rendering settings

load /mnt/MORIA/experiments/hla_b39/2022_10_14_bioe_leviathan_backup/T1DM_B39_traj_analysis/gmx_traj_data/B3801_ALL_1/grinn_output_skip10/system_dry.pdb, structure

# ===========================================
# Publication-quality Goodsell-style settings
# ===========================================

# Background
set bg_color, white

# Ambient occlusion for soft shadows
set ambient_occlusion_mode, 1
set ambient_occlusion_scale, 25
set ambient_occlusion_smooth, 10

# Matte/flat shading (no specular highlights)
set specular, 0
set spec_reflect, 0
set shininess, 0

# Lighting for soft, diffuse look
set ambient, 0.4
set direct, 0.6
set reflect, 0.0

# Ray tracing settings (Goodsell-style)
set ray_trace_mode, 1
set ray_trace_gain, 0
set ray_trace_color, black
set ray_shadow, 1

# High quality rendering
set antialias, 2
set ray_trace_fog, 0
set depth_cue, 0

# Set up cartoon representation
hide all
show cartoon, structure
color gray80, structure
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1

# Select and highlight residue pairs
select pair_residues, (chain M and resi 217) or (chain M and resi 259)

# Show spheres and apply color
show spheres, pair_residues
color blue, pair_residues
set sphere_scale, 1.0, pair_residues
set sphere_mode, 0

# Zoom to highlighted residues
zoom pair_residues

# Deselect to clean up view
deselect

# ===========================================
# To render publication-quality image:
# ray 2400, 2400
# png output.png, dpi=300
# ===========================================