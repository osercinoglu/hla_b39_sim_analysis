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
select pair_residues, (chain M and resi 78) or (chain M and resi 95) or (chain M and resi 109) or (chain M and resi 130) or (chain B and resi 62) or (chain M and resi 96) or (chain M and resi 116) or (chain M and resi 123) or (chain M and resi 81) or (chain M and resi 118) or (chain M and resi 103) or (chain M and resi 110) or (chain M and resi 124) or (chain M and resi 135) or (chain M and resi 169) or (chain M and resi 173) or (chain M and resi 217) or (chain M and resi 246) or (chain M and resi 33) or (chain M and resi 52) or (chain M and resi 85) or (chain M and resi 58) or (chain M and resi 62) or (chain M and resi 125) or (chain M and resi 133) or (chain M and resi 127) or (chain M and resi 134) or (chain M and resi 99) or (chain M and resi 113) or (chain M and resi 10) or (chain B and resi 56) or (chain M and resi 261) or (chain M and resi 272) or (chain M and resi 108) or (chain M and resi 165) or (chain M and resi 49) or (chain M and resi 228) or (chain M and resi 248) or (chain M and resi 147) or (chain M and resi 32) or (chain M and resi 239) or (chain M and resi 188) or (chain M and resi 204) or (chain M and resi 114) or (chain M and resi 156) or (chain M and resi 227) or (chain M and resi 2) or (chain M and resi 105) or (chain M and resi 194) or (chain M and resi 168) or (chain M and resi 200) or (chain M and resi 247) or (chain M and resi 224) or (chain M and resi 126) or (chain M and resi 131) or (chain M and resi 153) or (chain M and resi 184) or (chain M and resi 267) or (chain M and resi 198) or (chain M and resi 11) or (chain M and resi 218) or (chain M and resi 260) or (chain M and resi 231) or (chain M and resi 244) or (chain M and resi 215) or (chain M and resi 104) or (chain M and resi 189) or (chain M and resi 55) or (chain M and resi 59) or (chain M and resi 25) or (chain M and resi 34) or (chain M and resi 102) or (chain M and resi 172) or (chain M and resi 30) or (chain M and resi 211) or (chain M and resi 82) or (chain M and resi 91) or (chain M and resi 5) or (chain B and resi 8) or (chain M and resi 234) or (chain M and resi 208) or (chain M and resi 242) or (chain M and resi 185) or (chain M and resi 191) or (chain M and resi 201) or (chain M and resi 35) or (chain B and resi 53) or (chain B and resi 10) or (chain B and resi 54) or (chain B and resi 60) or (chain M and resi 120) or (chain M and resi 229) or (chain M and resi 206) or (chain M and resi 98) or (chain M and resi 97) or (chain M and resi 115) or (chain M and resi 144) or (chain M and resi 112) or (chain M and resi 129) or (chain M and resi 137) or (chain M and resi 255) or (chain M and resi 273) or (chain M and resi 259) or (chain M and resi 28) or (chain M and resi 171) or (chain M and resi 187) or (chain M and resi 15) or (chain M and resi 262) or (chain M and resi 167) or (chain M and resi 199) or (chain M and resi 157) or (chain M and resi 111) or (chain M and resi 128) or (chain M and resi 236) or (chain M and resi 202) or (chain M and resi 13) or (chain M and resi 151) or (chain M and resi 155) or (chain M and resi 141) or (chain M and resi 94) or (chain M and resi 257) or (chain M and resi 93) or (chain M and resi 74) or (chain M and resi 243) or (chain M and resi 9) or (chain M and resi 24) or (chain M and resi 205) or (chain M and resi 76) or (chain M and resi 80) or (chain M and resi 258) or (chain M and resi 66) or (chain M and resi 70) or (chain M and resi 23) or (chain M and resi 4) or (chain M and resi 216) or (chain M and resi 274) or (chain M and resi 3) or (chain M and resi 230) or (chain M and resi 245) or (chain M and resi 26) or (chain M and resi 213) or (chain M and resi 84) or (chain M and resi 139) or (chain M and resi 117) or (chain M and resi 122) or (chain M and resi 192) or (chain M and resi 161) or (chain M and resi 176) or (chain M and resi 180) or (chain M and resi 72) or (chain M and resi 266) or (chain M and resi 270) or (chain M and resi 214) or (chain M and resi 87) or (chain M and resi 22) or (chain M and resi 38) or (chain M and resi 249) or (chain M and resi 145) or (chain M and resi 64) or (chain M and resi 68) or (chain M and resi 8) or (chain M and resi 136) or (chain M and resi 6) or (chain M and resi 29) or (chain M and resi 159) or (chain B and resi 12) or (chain M and resi 100) or (chain M and resi 36) or (chain M and resi 27) or (chain M and resi 19) or (chain M and resi 75) or (chain M and resi 160) or (chain M and resi 7) or (chain M and resi 20) or (chain M and resi 101) or (chain M and resi 271) or (chain M and resi 219) or (chain M and resi 142) or (chain M and resi 146) or (chain M and resi 254) or (chain M and resi 71) or (chain M and resi 140) or (chain M and resi 88) or (chain M and resi 203) or (chain M and resi 186) or (chain M and resi 67) or (chain M and resi 241) or (chain M and resi 31) or (chain M and resi 233) or (chain M and resi 14) or (chain M and resi 21) or (chain M and resi 39) or (chain M and resi 121) or (chain M and resi 164) or (chain M and resi 177) or (chain M and resi 181) or (chain B and resi 67) or (chain M and resi 237) or (chain M and resi 152) or (chain M and resi 143) or (chain M and resi 63) or (chain M and resi 12) or (chain M and resi 51) or (chain M and resi 190) or (chain B and resi 65) or (chain M and resi 235) or (chain M and resi 18) or (chain M and resi 174) or (chain M and resi 163) or (chain M and resi 46) or (chain M and resi 37) or (chain M and resi 77) or (chain M and resi 179) or (chain M and resi 148) or (chain M and resi 60) or (chain M and resi 269) or (chain M and resi 1) or (chain M and resi 61) or (chain M and resi 65) or (chain M and resi 89) or (chain M and resi 69) or (chain M and resi 250) or (chain M and resi 210) or (chain M and resi 263) or (chain M and resi 132)

# Show spheres and apply color
show spheres, pair_residues
color red, pair_residues
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