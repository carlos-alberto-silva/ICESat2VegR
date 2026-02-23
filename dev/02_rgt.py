import numpy as np
import matplotlib.pyplot as plt

# Parameters for the Panorama View
# Units in Meters
swathe_width = 6600            # Total width (~6.6 km)
pair_separation = 3300         # Distance between pairs (3.3 km)
beam_separation = 90           # Distance within a pair (90 m)
along_track_length = 1000      # Show 1km of flight
icesat2_green = '#00FF00'

fig, ax = plt.subplots(figsize=(14, 8), facecolor='#0e0e0e')
ax.set_facecolor('#0e0e0e')

# Beam pair offsets relative to center (RGT)
# Pair 1 (Left), Pair 2 (Center/RGT), Pair 3 (Right)
pair_offsets = [-pair_separation, 0, pair_separation]

for i, offset in enumerate(pair_offsets):
    # Strong Beam (slightly thicker line/higher alpha)
    ax.plot([0, along_track_length], [offset, offset], 
            color=icesat2_green, alpha=0.8, lw=1.5, label=f"Pair {i+1} Strong" if i==0 else "")
    
    # Weak Beam (90m away from its strong partner)
    weak_y = offset + beam_separation
    ax.plot([0, along_track_length], [weak_y, weak_y], 
            color=icesat2_green, alpha=0.4, lw=1, label=f"Pair {i+1} Weak" if i==0 else "")

# --- PANORAMA ANNOTATIONS ---

# 1. The RGT (Reference Ground Track)
ax.axhline(0, color='white', linestyle='--', alpha=0.3, lw=1)
ax.text(along_track_length + 50, 0, "Reference Ground\nTrack (RGT)", 
        color='white', va='center', fontsize=9, style='italic')

# 2. Total Swathe Width (6.6 km)
ax.annotate('', xy=(-100, -pair_separation), xytext=(-100, pair_separation + beam_separation),
            arrowprops=dict(arrowstyle='<->', color='yellow', lw=2))
ax.text(-150, 0, "6.6 km Total Swathe", color='yellow', 
        ha='right', va='center', rotation=90, fontweight='bold')

# 3. Inter-pair spacing (3.3 km)
ax.annotate('', xy=(700, 0), xytext=(700, pair_separation),
            arrowprops=dict(arrowstyle='<->', color='cyan', lw=1.5))
ax.text(720, pair_separation/2, "3.3 km spacing", color='cyan', va='center')

# 4. Zoom-in callout for the 90m pair
ax.annotate('Each pair is 90m apart', 
            xy=(100, beam_separation), xytext=(100, 1000),
            color='lightgray', arrowprops=dict(arrowstyle='->', color='lightgray', connectionstyle="arc3,rad=-0.2"))

# Formatting
ax.set_xlabel("Along-track distance (m)", color='white')
ax.set_ylabel("Across-track distance (m)", color='white')
ax.set_title("ICESat-2 Global Beam Configuration (All 3 Pairs)", color='white', fontsize=16, pad=20)

# Set limits to show the full scale
ax.set_xlim(-250, 1250)
ax.set_ylim(-4000, 4500)

ax.tick_params(colors='white')
for spine in ax.spines.values():
    spine.set_edgecolor('#444444')

plt.tight_layout()
plt.show()