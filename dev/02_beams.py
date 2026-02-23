import numpy as np
import matplotlib.pyplot as plt

# Parameters
footprint_diameter = 11        # Updated to 11m (ATLAS standard)
footprint_radius = footprint_diameter / 2
pulse_spacing = 0.7            # 0.7m spacing along-track
beam_separation = 90           # 90m between the pair
num_pulses = 120
icesat2_green = '#00FF00'

# Generate along-track positions
x_positions = np.arange(0, num_pulses * pulse_spacing, pulse_spacing)
last_x = x_positions[-1]

fig, ax = plt.subplots(figsize=(14, 8), facecolor='#111111')
ax.set_facecolor('#111111')

def draw_icesat2_beam(y_offset, alpha_val, label):
    for x in x_positions:
        circle = plt.Circle((x, y_offset), footprint_radius,
                            color=icesat2_green, alpha=alpha_val, lw=0)
        ax.add_patch(circle)
    ax.plot(x_positions, [y_offset]*len(x_positions), color=icesat2_green,
            linestyle='dotted', linewidth=0.5, alpha=0.5)
    ax.text(last_x + 10, y_offset, label, color='white', va='center', fontsize=10)

# Draw the pair
draw_icesat2_beam(0, 0.08, "Strong Beam")
draw_icesat2_beam(beam_separation, 0.03, "Weak Beam")

# --- ANNOTATIONS ---

# 1. Across-track distance (90m)
ax.annotate('', xy=(20, 0), xytext=(20, beam_separation),
            arrowprops=dict(arrowstyle='<->', color='yellow', lw=1.5))
ax.text(22, beam_separation/2, "90m Across-Track\nSeparation",
        color='yellow', va='center', fontweight='bold')

# 2. Along-track pulse spacing (0.7m) - Zoom callout at the start
ax.annotate('0.7m pulse interval',
            xy=(0.35, -2), xytext=(5, -25),
            color='cyan', arrowprops=dict(arrowstyle='->', color='cyan', connectionstyle="arc3,rad=.2"))
ax.plot([0, pulse_spacing], [-2, -2], color='cyan', lw=2)

# 3. Footprint Diameter (Moved to the end)
# Draw the "caliper" lines for the diameter
diag_y = 0
ax.plot([last_x - footprint_radius, last_x + footprint_radius], [diag_y, diag_y],
        color='#FF8C00', lw=2, label="Diameter")
# Tick marks at edges
ax.plot([last_x - footprint_radius, last_x - footprint_radius], [-2, 2], color='#FF8C00', lw=1.5)
ax.plot([last_x + footprint_radius, last_x + footprint_radius], [-2, 2], color='#FF8C00', lw=1.5)

ax.text(last_x, -8, "11m Diameter", color='#FF8C00', ha='center', fontweight='bold')

# Formatting
ax.set_aspect('equal')
ax.set_xlabel("Along-track distance (m)", color='white')
ax.set_ylabel("Across-track distance (m)", color='white')
ax.set_title("ICESat-2 ATLAS Geometry: Footprint & Beam Scales", color='white', fontsize=16, pad=20)

ax.tick_params(colors='white')
for spine in ax.spines.values():
    spine.set_edgecolor('#444444')

ax.set_xlim(-15, last_x + 60)
ax.set_ylim(-35, beam_separation + 25)

plt.tight_layout()
plt.show()