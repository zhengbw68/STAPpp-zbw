"""Extract displacement results from convergence .out files."""
import os
import re
import math

DATA = os.path.dirname(os.path.abspath(__file__))

def extract_displacements(outpath):
    """Extract nodal displacements from a .out file."""
    nodes = {}
    in_disps = False
    try:
        with open(outpath) as f:
            for line in f:
                if 'X-DISPLACEMENT' in line and 'Y-DISPLACEMENT' in line:
                    in_disps = True
                    continue
                if in_disps:
                    if not line.strip() or line.strip().startswith('-') or 'TIME' in line or 'TOTAL' in line or 'NUMBER' in line:
                        if 'TIME' in line or 'TOTAL' in line:
                            break
                        continue
                    try:
                        parts = line.split()
                        if len(parts) >= 7:
                            nid = int(parts[0])
                            dx, dy, dz = float(parts[1]), float(parts[2]), float(parts[3])
                            nodes[nid] = (dx, dy, dz)
                    except (ValueError, IndexError):
                        pass
    except Exception:
        pass
    return nodes

# ==============================
# BAR convergence results
# ==============================
print("=" * 60)
print("BAR (T3D2) CONVERGENCE - Tip load P=1000N, L=5.5m")
print("Exact: u_tip = PL/EA = 1000*5.5/(2E11*0.01) = 2.750e-06 m")
print("=" * 60)
print(f"{'N_elements':>10} {'u_tip_(m)':>14} {'Error':>14} {'Rate':>8}")
print("-" * 50)

prev_error = None
for n in [1, 2, 4, 8, 16]:
    outpath = os.path.join(DATA, 'test-bar', f'bar_convergence_{n}.out')
    nodes = extract_displacements(outpath)
    tip_nid = n + 1
    if tip_nid in nodes:
        u_tip = nodes[tip_nid][0]  # x-displacement
        error = abs(u_tip - 2.750e-6)
        if prev_error and prev_error > 0:
            rate = math.log(error / prev_error) / math.log(0.5)
        else:
            rate = float('nan')
        print(f"{n:>10} {u_tip:>14.3e} {error:>14.3e} {rate:>8.3f}")
        prev_error = error

# ==============================
# BEAM convergence results
# ==============================
print()
print("=" * 60)
print("BEAM (B31) CONVERGENCE - Cantilever, tip load P=-1000N, L=1.0m")
# I = (ab^3 - (a-2t)(b-2t)^3)/12 for box section
a, b = 0.3, 0.3
t1 = t2 = t3 = t4 = 0.05
I = (a*b**3 - (a-t1-t3)*(b-t2-t4)**3) / 12.0
E = 2e11
P = -1000.0
L = 1.0
exact_dz = P*L**3/(3*E*I)
exact_ry = P*L**2/(2*E*I)
print(f"Box section I = {I:.6e} m^4")
print(f"Exact: dz = {exact_dz:.3e} m, ry = {exact_ry:.3e} rad")
print("=" * 60)
print(f"{'N_elements':>10} {'dz_(m)':>14} {'Error_dz':>14} {'ry_(rad)':>14} {'Error_ry':>14}")
print("-" * 70)

for n in [1, 2, 4, 8]:
    outpath = os.path.join(DATA, 'test-beam', f'beam_convergence_{n}.out')
    nodes = extract_displacements(outpath)
    tip_nid = n + 1
    if tip_nid in nodes:
        dz = nodes[tip_nid][2]  # z-displacement
        error_dz = abs(dz - exact_dz)
        # Extract rotation too from full .out
        ry = 0.0
        error_ry = 0.0
        print(f"{n:>10} {dz:>14.3e} {error_dz:>14.3e} {'-':>14} {'-':>14}")

# ==============================
# 8H convergence — more detailed extraction
# ==============================
print()
print("=" * 60)
print("8H (C3D8R) CONVERGENCE - Block 1x1x1, self-weight, bottom fixed")
print("Exact: u_z,top = -7800*10*1^2/(2*2e11) = -1.950e-07 m")
print("Poisson lateral: u_x,top = 0.3*1.950e-07 = 5.850e-08 m")
print("=" * 60)
print(f"{'Mesh':>8} {'DOFs':>6} {'u_z,top':>14} {'u_x,corner':>14} {'Error_u_z':>14}")
print("-" * 60)

for n in [1, 2, 4, 8]:
    outpath = os.path.join(DATA, 'test-8h', f'8h_convergence_{n}.out')
    nodes = extract_displacements(outpath)
    nn = n + 1
    # Top corner node at (1, 1, 1) has highest node number
    top_corner = nn * nn * nn

    if top_corner in nodes:
        dz = nodes[top_corner][2]
        dx = nodes[top_corner][0]
        error = abs(dz - (-1.950e-7))
        neq = 0
        try:
            with open(outpath, encoding='utf-8', errors='replace') as f:
                for line in f:
                    m = re.search(r'NUMBER OF EQUATIONS.*NEQ\)\s*=\s*(\d+)', line)
                    if m:
                        neq = int(m.group(1))
                        break
        except:
            pass
        print(f"{n}x{n}x{n:>4} {neq:>6} {dz:>14.3e} {dx:>14.3e} {error:>14.3e}")

print()
print("NOTE: Bar and Beam elements give exact nodal displacements with any")
print("mesh (their shape functions are exact for the governing equation).")
print("Only 8H (linear hexahedra) shows convergence behavior for this problem.")
