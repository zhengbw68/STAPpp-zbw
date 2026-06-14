"""Generate convergence test .dat files - using exact Fortran fixed-width format.

Format (from reverse-engineering test_bar.dat and test_beam.dat):
  Control:       4I8  (cols 1-8, 9-16, 17-24, 25-32)
  Node:          I8,3I4,3F12.6  (8 + 4+4+4 + 12+12+12 = 56 chars)
  Load header:   2I8
  Load line:     I8,I4,F20.6
  Element type:  3I8
  Bar material:  I8,2F...
  Beam material: I8,10F...
  8H material:   I8,2F...
  Connectivity:  depends on element type
"""
import os

DATA = os.path.dirname(os.path.abspath(__file__))

def write_dat(path, lines):
    with open(path, 'w') as f:
        for line in lines:
            f.write(line.rstrip() + '\r\n')
        # blank line at end
        f.write('\r\n')
    print(f"  Created: {path}")

# ============================================================
# 8H Convergence - Block under self-weight with NxNxN mesh
# ============================================================
def gen_8h_convergence():
    for n in [1, 2, 4, 8]:
        nn = n + 1
        numnp = nn * nn * nn
        numeg = 1
        nlcase = 1
        modex = 1

        lines = [f"8H convergence - {n}x{n}x{n} mesh under self-weight"]
        lines.append("")

        # Control: 4I8
        lines.append(f"{numnp:8d}{numeg:8d}{nlcase:8d}{modex:8d}")
        lines.append("")

        E = 2e11
        nu = 0.3
        rho = 7800.0
        L = 1.0
        g = 10.0

        dx = L / n
        nid = 0
        node_xyz = {}
        node_bc = {}
        for iz in range(nn):
            z = iz * dx
            for iy in range(nn):
                y = iy * dx
                for ix in range(nn):
                    x = ix * dx
                    nid += 1
                    node_xyz[nid] = (x, y, z)
                    node_bc[nid] = (1, 1, 1) if iz == 0 else (0, 0, 0)

        # I8,3I4,3F12.6 = 56 chars per line
        for nid in range(1, numnp + 1):
            x, y, z = node_xyz[nid]
            b = node_bc[nid]
            lines.append(f"{nid:8d}{b[0]:4d}{b[1]:4d}{b[2]:4d}{x:12.6f}{y:12.6f}{z:12.6f}")
        lines.append("")

        # Body force - element-by-element consistent nodal forces
        vol_elem = dx**3
        elem_force_per_node = -rho * g * vol_elem / 8.0

        # Count how many elements share each node
        node_elem_count = {nid: 0 for nid in range(1, numnp + 1)}
        for iz in range(n):
            for iy in range(n):
                for ix in range(n):
                    for dz in [0, 1]:
                        for dy in [0, 1]:
                            for dxi in [0, 1]:
                                ni = (iz+dz)*nn*nn + (iy+dy)*nn + (ix+dxi) + 1
                                node_elem_count[ni] += 1

        # Apply forces only to non-fixed nodes (those with bcode[2] == 0)
        num_active_loads = sum(1 for nid in range(1, numnp+1)
                               if node_elem_count[nid] > 0 and node_bc[nid][2] == 0)
        lines.append(f"{1:8d}{num_active_loads:8d}")
        for nid in range(1, numnp + 1):
            if node_elem_count[nid] > 0 and node_bc[nid][2] == 0:
                force = elem_force_per_node * node_elem_count[nid]
                lines.append(f"{nid:8d}{3:4d}{force:20.6f}")

        # Element group: 3I8
        nel = n**3
        lines.append(f"{4:8d}{nel:8d}{1:8d}")

        # Material: I8, 2F...
        lines.append(f"{1:8d}{E:20.6f}{nu:14.6f}")

        # Connectivity: I8, 8I8, I8  (8 nodes + material set)
        eid = 0
        for iz in range(n):
            for iy in range(n):
                for ix in range(n):
                    eid += 1
                    n1 = iz*nn*nn + iy*nn + ix + 1
                    n2 = iz*nn*nn + iy*nn + (ix+1) + 1
                    n4 = iz*nn*nn + (iy+1)*nn + ix + 1
                    n3 = iz*nn*nn + (iy+1)*nn + (ix+1) + 1
                    n5 = (iz+1)*nn*nn + iy*nn + ix + 1
                    n6 = (iz+1)*nn*nn + iy*nn + (ix+1) + 1
                    n8 = (iz+1)*nn*nn + (iy+1)*nn + ix + 1
                    n7 = (iz+1)*nn*nn + (iy+1)*nn + (ix+1) + 1
                    lines.append(f"{eid:8d}{n1:8d}{n2:8d}{n3:8d}{n4:8d}{n5:8d}{n6:8d}{n7:8d}{n8:8d}{1:8d}")

        write_dat(os.path.join(DATA, 'test-8h', f'8h_convergence_{n}.dat'), lines)

# ============================================================
# Bar Convergence
# ============================================================
def gen_bar_convergence():
    L = 5.5
    E = 2e11
    A = 0.01
    P = 1000.0

    for n in [1, 2, 4, 8, 16]:
        numnp = n + 1
        lines = [f"Bar convergence - {n} elements, tip load"]
        lines.append("")

        lines.append(f"{numnp:8d}{1:8d}{1:8d}{1:8d}")
        lines.append("")

        dx = L / n
        for i in range(numnp):
            x = i * dx
            b = (1, 1, 1) if i == 0 else (0, 1, 1)
            lines.append(f"{i+1:8d}{b[0]:4d}{b[1]:4d}{b[2]:4d}{x:12.6f}{0:12.6f}{0:12.6f}")
        lines.append("")

        lines.append(f"{1:8d}{1:8d}")
        lines.append(f"{numnp:8d}{1:4d}{P:20.6f}")

        lines.append(f"{1:8d}{n:8d}{1:8d}")
        lines.append(f"{1:8d}{E:20.6f}{A:14.6f}")
        for i in range(n):
            lines.append(f"{i+1:8d}{i+1:8d}{i+2:8d}{1:8d}")

        write_dat(os.path.join(DATA, 'test-bar', f'bar_convergence_{n}.dat'), lines)

# ============================================================
# Beam Convergence
# ============================================================
def gen_beam_convergence():
    L = 1.0
    E = 2e11
    nu = 0.3
    sec_a = 0.3
    sec_b = 0.3
    t1 = 0.05
    t2 = 0.05
    t3 = 0.05
    t4 = 0.05
    n1 = 0.0
    n2 = 1.0
    n3 = 0.0
    P = -1000.0

    for n_el in [1, 2, 4, 8]:
        numnp = n_el + 1
        lines = [f"Beam convergence - {n_el} elements, tip load"]
        lines.append("")

        lines.append(f"{numnp:8d}{1:8d}{1:8d}{1:8d}")
        lines.append("")

        dx = L / n_el
        for i in range(numnp):
            x = i * dx
            if i == 0:
                bcode = (1, 1, 1, 1, 1, 1)
            else:
                bcode = (0, 0, 0, 0, 0, 0)
            lines.append(f"{i+1:8d}{bcode[0]:4d}{bcode[1]:4d}{bcode[2]:4d}{bcode[3]:4d}{bcode[4]:4d}{bcode[5]:4d}{x:12.6f}{0:12.6f}{0:12.6f}")
        lines.append("")

        lines.append(f"{1:8d}{1:8d}")
        lines.append(f"{numnp:8d}{3:4d}{P:20.6f}")

        lines.append(f"{5:8d}{n_el:8d}{1:8d}")
        lines.append(f"{1:8d}{E:20.6f}{nu:14.6f}{sec_a:14.6f}{sec_b:14.6f}{t1:14.6f}{t2:14.6f}{t3:14.6f}{t4:14.6f}{n1:14.6f}{n2:14.6f}{n3:14.6f}")
        for i in range(n_el):
            lines.append(f"{i+1:8d}{i+1:8d}{i+2:8d}{1:8d}")

        write_dat(os.path.join(DATA, 'test-beam', f'beam_convergence_{n_el}.dat'), lines)

if __name__ == '__main__':
    print("Generating convergence test files...")
    gen_8h_convergence()
    gen_bar_convergence()
    gen_beam_convergence()
    print("Done.")
