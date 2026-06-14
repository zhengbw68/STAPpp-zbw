
#include "PostOutputter.h"
#include "Domain.h"
#include "Elements/8H.h"
#include <iomanip>

const double coeff = 107.7;

PostOutputter* PostOutputter::_instance = nullptr;

PostOutputter::PostOutputter(string FileName)
{
    OutputFile.open(FileName);

    if (!OutputFile)
    {
        cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
        exit(3);
    }
}

// Return the single instance of the class
PostOutputter* PostOutputter::Instance(string FileName)
{
    if (!_instance)
        _instance = new PostOutputter(FileName);
    return _instance;
}

// ---------------------------------------------------------------------------
//  Helper: VTK cell type constants
//    VTK_HEXAHEDRON = 12  (8 nodes)
// ---------------------------------------------------------------------------
static int GetVTKCellType(ElementTypes type)
{
    switch (type)
    {
    case ElementTypes::Bar:          return 12;   // VTK_HEXAHEDRON
    case ElementTypes::Beam:         return 12;   // VTK_HEXAHEDRON
    case ElementTypes::Hexahedron:   return 12;   // VTK_HEXAHEDRON
    default:                         return 0;
    }
}

static int GetNumVizNodes(ElementTypes type)
{
    switch (type)
    {
    case ElementTypes::Bar:          return 8;
    case ElementTypes::Beam:         return 8;
    case ElementTypes::Hexahedron:   return 8;
    default:                         return 0;
    }
}

// ---------------------------------------------------------------------------
//  Stress invariants (same formulas as the original Tecplot output)
// ---------------------------------------------------------------------------
static double InvariantI1(const double s[6])
{
    return s[0] + s[1] + s[2];
}
static double InvariantI2(const double s[6])
{
    return s[0] * s[1] - s[3] * s[3]
         + s[0] * s[2] - s[5] * s[5]
         + s[1] * s[2] - s[4] * s[4];
}
static double InvariantI3(const double s[6])
{
    return s[0] * s[1] * s[2]
         + s[3] * s[4] * s[5] * 2.0
         - s[1] * s[5] * s[5]
         - s[2] * s[3] * s[3]
         - s[4] * s[4] * s[0];
}

// ---------------------------------------------------------------------------
//  Main: collect all element data, then write a single VTK Legacy Unstructured
//  Grid file that ParaView can read directly.
// ---------------------------------------------------------------------------
void PostOutputter::OutputElementStress()
{
    CDomain* FEMData = CDomain::Instance();

    double* Displacement = FEMData->GetDisplacement();

    const unsigned int NUMEG = FEMData->GetNUMEG(); // Number of element groups

    // -----------------------------------------------------------------------
    //  Data buffers — buffer all node / cell / stress data in one pass
    // -----------------------------------------------------------------------
    struct PointRec {
        double x, y, z;
        double sxx, syy, szz, sxy, syz, szx;
    };
    vector<PointRec> allPoints;
    vector<vector<int>> cellConnectivity;
    vector<int> cellTypes;

    // Loop for each element group
    for (unsigned int eg = 0; eg < NUMEG; eg++)
    {
        CElementGroup& group = FEMData->GetEleGrpList()[eg];

        ElementTypes elemType = group.GetElementType();
        unsigned int NUME = group.GetNUME();

        int nViz = GetNumVizNodes(elemType);
        int vtkCellType = GetVTKCellType(elemType);
        if (nViz == 0)      // unsupported element type → skip
            continue;

        switch (elemType)
        {
        // ===================================================================
        case ElementTypes::Bar:
        // ===================================================================
        {
            for (unsigned int ele = 0; ele < NUME; ele++)
            {
                CElement& elem = group.GetElement(ele);

                double stress[48];      // 8 nodes × 6 components
                double prePos[24];      // 8 nodes × 3 coordinates
                double postPos[24];     // 8 nodes × 3 coordinates

                elem.ElementPostInfo(stress, Displacement, prePos, postPos);

                vector<int> conn(static_cast<size_t>(nViz));
                for (int n = 0; n < nViz; n++)
                {
                    PointRec pt;
                    pt.x   = (1.0 - coeff) * prePos[3 * n + 0] + coeff * postPos[3 * n + 0];
                    pt.y   = (1.0 - coeff) * prePos[3 * n + 1] + coeff * postPos[3 * n + 1];
                    pt.z   = (1.0 - coeff) * prePos[3 * n + 2] + coeff * postPos[3 * n + 2];
                    pt.sxx = stress[6 * n + 0];
                    pt.syy = stress[6 * n + 1];
                    pt.szz = stress[6 * n + 2];
                    pt.sxy = stress[6 * n + 3];
                    pt.syz = stress[6 * n + 4];
                    pt.szx = stress[6 * n + 5];

                    conn[static_cast<size_t>(n)] = static_cast<int>(allPoints.size());
                    allPoints.push_back(pt);
                }
                cellConnectivity.push_back(conn);
                cellTypes.push_back(vtkCellType);
            }
            break;
        }

        // ===================================================================
        case ElementTypes::Beam:
        // ===================================================================
        {
            for (unsigned int ele = 0; ele < NUME; ele++)
            {
                CElement& elem = group.GetElement(ele);

                double stress[48];      // 8 nodes × 6 components
                double prePos[24];
                double postPos[24];

                elem.ElementPostInfo(stress, Displacement, prePos, postPos);

                vector<int> conn(static_cast<size_t>(nViz));
                for (int n = 0; n < nViz; n++)
                {
                    PointRec pt;
                    pt.x   = (1.0 - coeff) * prePos[3 * n + 0] + coeff * postPos[3 * n + 0];
                    pt.y   = (1.0 - coeff) * prePos[3 * n + 1] + coeff * postPos[3 * n + 1];
                    pt.z   = (1.0 - coeff) * prePos[3 * n + 2] + coeff * postPos[3 * n + 2];
                    pt.sxx = stress[6 * n + 0];
                    pt.syy = stress[6 * n + 1];
                    pt.szz = stress[6 * n + 2];
                    pt.sxy = stress[6 * n + 3];
                    pt.syz = stress[6 * n + 4];
                    pt.szx = stress[6 * n + 5];

                    conn[static_cast<size_t>(n)] = static_cast<int>(allPoints.size());
                    allPoints.push_back(pt);
                }
                cellConnectivity.push_back(conn);
                cellTypes.push_back(vtkCellType);
            }
            break;
        }

        // ===================================================================
        case ElementTypes::Hexahedron:
        // ===================================================================
        {
            // Temporary arrays for all hex elements in this group
            double* stressHex = new double[NUME * 48];
            double* prePos8H  = new double[NUME * 24];
            double* pos8H     = new double[NUME * 24];

            for (unsigned int ele = 0; ele < NUME; ele++)
            {
                group.GetElement(ele).ElementPostInfo(
                    &stressHex[48 * ele], Displacement,
                    &prePos8H[24 * ele], &pos8H[24 * ele]);
            }

            for (unsigned int ele = 0; ele < NUME; ele++)
            {
                vector<int> conn(static_cast<size_t>(nViz));
                for (int n = 0; n < nViz; n++)
                {
                    PointRec pt;
                    pt.x   = prePos8H[24 * ele + 3 * n + 0]
                           + coeff * (pos8H[24 * ele + 3 * n + 0] - prePos8H[24 * ele + 3 * n + 0]);
                    pt.y   = prePos8H[24 * ele + 3 * n + 1]
                           + coeff * (pos8H[24 * ele + 3 * n + 1] - prePos8H[24 * ele + 3 * n + 1]);
                    pt.z   = prePos8H[24 * ele + 3 * n + 2]
                           + coeff * (pos8H[24 * ele + 3 * n + 2] - prePos8H[24 * ele + 3 * n + 2]);
                    pt.sxx = stressHex[48 * ele + 6 * n + 0];
                    pt.syy = stressHex[48 * ele + 6 * n + 1];
                    pt.szz = stressHex[48 * ele + 6 * n + 2];
                    pt.sxy = stressHex[48 * ele + 6 * n + 3];
                    pt.syz = stressHex[48 * ele + 6 * n + 4];
                    pt.szx = stressHex[48 * ele + 6 * n + 5];

                    conn[static_cast<size_t>(n)] = static_cast<int>(allPoints.size());
                    allPoints.push_back(pt);
                }
                cellConnectivity.push_back(conn);
                cellTypes.push_back(vtkCellType);
            }

            delete[] stressHex;
            delete[] prePos8H;
            delete[] pos8H;
            break;
        }

        // ===================================================================
        default:
            cerr << "*** Error *** Element type " << static_cast<int>(elemType)
                 << " has not been implemented for post-processing.\n\n";
            break;
        } // end switch
    } // end for each group

    // ========================================================================
    //  Write VTK Legacy Unstructured Grid file
    // ========================================================================
    int nPoints = static_cast<int>(allPoints.size());
    int nCells  = static_cast<int>(cellConnectivity.size());

    // Compute total connectivity array size (CELLS section)
    int connSize = nCells;   // one leading "n" per cell
    for (const auto& conn : cellConnectivity)
        connSize += static_cast<int>(conn.size());

    // --- Header ---
    *this << "# vtk DataFile Version 3.0" << endl;
    *this << "STAPpp FEM Output" << endl;
    *this << "ASCII" << endl;
    *this << "DATASET UNSTRUCTURED_GRID" << endl;
    *this << endl;

    // --- POINTS ---
    *this << "POINTS " << nPoints << " double" << endl;
    for (const auto& pt : allPoints)
    {
        *this << scientific << setprecision(12)
              << pt.x << " " << pt.y << " " << pt.z << endl;
    }
    *this << endl;

    // --- CELLS ---
    *this << "CELLS " << nCells << " " << connSize << endl;
    for (const auto& conn : cellConnectivity)
    {
        *this << static_cast<int>(conn.size());
        for (int id : conn)
            *this << " " << id;
        *this << endl;
    }
    *this << endl;

    // --- CELL_TYPES ---
    *this << "CELL_TYPES " << nCells << endl;
    for (int type : cellTypes)
        *this << type << endl;
    *this << endl;

    // --- POINT_DATA ---
    *this << "POINT_DATA " << nPoints << endl;

    // Helper: pack a PointRec into a double[6] for invariant functions
    auto pack6 = [](const PointRec& p) -> const double* {
        static double s[6];
        s[0] = p.sxx; s[1] = p.syy; s[2] = p.szz;
        s[3] = p.sxy; s[4] = p.syz; s[5] = p.szx;
        return s;
    };

    // STRESS_I  (first invariant I₁ = σ_xx + σ_yy + σ_zz)
    *this << "SCALARS STRESS_I double 1" << endl;
    *this << "LOOKUP_TABLE default" << endl;
    for (const auto& pt : allPoints)
        *this << scientific << setprecision(12) << InvariantI1(pack6(pt)) << endl;

    // STRESS_II
    *this << "SCALARS STRESS_II double 1" << endl;
    *this << "LOOKUP_TABLE default" << endl;
    for (const auto& pt : allPoints)
        *this << scientific << setprecision(12) << InvariantI2(pack6(pt)) << endl;

    // STRESS_III
    *this << "SCALARS STRESS_III double 1" << endl;
    *this << "LOOKUP_TABLE default" << endl;
    for (const auto& pt : allPoints)
        *this << scientific << setprecision(12) << InvariantI3(pack6(pt)) << endl;

    // STRESS_VONMISES  (computed as sqrt(I₁² − I₂), matching original code)
    *this << "SCALARS STRESS_VONMISES double 1" << endl;
    *this << "LOOKUP_TABLE default" << endl;
    for (const auto& pt : allPoints)
    {
        const double* s = pack6(pt);
        double i1 = InvariantI1(s);
        double i2 = InvariantI2(s);
        *this << scientific << setprecision(12) << sqrt(i1 * i1 - i2) << endl;
    }

    // STRESS_XX
    *this << "SCALARS STRESS_XX double 1" << endl;
    *this << "LOOKUP_TABLE default" << endl;
    for (const auto& pt : allPoints)
        *this << scientific << setprecision(12) << pt.sxx << endl;

    // STRESS_YY
    *this << "SCALARS STRESS_YY double 1" << endl;
    *this << "LOOKUP_TABLE default" << endl;
    for (const auto& pt : allPoints)
        *this << scientific << setprecision(12) << pt.syy << endl;

    // STRESS_ZZ
    *this << "SCALARS STRESS_ZZ double 1" << endl;
    *this << "LOOKUP_TABLE default" << endl;
    for (const auto& pt : allPoints)
        *this << scientific << setprecision(12) << pt.szz << endl;

    // STRESS_XY
    *this << "SCALARS STRESS_XY double 1" << endl;
    *this << "LOOKUP_TABLE default" << endl;
    for (const auto& pt : allPoints)
        *this << scientific << setprecision(12) << pt.sxy << endl;

    // STRESS_YZ
    *this << "SCALARS STRESS_YZ double 1" << endl;
    *this << "LOOKUP_TABLE default" << endl;
    for (const auto& pt : allPoints)
        *this << scientific << setprecision(12) << pt.syz << endl;

    // STRESS_ZX
    *this << "SCALARS STRESS_ZX double 1" << endl;
    *this << "LOOKUP_TABLE default" << endl;
    for (const auto& pt : allPoints)
        *this << scientific << setprecision(12) << pt.szx << endl;
}
