#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <iostream>

#include "model.h"

inline void writeVTU(
    const std::string& filename,
    const Grid& grid,
    const std::vector<double>& temperature
)
{
    if (temperature.size() != grid.nodes.size()) {
        std::cerr << "VTU error: temperature vector size mismatch\n";
        return;
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "VTU error: cannot open file\n";
        return;
    }

    const int numPoints = grid.nodes.size();
    const int numCells  = grid.elements.size();

    // ===== XML HEADER =====
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << numPoints
         << "\" NumberOfCells=\"" << numCells << "\">\n";

    // ===== POINT DATA =====
    file << "      <PointData Scalars=\"Temperature\">\n";
    file << "        <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n";
    for (double T : temperature)
        file << "          " << T << "\n";
    file << "        </DataArray>\n";
    file << "      </PointData>\n";

    // ===== CELL DATA (empty) =====
    file << "      <CellData>\n";
    file << "      </CellData>\n";

    // ===== POINTS =====
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& n : grid.nodes)
        file << "          " << n.x << " " << n.y << " 0.0\n";
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // ===== CELLS =====
    file << "      <Cells>\n";

    // -- connectivity
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& e : grid.elements) {
        file << "          ";
        for (int i = 0; i < 4; ++i)
            file << (e.node_ids[i] - 1) << " ";
        file << "\n";
    }
    file << "        </DataArray>\n";

    // -- offsets
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (int i = 0; i < numCells; ++i) {
        offset += 4;
        file << "          " << offset << "\n";
    }
    file << "        </DataArray>\n";

    // -- types
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < numCells; ++i)
        file << "          9\n";  // VTK_QUAD
    file << "        </DataArray>\n";

    file << "      </Cells>\n";

    // ===== FOOTER =====
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
}

inline void writePVD(
    const std::string& filename,
    int numSteps,
    double dt
)
{
    std::ofstream file(filename);
    if (!file.is_open()) return;

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <Collection>\n";

    for (int i = 0; i < numSteps; ++i) {
        file << "    <DataSet timestep=\""
             << i * dt
             << "\" file=\"step_"
             << std::setw(3) << std::setfill('0') << i
             << ".vtu\"/>\n";
    }

    file << "  </Collection>\n";
    file << "</VTKFile>\n";

    file.close();
}