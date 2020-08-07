#include "euler.h"

void write_results(euler &e, int step)
{

    std::string file_name;

    file_name = e.case_name + "_Asii." + std::to_string(step);

    e.file.open(file_name + ".vtu", std::fstream::out);
    e.file << std::fixed;
    e.file << std::setprecision(10);
    e.file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    e.file << "\t<UnstructuredGrid>\n";
    e.file << "\t\t<Piece NumberOfPoints=\"" << e.N_vertices << "\" NumberOfCells=\"" << e.N_cells << "\">\n";
    e.file << "\t\t\t<CellData>\n";
    e.file << "\t\t\t\t<DataArray type=\"Int32\" Name=\"type\" NumberOfComponents=\"1\" Format=\"ascii\">\n";
    for (int i = 0; i < e.N_cells; i++)
    {
        e.file << e.cells[i].type << "\n";
    }
    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"S\" NumberOfComponents=\"1\" Format=\"ascii\">\n";
    for (int i = 0; i < e.N_cells; i++)
    {
        e.file << e.cells[i].S << "\n";
    }
    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"rho\" NumberOfComponents=\"1\" Format=\"ascii\">\n";
    for (int i = 0; i < e.N_cells; i++)
    {
        e.file << e.cells[i].Q(0) << "\n";
    }
    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"u\" NumberOfComponents=\"1\" Format=\"ascii\">\n";
    for (int i = 0; i < e.N_cells; i++)
    {
        e.file << e.cells[i].Q(1) / e.cells[i].Q(0) << "\n";
    }
    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"1\" Format=\"ascii\">\n";
    for (int i = 0; i < e.N_cells; i++)
    {
        e.file << e.cells[i].Q(2) / e.cells[i].Q(0) << "\n";
    }
    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"p\" NumberOfComponents=\"1\" Format=\"ascii\">\n";
    for (int i = 0; i < e.N_cells; i++)
    {
        double rho, u, v, E, p;
        rho = e.cells[i].Q(0);
        u = e.cells[i].Q(1) / rho;
        v = e.cells[i].Q(2) / rho;
        E = e.cells[i].Q(3) / rho;
        p = (_gamma - 1.0) * rho * (E - 0.5 * (u * u + v * v));
        e.file << p << "\n";
    }
    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"dtau\" NumberOfComponents=\"1\" Format=\"ascii\">\n";
    for (int i = 0; i < e.N_cells; i++)
    {
        e.file << e.cells[i].dtau << "\n";
    }
    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"R1\" NumberOfComponents=\"1\" Format=\"ascii\">\n";
    for (int i = 0; i < e.N_cells; i++)
    {
        e.file << e.cells[i].R(0) << "\n";
    }
    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"R2\" NumberOfComponents=\"1\" Format=\"ascii\">\n";
    for (int i = 0; i < e.N_cells; i++)
    {
        e.file << e.cells[i].R(1) << "\n";
    }
    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"R3\" NumberOfComponents=\"1\" Format=\"ascii\">\n";
    for (int i = 0; i < e.N_cells; i++)
    {
        e.file << e.cells[i].R(2) << "\n";
    }
    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"R4\" NumberOfComponents=\"1\" Format=\"ascii\">\n";
    for (int i = 0; i < e.N_cells; i++)
    {
        e.file << e.cells[i].R(3) << "\n";
    }
    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t</CellData>\n";
    e.file << "\t\t\t<Points>\n";
    e.file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"position\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
    for (int i = 0; i < e.N_vertices; i++)
    {
        e.file << e.vertices[i](0) << "\t";
        e.file << e.vertices[i](1) << "\t";
        e.file << e.vertices[i](2) << "\n";
    }

    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t</Points>\n";
    e.file << "\t\t\t<Cells>\n";
    e.file << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";

    for (int i = 0; i < e.N_cells; i++)
    {
        int N_vertices = e.cells[i].vertices.size();
        for (int j = 0; j < N_vertices; j++)
        {
            if (j < (N_vertices - 1))
            {
                e.file << e.cells[i].vertices[j] << "\t";
            }
            else
            {
                e.file << e.cells[i].vertices[j] << "\n";
            }
        }
    }

    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";

    for (int i = 0; i < e.N_cells; i++)
    {
        int N_vertices = e.cells[i].vertices.size();
        e.file << N_vertices * (i + 1) << "\n";
    }

    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";

    for (int i = 0; i < e.N_cells; i++)
    {
        int N_vertices = e.cells[i].vertices.size();
        if (N_vertices == 3)
        {
            e.file << "5\n";
        }
        else
        {
            e.file << "9\n";
        }
    }

    e.file << "\t\t\t\t</DataArray>\n";
    e.file << "\t\t\t</Cells>\n";
    e.file << "\t\t</Piece>\n";
    e.file << "\t</UnstructuredGrid>\n";
    e.file << "</VTKFile>\n";
    e.file.close();
}