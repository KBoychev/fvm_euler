#include "euler.h"

//Read grid
void read_grid(euler &e)
{
    // Read stl and remove duplicate vertices
    ////////////////////////////////////////////////////////////////////

    std::string file_line;
    std::stringstream file_line_stream;

    e.file.open(e.case_name + ".grd", std::fstream::in);

    int j = 0;
    if (e.file.is_open())
    {
        std::getline(e.file, file_line);
        file_line_stream.clear();
        file_line_stream.str(file_line);
        file_line_stream >> e.N_vertices >> e.N_cells;

        e.vertices.resize(e.N_vertices);

        for (int i = 0; i < e.N_vertices; i++)
        {
            Eigen::Vector3d vertex;

            std::getline(e.file, file_line);
            file_line_stream.clear();
            file_line_stream.str(file_line);
            file_line_stream >> vertex(0) >> vertex(1) >> vertex(2);
            e.vertices[i] = vertex;
        }

        e.cells.resize(e.N_cells);

        for (int i = 0; i < e.N_cells; i++)
        {
            int vertex;
            std::vector<int> vertices;

            std::getline(e.file, file_line);
            file_line_stream.clear();
            file_line_stream.str(file_line);
            while (file_line_stream >> vertex)
            {
                if (vertex < 0)
                {
                    std::cout << "Cell " << i << " vertex id is negative! Vertex ids start from 0!" << std::endl;
                    std::exit(-1);
                }
                vertices.push_back(vertex);
            }
            if (vertices.size() > 2)
            {
                cell cell(vertices);
                e.cells[i] = cell;
            }
            else
            {
                std::cout << "Cell " << i << " has only 2 vertices. Check grd file." << std::endl;
                std::exit(-1);
            }
        }

        e.file.close();

        e.N_vertices = e.vertices.size();
        e.N_cells = e.cells.size();

        // Create edges
        ////////////////////////////////////////////////////////////////////

        for (int i = 0; i < e.N_cells; i++)
        {
            int N_vertices = e.cells[i].vertices.size();

            for (int j = 0; j < N_vertices; j++)
            {
                if (j == (N_vertices - 1))
                {
                    edge edge(e.cells[i].vertices[j], e.cells[i].vertices[0], i, -1);
                    e.edges.push_back(edge);
                }
                else
                {
                    edge edge(e.cells[i].vertices[j], e.cells[i].vertices[j + 1], i, -1);
                    e.edges.push_back(edge);
                }
            }
        }

        e.N_edges = e.edges.size();

        //Remove duplicate edges
        ////////////////////////////////////////////////////////////////////

        for (int i = 0; i < e.N_edges; i++)
        {
            for (int j = 0; j < e.N_edges; j++)
            {
                if (e.edges[i].unique == 1)
                {
                    if (e.edges[i].vertex1 == e.edges[j].vertex2 && e.edges[i].vertex2 == e.edges[j].vertex1)
                    {
                        //Edge j is duplicate of edge i
                        e.edges[i].cellr = e.edges[j].celll;
                        e.edges[j].unique = 0;
                        break;
                    }
                }
            }
        }

        e.edges.erase(std::remove_if(e.edges.begin(), e.edges.end(), [](const edge &edge) {
                          return edge.unique == 0;
                      }),
                      e.edges.end());

        e.N_edges = e.edges.size();

        // Calculate cell normals, centroids, and area
        ////////////////////////////////////////////////////////////////////
        for (int i = 0; i < e.N_cells; i++)
        {
            int N_vertices;
            Eigen::Vector3d r;
            std::vector<Eigen::Vector3d> ri;

            N_vertices = e.cells[i].vertices.size();
            r = Eigen::Vector3d::Zero();

            for (int j = 0; j < N_vertices; j++)
            {
                r += e.vertices[e.cells[i].vertices[j]];
                ri.push_back(e.vertices[e.cells[i].vertices[j]]);
            }

            r /= N_vertices;
            e.cells[i].r = r;

            double S;
            Eigen::Vector3d d1, d2, n;

            if (N_vertices == 3)
            {
                d1 = ri[1] - ri[0];
                d2 = ri[2] - ri[0];
                n = d1.cross(d2);
            }

            if (N_vertices == 4)
            {
                d1 = ri[2] - ri[0];
                d2 = ri[3] - ri[1];
                n = d1.cross(d2);
            }
            S = 0.5 * n.norm();
            if (std::isnan(S) || S <= 0)
            {
                std::cout << "Cell " << i << " area is nan or negative!" << std::endl;
                std::cout << "Cell vertices = ";
                for (int j = 0; j < N_vertices; j++)
                {
                    std::cout << e.cells[i].vertices[j] << ", ";
                }
                std::cout << std::endl;
                std::exit(-1);
            }
            e.cells[i].n = n / n.norm();
            e.cells[i].S = S;
        }

        // Calculate edge normals, centroids, and length
        ////////////////////////////////////////////////////////////////////

        for (int i = 0; i < e.N_edges; i++)
        {

            int celll, cellr;
            Eigen::Vector3d r1, r2, d, n;

            celll = e.edges[i].celll;
            cellr = e.edges[i].cellr;
            r1 = e.vertices[e.edges[i].vertex1];
            r2 = e.vertices[e.edges[i].vertex2];
            d = r2 - r1;
            n = d.cross(e.cells[celll].n);
            e.edges[i].n = n / n.norm();
            e.edges[i].l = d.norm();
            e.edges[i].r = 0.5 * (r1 + r2);
        }

        // Mark boundaries
        ////////////////////////////////////////////////////////////////////
        for (int i = 0; i < e.N_edges; i++)
        {
            int celll, cellr;
            celll = e.edges[i].celll;
            cellr = e.edges[i].cellr;
            if (cellr == -1)
            {
                if (e.cells[celll].S > 0.01)
                {
                    e.cells[celll].type = 1;
                }
                else
                {
                    e.cells[celll].type = 2;
                }
            }
        }

        // Mark boundaries
        ////////////////////////////////////////////////////////////////////
        double S = 0.0;
        for (int i = 0; i < e.N_cells; i++)
        {
            S += e.cells[i].S;
        }
        double sum_nx = 0.0, sum_ny = 0.0;
        for (int i = 0; i < e.N_edges; i++)
        {
            int cellr;
            cellr = e.edges[i].celll;
            if (cellr == -1)
            {
                sum_nx += e.edges[i].n(0);
                sum_ny += e.edges[i].n(1);
            }
        }
        std::cout << std::scientific;
        std::cout << std::setprecision(10);
        std::cout << "Total volume of the domain = " << S << std::endl;
        std::cout << "Sum of boundary face normal (nx) = " << sum_nx << std::endl;
        std::cout << "Sum of boundary face normal (ny) = " << sum_ny << std::endl;
    }
    else
    {
        std::cout << "Failed reading grid file!" << std::endl;
        std::exit(-1);
    }
}
