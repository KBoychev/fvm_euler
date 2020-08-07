#ifndef _edge_h
#define _edge_h

#include <Eigen/Dense>

// Edge class

class edge
{
public:
    int vertex1;
    int vertex2;
    int celll;
    int cellr;
    double l;
    int unique;
    Eigen::Vector3d n;
    Eigen::Vector3d r;
    Eigen::VectorXd Ql;
    Eigen::VectorXd Qr;
    edge(int, int, int, int);
    ~edge(void);
};
#endif
