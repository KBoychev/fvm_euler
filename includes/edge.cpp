#include "edge.h"

edge::edge(int vertex1, int vertex2, int celll, int cellr)
{
    this->vertex1 = vertex1;
    this->vertex2 = vertex2;
    this->celll = celll;
    this->cellr = cellr;
    this->Ql = Eigen::VectorXd::Zero(4);
    this->Qr = Eigen::VectorXd::Zero(4);
    this->unique = 1;
}

edge::~edge(void) {}