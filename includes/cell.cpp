
#include "cell.h"

cell::cell(void)
{
}
cell::cell(std::vector<int> vertices)
{
    this->vertices = vertices;

    this->Q = Eigen::VectorXd::Zero(4);
    this->R = Eigen::VectorXd::Zero(4);
    this->type = 0;
    this->lambda_max = 0;
    this->dtau = 0;
}

cell::cell(int vertex1, int vertex2, int vertex3)
{
    this->vertices.push_back(vertex1);
    this->vertices.push_back(vertex2);
    this->vertices.push_back(vertex3);

    this->Q = Eigen::VectorXd::Zero(4);
    this->R = Eigen::VectorXd::Zero(4);
    this->type = 0;
    this->lambda_max = 0;
    this->dtau = 0;
}

cell::cell(int vertex1, int vertex2, int vertex3, int vertex4)
{
    this->vertices.push_back(vertex1);
    this->vertices.push_back(vertex2);
    this->vertices.push_back(vertex3);
    this->vertices.push_back(vertex4);

    this->Q = Eigen::VectorXd::Zero(4);
    this->R = Eigen::VectorXd::Zero(4);
    this->type = 0;
    this->lambda_max = 0;
    this->dtau = 0;
}

cell::~cell(void) {}