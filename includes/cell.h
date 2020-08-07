#ifndef _cell_h
#define _cell_h

#include <vector>
#include <Eigen/Dense>

// Cell class

class cell
{

public:
  // int vertex1;
  // int vertex2;
  // int vertex3;
  std::vector<int> vertices;
  Eigen::Vector3d r;
  Eigen::Vector3d n;
  double S;
  int type;
  Eigen::VectorXd Q;
  Eigen::VectorXd R;
  double lambda_max;
  double dtau;
  cell(void);
  cell(std::vector<int>);
  cell(int, int, int);
  cell(int, int, int, int);
  ~cell(void);
};

#endif
