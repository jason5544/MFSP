#pragma once

#include <stack>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <math.h>
using namespace Eigen;  
using namespace Eigen::internal;  
using namespace Eigen::Architecture;
using namespace std;


class process
{
public:
	static MatrixXd get_descendant(const MatrixXd& ex);

	static MatrixXd cal_dis_sim (const MatrixXd& ex);

	static MatrixXd cal_Sim (const MatrixXd & MD, const MatrixXd& DD, const int& iter_num, const double& b);
	
	static MatrixXd  merge_path(const MatrixXd & MD, const MatrixXd& DD, const int& flag);

	process(void);
	~process(void);
};

