// We implemented MFSP in C++ with visual studio 2010 based on Windows7,
// and used the eigen3 (downloaded from http://eigen.tuxfamily.org/)
// submitted by Pingjian Ding
// email:dpjhnu@qq.com

//#include "stdafx.h"
#include "process.h"

double a = 0.6; //weight ratio
int b = 2; //maximum transferring times

int _tmain(int argc, _TCHAR* argv[])
{
	//if the node i is the parent of node j, the element (i,j)-th of matrix DAG is 1, 0 otherwise. 
	MatrixXd DAG(6,6);
	DAG << 0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0;

	MatrixXd RMD(2,6);
	RMD << 1,1,1,1,1,0,1,1,0,0,1,1;

	MatrixXd decens = process::get_descendant(DAG);	
	MatrixXd RDD = process::cal_dis_sim(decens);
	cout << "the similarity matrix of disease:" <<endl << RDD << endl << endl;

	MatrixXd RMM = process::cal_Sim(RMD,RDD,b,a);
	cout << "the similarity matrix of miRNA:" <<endl << RMM << endl;

	return 0;
}

