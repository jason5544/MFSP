#include "StdAfx.h"
#include "process.h"


process::process(void)
{
}


process::~process(void)
{
}


//to obtain the length of shortest path between two diseases
MatrixXd process::get_descendant(const MatrixXd& ex)
{
	int dis_num = ex.rows();
	MatrixXd output=MatrixXd::Zero(dis_num,dis_num);
	for (int i=0; i!=dis_num; i++)
	{
		stack <int> st;
		st.push(i);
		int j =0;

		while (!st.empty())
		{
			while (j!=dis_num)
			{
				if (ex(st.top(),j)!=0)
				{
					st.push(j);
					int temp = st.size();
					if ( output(i,j) > temp ||output(i,j)==0 )
					{
						output(i,j) = temp;

					}
					j = 0;
				}
				else
				{
					j ++;
				}
			}
			j = st.top()+1;
			st.pop();
		}
		output(i,i)  =1;
	}
	return output;
}


//to calculate the similarity of disease
MatrixXd process::cal_dis_sim (const MatrixXd& ex)
{
	int dis_num = ex.cols();
	MatrixXd temp = MatrixXd::Zero(dis_num,dis_num);
	MatrixXd output = MatrixXd::Zero(dis_num,dis_num);
	for (int i=0; i!=dis_num; i++)
	{
		for (int j=0; j!=dis_num; j++)
		{
			if (ex(i,j)!=0)
			{
				temp(i,j) = 1.0 / ex(i,j);
			}	
		}
	}

	for (int i=0; i!=dis_num; i++)
	{
		if (temp.col(i).sum()!=0)
		{
			for (int j=i; j!=dis_num; j++)
			{
				if (temp.col(j).sum()!=0)
				{
					double fenzi = temp.col(i).transpose()*temp.col(j);
					double fenmu1 = sqrt(static_cast <double> (temp.col(i).transpose()*temp.col(i)));
					double fenmu2 = sqrt(static_cast <double> (temp.col(j).transpose()*temp.col(j)));
					output(i,j) = fenzi / (fenmu1*fenmu2);
				}
				else
				{
					cout << " error" << endl;

				}
			}
		}
		else
		{
			cout << " error" << endl;
		}
	}

	for (int i=0; i!=dis_num; i++)
	{
		for (int j=0; j!=i; j++)
		{
			output(i,j) = output(j,i);
		}
	}
	return output;
}


//to calculate the similarity of miRNA
MatrixXd process::cal_Sim(const MatrixXd & MD, const MatrixXd& DD, const int& iter_num, const double& b)
{
	int miRNA_num = MD.rows();
	MatrixXd out = MatrixXd::Zero(miRNA_num,miRNA_num);
	double out_sum = 0;
	for (int i=0; i<=iter_num; i++)
	{
		MatrixXd out_temp = merge_path(MD,DD,i);
		double d_temp = pow(b,i);
		out += d_temp * out_temp;
		out_sum += d_temp;
	}
	out = out * (1/out_sum);

	MatrixXd sim = MatrixXd::Zero(miRNA_num,miRNA_num);
	for (int i=0; i!=miRNA_num; i++)
	{
		for (int j=0; j!=miRNA_num; j++)
		{
			sim(i,j) = 2*out(i,j)/(out(i,i)+out(j,j));
		}
	}
	return sim;
}


MatrixXd process::merge_path (const MatrixXd & MD, const MatrixXd& DD, const int& flag)
{
	MatrixXd out;
	switch (flag)
	{
	case 0: out = MD * MD.transpose();
		break;
	case 1: out = MD * DD * MD.transpose();
		break;
	case 2: out = MD * DD * DD * MD.transpose();
		break;
	case 3: out = MD * DD * DD * DD * MD.transpose();
		break;
	case 4: out = MD * DD * DD * DD * DD * MD.transpose();
		break;
	case 5: out = MD * DD * DD * DD * DD * DD * MD.transpose();
		break;
	case 6: out = MD * DD * DD * DD * DD * DD * DD * MD.transpose();
		break;
	case 7: out = MD * DD * DD * DD * DD * DD * DD * DD * MD.transpose();
		break;
	case 8: out = MD * DD * DD * DD * DD * DD * DD * DD * DD * MD.transpose();
		break;
	case 9:out = MD * DD * DD * DD * DD * DD * DD * DD * DD * DD * MD.transpose();
		break;
	case 10:out = MD * DD * DD * DD * DD * DD * DD * DD * DD * DD * DD * MD.transpose();
		break;
	default: cout << "error" << endl;
	}


	return out;
}
