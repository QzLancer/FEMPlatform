#include "FEM2DNRSolver.h"

#include <fstream>
void FEM2DNRSolver::solveStatic()
{
	//计算三角形单元几何部分
	makeTrangle();
	//处理边界条件、材料和负载
	processBoundaryCondition();
	processMaterial();
	processLoad();

	clock_t start, end;
	start = clock();
	if (dimension == FEMModel::DIMENSION::D2AXISM) {
		solve2DAxim();
	}
	else if (dimension == FEMModel::DIMENSION::D2PLANE) {
		solve2DPlane();
	}
	end = clock();
	cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
}

void FEM2DNRSolver::solve2DAxim()
{
	//考虑第一类边界条件的装配
	std::vector<std::vector<int>> locs(2, std::vector<int>(9 * (size_t)(m_num_triele)));
	std::vector<double> vals(9 * (size_t)(m_num_triele));
	std::vector<double> F(num_freenodes);
	int pos = 0;
	int linearelesize = 0, nonlinearelesize = 0;

	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		CTriElement triele = mp_triele[i_tri];
		if (triele.material->getLinearFlag() == true) linearelesize++;
		else nonlinearelesize++;
		for (int i = 0; i < 3; ++i) {
			int n1 = triele.n[i];
			for (int j = 0; j < 3; ++j) {
				int n2 = triele.n[j];
				if (triele.material->getLinearFlag() == true) {
					if (mp_node[n1].bdr != 1 && mp_node[n2].bdr != 1) {
						double mu = triele.material->getMu();	//存在mu=0的情况
						double mut = mu * triele.xdot;
						double Se = triele.C[i][j] / mut;
						locs[0][pos] = node_pos[n1];
						locs[1][pos] = node_pos[n2];
						vals[pos] = Se;
						++pos;
					}
				}
			}
			if (mp_node[n1].bdr != 1) {
				double Fe = triele.J * triele.area / 3;
				F[node_pos[n1]] += Fe;
				//永磁计算
				double h_c = triele.material->getH_c();
				double theta_m = triele.material->getTheta_m();
				F[node_pos[n1]] += h_c / 2 * (triele.R[i] * cos(theta_m) - triele.Q[i] * sin(theta_m));
			}
		}
	}

	//for (int i = 0; i < num_freenodes; ++i) {
	//	if (node_pos[i] < num_freenodes && F[node_pos[i]] != 0) {
	//		cout << "n: " << i << ", F: " << F[node_pos[i]] << endl;
	//	}
	//}

	//非线性部分迭代
	int pos1 = pos;
	std::vector<double> F1 = F;
	std::vector<double> A_old(m_num_nodes, 0);
	for (int step = 0; step < maxitersteps; ++step) {
		cout << "Iteration step " << step + 1 << " start." << endl;
		pos = pos1;
		F = F1;
		for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
			CTriElement triele = mp_triele[i_tri];
			if (triele.material->getLinearFlag() == false) {
				//计算单元Jacobi矩阵
				vector<vector<double>> J(3, vector<double>(3, 0));	//单元Jacobi矩阵
				vector<double> Fj(3, 0);	//新增的右侧项
				double mu, mut, dvdb, dvdbt, Bt, sigmai, sigmaj;
				mu = triele.material->getMu(mp_triele[i_tri].B);
				//cout << "mu: " << mu << endl;
				mut = mu * triele.xdot;
				dvdb = triele.material->getdvdB(mp_triele[i_tri].B);
				dvdbt = dvdb / triele.xdot / triele.xdot;
				Bt = mp_triele[i_tri].B * triele.xdot;
				vector<int> n(3);
				n[0] = triele.n[0], n[1] = triele.n[1], n[2] = triele.n[2];
				for (int i = 0; i < 3; ++i) {
					sigmai = (triele.C[i][0] * mp_node[n[0]].At + triele.C[i][1] * mp_node[n[1]].At + triele.C[i][2] * mp_node[n[2]].At);
					for (int j = 0; j < 3; ++j) {
						sigmaj = (triele.C[j][0] * mp_node[n[0]].At + triele.C[j][1] * mp_node[n[1]].At + triele.C[j][2] * mp_node[n[2]].At);
						if (Bt != 0) {
							J[i][j] = triele.C[i][j] / mut + dvdbt * sigmai * sigmaj / Bt / triele.area;
						}
						else {
							J[i][j] = triele.C[i][j] / mut;
						}
						Fj[i] += (J[i][j] - triele.C[i][j] / mut) * mp_node[n[j]].At;
					}
				}
				//装配
				for (int i = 0; i < 3; ++i) {
					int n1 = triele.n[i];
					for (int j = 0; j < 3; ++j) {
						int n2 = triele.n[j];
						if (mp_node[n1].bdr != 1 && mp_node[n2].bdr != 1) {
							locs[0][pos] = node_pos[n1];
							locs[1][pos] = node_pos[n2];
							vals[pos] = J[i][j];
							++pos;
						}
					}
					if (mp_node[n1].bdr != 1) {
						F[node_pos[n1]] += Fj[i];
					}
				}
			}
		}
		//求解
		locs[0].resize(pos);
		locs[1].resize(pos);

		vector<double> res1 = matsolver->solveMatrix(locs, vals, F, pos, num_freenodes);
		for (int i = 0; i < num_freenodes; ++i) {
			int index = node_reorder[i];
			mp_node[index].At = res1[i];
			mp_node[index].A = mp_node[index].At / mp_node[index].x;
		}
		//更新磁场结果
		updateB();
		if (nonlinearelesize == 0) {
			cout << "Linear problem solved!\n";
			return;
		}
		//判断收敛性
		double error = 0, a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (mp_node[i].A - A_old[i]) * (mp_node[i].A - A_old[i]);
			b += mp_node[i].A * mp_node[i].A;
		}
		error = sqrt(a) / sqrt(b);
		cout << "Relative error: " << error << endl;
		if (error > maxerror) {
			for (int i = 0; i < m_num_nodes; ++i) {
				A_old[i] = mp_node[i].A;
			}
		}
		else {
			cout << "Nonlinear iteration finish.\n";
			return;
		}
	}

	cout << "Warning: Number of iterations out of limit.\n";
}

void FEM2DNRSolver::solve2DPlane()
{
	//考虑第一类边界条件的装配
	std::vector<std::vector<int>> locs(2, std::vector<int>(9 * (size_t)(m_num_triele)));
	std::vector<double> vals(9 * (size_t)(m_num_triele));
	std::vector<double> F(num_freenodes);
	int pos = 0;
	int linearelesize = 0, nonlinearelesize = 0;

	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		CTriElement triele = mp_triele[i_tri];
		if (triele.material->getLinearFlag() == true) linearelesize++;
		else nonlinearelesize++;
		for (int i = 0; i < 3; ++i) {
			int n1 = triele.n[i];
			for (int j = 0; j < 3; ++j) {
				int n2 = triele.n[j];
				if (triele.material->getLinearFlag() == true) {
					if (mp_node[n1].bdr != 1 && mp_node[n2].bdr != 1) {
						double mu = triele.material->getMu();	//存在mu=0的情况
						double Se = triele.C[i][j] / mu;
						locs[0][pos] = node_pos[n1];
						locs[1][pos] = node_pos[n2];
						vals[pos] = Se;
						++pos;
					}
				}
			}
			if (mp_node[n1].bdr != 1) {
				double Fe = triele.J * triele.area / 3;
				F[node_pos[n1]] += Fe;
				//永磁计算
				double h_c = triele.material->getH_c();
				double theta_m = triele.material->getTheta_m();
				F[node_pos[n1]] += h_c / 2 * (triele.R[i] * cos(theta_m) - triele.Q[i] * sin(theta_m));
			}
		}
	}

	//for (int i = 0; i < num_freenodes; ++i) {
	//	if (node_pos[i] < num_freenodes && F[node_pos[i]] != 0) {
	//		cout << "n: " << i << ", F: " << F[node_pos[i]] << endl;
	//	}
	//}

	//非线性部分迭代
	int pos1 = pos;
	std::vector<double> F1 = F;
	std::vector<double> A_old(m_num_nodes, 0);
	for (int step = 0; step < maxitersteps; ++step) {
		cout << "Iteration step " << step + 1 << " start." << endl;
		pos = pos1;
		F = F1;
		for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
			CTriElement triele = mp_triele[i_tri];
			if (triele.material->getLinearFlag() == false) {
				//计算单元Jacobi矩阵
				vector<vector<double>> J(3, vector<double>(3, 0));	//单元Jacobi矩阵
				vector<double> Fj(3, 0);	//新增的右侧项
				double mu, dvdb, B, sigmai, sigmaj;
				//mu = triele.material->getMu(B[i_tri]);
				B = mp_triele[i_tri].B;
				//cout << "i_tri: " << i_tri << ", B: " << B << endl;
				mu = triele.material->getMu(B);
				//cout << "mu: " << mu << endl;
				dvdb = triele.material->getdvdB(B);
				vector<int> n(3);
				n[0] = triele.n[0], n[1] = triele.n[1], n[2] = triele.n[2];
				for (int i = 0; i < 3; ++i) {
					sigmai = (triele.C[i][0] * mp_node[n[0]].A + triele.C[i][1] * mp_node[n[1]].A + triele.C[i][2] * mp_node[n[2]].A);
					for (int j = 0; j < 3; ++j) {
						sigmaj = (triele.C[j][0] * mp_node[n[0]].A + triele.C[j][1] * mp_node[n[1]].A + triele.C[j][2] * mp_node[n[2]].A);
						if (B != 0) {
							J[i][j] = triele.C[i][j] / mu + dvdb * sigmai * sigmaj / B / triele.area;
						}
						else {
							J[i][j] = triele.C[i][j] / mu;
						}
						Fj[i] += (J[i][j] - triele.C[i][j] / mu) * mp_node[n[j]].A;
					}
				}
				//装配
				for (int i = 0; i < 3; ++i) {
					int n1 = triele.n[i];
					for (int j = 0; j < 3; ++j) {
						int n2 = triele.n[j];
						if (mp_node[n1].bdr != 1 && mp_node[n2].bdr != 1) {
							locs[0][pos] = node_pos[n1];
							locs[1][pos] = node_pos[n2];
							vals[pos] = J[i][j];
							++pos;
						}
					}
					if (mp_node[n1].bdr != 1) {
						F[node_pos[n1]] += Fj[i];
					}
				}
			}
		}
		//求解
		locs[0].resize(pos);
		locs[1].resize(pos);

		vector<double> res1 = matsolver->solveMatrix(locs, vals, F, pos, num_freenodes);
		for (int i = 0; i < num_freenodes; ++i) {
			int index = node_reorder[i];
			mp_node[index].A = res1[i];
			//cout << "mp_node[index].A: " << mp_node[index].A << endl;
		}




		//更新磁场结果
		updateB();
		if (nonlinearelesize == 0) {
			cout << "Linear problem solved!\n";
			return;
		}
		//判断收敛性
		double error = 0, a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (mp_node[i].A - A_old[i]) * (mp_node[i].A - A_old[i]);
			b += mp_node[i].A * mp_node[i].A;
		}
		error = sqrt(a) / sqrt(b);
		cout << "Relative error: " << error << endl;
		if (error > maxerror) {
			for (int i = 0; i < m_num_nodes; ++i) {
				A_old[i] = mp_node[i].A;
			}
		}
		else {
			cout << "Nonlinear iteration finish.\n";
			return;
		}
	}

	cout << "Warning: Number of iterations out of limit.\n";
}


