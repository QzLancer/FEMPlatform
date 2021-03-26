#include "FEM2DNDDRSolver.h"

#include <chrono>
#include <omp.h>
void FEM2DNDDRSolver::solve()
{
	A.resize(m_num_nodes);
	At.resize(m_num_nodes);
	Bx.resize(m_num_triele);
	By.resize(m_num_triele);
	B.resize(m_num_triele);
	//计算三角形单元几何部分
	makeTrangle();
	//处理边界条件、材料和负载
	processBoundaryCondition();
	processMaterial();
	processLoad();
	//处理NDDR的节点信息
	nddrnode.resize(m_num_nodes);
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		for (int j = 0; j < 3; ++j) {
			int n = mp_triele[i_tri].n[j];
			nddrnode[n].NumberofNeighbourElement++;
			nddrnode[n].NeighbourElementId.push_back(i_tri);
			nddrnode[n].NeighbourElementNumber.push_back(j);
		}
	}
	////测试
	//for (int i = 0; i < m_num_nodes; ++i) {
	//	cout << "node: " << i << ", Number of Neighbour Element: " << nddrnode[i].NumberofNeighbourElement << endl;
	//	for (int j = 0; j < nddrnode[i].NumberofNeighbourElement; ++j) {
	//		cout << "NeighbourElementId: " << nddrnode[i].NeighbourElementId[j] << endl;
	//		cout << "NeighbourElementNumber: " << nddrnode[i].NeighbourElementNumber[j] << endl;
	//	}
	//}

	//NDDR迭代过程
	//周围所有节点都当成第一类边界条件处理
	//先把子域当成线性处理，能够收敛，但是求解效率很低
	//maxitersteps = 20000;
	//vector<double> At_old(m_num_nodes, 0);
	//for (int iter = 0; iter < maxitersteps; ++iter) {
	//	cout << "Iteration step " << iter + 1 << " start." << endl;
	//	for (int n = 0; n < m_num_nodes; ++n) { 
	//		if (mp_node[n].bdr == 1) {
	//			continue;
	//		}
	//		double S = 0, F = 0;
	//		for (int k = 0; k < nddrnode[n].NumberofNeighbourElement; ++k) {
	//			int i_tri = nddrnode[n].NeighbourElementId[k];
	//			CTriElement triele = mp_triele[i_tri];
	//			int nodenumber = nddrnode[n].NeighbourElementNumber[k];
	//			double mut = triele.material->getMu(B[i_tri]) * triele.xdot;
				//for (int i = 0; i < 3; ++i) {
				//	double Se = triele.C[nodenumber][i] / mut;
				//	if (nodenumber == i) {
				//		S += Se;
				//		F += triele.J * triele.area / 3;
				//	}
				//	else {
				//		F -= Se * At_old[triele.n[i]];
				//	}
				//}
	//		}
	//		At[n] = F / S;
	//		A[n] = At[n] / mp_node[n].x;
	//	}
	//	updateB();
		////判断收敛性
		//double error = 0, a = 0, b = 0;
		//for (int i = 0; i < m_num_nodes; ++i) {
		//	a += (At[i] - At_old[i]) * (At[i] - At_old[i]);
		//	b += At[i] * At[i];
		//}
		//error = sqrt(a) / sqrt(b);
		//cout << "Relative error: " << error << endl;
		//if (error > maxerror) {
		//	At_old = At;
		//}
		//else {
		//	cout << "Nonlinear NDDR iteration finish.\n";
		//	return;
		//}
	//}

	//子域内部采用牛顿迭代
	vector<double> At_old(m_num_nodes, 0);
	for (int iter = 0; iter < maxitersteps; ++iter) {
		//cout << "Iteration step " << iter + 1 << " start." << endl;
#pragma omp parallel for num_threads(8)
		for (int n = 0; n < m_num_nodes; ++n) {
			if (mp_node[n].bdr == 1) {
				continue;
			}
			//节点内部迭代过程
			int maxNRitersteps = 100;
			double Ati = 0;
			for (int NRiter = 0; NRiter < maxNRitersteps; ++NRiter) {
				double S = 0, F = 0;
				double J = 0, Fj = 0;
				//装配过程
				for (int k = 0; k < nddrnode[n].NumberofNeighbourElement; ++k) {
					//vector<vector<double>> J(3, vector<double>(3, 0));	//单元Jacobi矩阵
					//vector<double> Fj(3, 0);	//新增的右侧项
					int i_tri = nddrnode[n].NeighbourElementId[k];
					CTriElement triele = mp_triele[i_tri];
					int nodenumber = nddrnode[n].NeighbourElementNumber[k];
					double mut = triele.material->getMu(B[i_tri]) * triele.xdot;
					//处理线性单元
					if (triele.material->getLinearFlag() == true) {
						for (int i = 0; i < 3; ++i) {
							double Se = triele.C[nodenumber][i] / mut;
							if (nodenumber == i) {
								S += Se;
								F += triele.J * triele.area / 3;
							}
							else {
								F -= Se * At_old[triele.n[i]];
							}
						}
					}
					//处理非线性单元
					else {
						double mu, mut, dvdb, dvdbt, Bt, sigmai = 0, sigmaj = 0;
						mu = triele.material->getMu(B[i_tri]);
						mut = mu * triele.xdot;
						dvdb = triele.material->getdvdB(B[i_tri]);
						dvdbt = dvdb / triele.xdot / triele.xdot;
						Bt = B[i_tri] * triele.xdot;
						for (int i = 0; i < 3; ++i) {
							//sigmai和sigmaj中的At如何处理？
							for (int m = 0; m < 3; ++m) {
								if (m == nodenumber) {
									sigmai += triele.C[i][m] * Ati;
								}
								else {
									sigmai += triele.C[i][m] * At_old[triele.n[m]];
								}
							}
							for (int j = 0; j < 3; ++j) {
								for (int m = 0; m < 3; ++m) {
									if (m == nodenumber) {
										sigmaj += triele.C[j][m] * Ati;
									}
									else {
										sigmaj += triele.C[j][m] * At_old[triele.n[m]];
									}
								}
							}
						}
						for (int i = 0; i < 3; ++i) {
							if (Bt != 0) {
								J = triele.C[nodenumber][i] / mut + sigmai * sigmaj / Bt / triele.area;
							}
							else {
								J = triele.C[nodenumber][i] / mut;
							}
							if (nodenumber == i) {
								S += J;
								F += (J - triele.C[nodenumber][i] / mut) * Ati;
							}
							else {
								F += (J - triele.C[nodenumber][i] / mut) * At_old[triele.n[i]];
								F -= J * At_old[triele.n[i]];
							}
						}
					}
					
				}
				Ati = F / S;
				//NR迭代收敛性判断 
				double a = (Ati - At[n]) * (Ati - At[n]);
				double b = Ati * Ati;
				double NRerror = sqrt(a) / sqrt(b);
				if (Ati == 0) {
					break;
				}
				if (NRerror > maxerror) {
					At[n] = Ati;
					A[n] = At[n] / mp_node[n].x;
					for (int i = 0; i < nddrnode[n].NumberofNeighbourElement; ++i) {
						updateB(nddrnode[n].NeighbourElementId[i]);
					}
				}
				else {
					break;
				}
			}
		}
		//判断全局收敛性
		double error = 0, a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (At[i] - At_old[i]) * (At[i] - At_old[i]);
			b += At[i] * At[i];
		}
		error = sqrt(a) / sqrt(b);
		if ((iter + 1) % 100 == 0) {
			cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
		}
		if (error > maxerror) {
			At_old = At;
		}
		else {
			cout << "Iteration step: " << iter + 1 << endl;
			cout << "Nonlinear NDDR iteration finish.\n";
			return;
		}
	}
}
