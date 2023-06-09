#include "FEM2DNDDRSolver.h"

#include <chrono>
#include <omp.h>
void FEM2DNDDRSolver::solveStatic()
{
	//A.resize(m_num_nodes);
	//At.resize(m_num_nodes);
	//Bx.resize(m_num_triele);
	//By.resize(m_num_triele);
	//B.resize(m_num_triele);
	//计算三角形单元几何部分
	makeTrangle();
	//处理边界条件、材料和负载
	processBoundaryCondition();
	processMaterial();
	processLoad();
	//处理NDDR的节点信息
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		for (int j = 0; j < 3; ++j) {
			int n = mp_triele[i_tri].n[j];
			int id = mp_node[n].NumberofNeighbourElement;
			mp_node[n].NeighbourElementId[id] = i_tri;
			mp_node[n].NeighbourElementNumber[id] = j;
			mp_node[n].NumberofNeighbourElement++;
		}
	}
	//测试
	//for (int i = 0; i < m_num_nodes; ++i) {
	//	cout << "node: " << i << ", Number of Neighbour Element: " << nddrnode[i].NumberofNeighbourElement << endl;
	//	for (int j = 0; j < nddrnode[i].NumberofNeighbourElement; ++j) {
	//		cout << "NeighbourElementId: " << nddrnode[i].NeighbourElementId[j] << endl;
	//		cout << "NeighbourElementNumber: " << nddrnode[i].NeighbourElementNumber[j] << endl;
	//	}
	//	cout << endl;
	//}

	//NDDR迭代过程
	//周围所有节点都当成第一类边界条件处理
	//先把子域当成线性处理，能够收敛，但是求解效率很低
//	maxitersteps = 20000;
//	vector<double> At_old(m_num_nodes, 0);
//	for (int iter = 0; iter < maxitersteps; ++iter) {
//		cout << "Iteration step " << iter + 1 << " start." << endl << endl;
//#pragma omp parallel for num_threads(8)
//		for (int n = 0; n < m_num_nodes; ++n) { 
//			if (mp_node[n].bdr == 1) {
//				continue;
//			}
//			double S = 0, F = 0;
//			for (int k = 0; k < nddrnode[n].NumberofNeighbourElement; ++k) {
//				int i_tri = nddrnode[n].NeighbourElementId[k];
//				CTriElement triele = mp_triele[i_tri];
//				int nodenumber = nddrnode[n].NeighbourElementNumber[k];
//				double mut = triele.material->getMu(B[i_tri]) * triele.xdot;
//				for (int i = 0; i < 3; ++i) {
//					double Se = triele.C[nodenumber][i] / mut;
//					//cout << "Se: " << Se << endl;
//					if (nodenumber == i) {
//						//if (Se == 0) {
//						//	cout << "Se < 0, node: " << n << ", i:" << i << ", i_tri: " << i_tri << endl;
//						//}
//						S += Se;
//						F += triele.J * triele.area / 3;
//						//永磁计算
//						double h_c = triele.material->getH_c();
//						double theta_m = triele.material->getTheta_m();
//						F += h_c / 2 * (triele.R[i] * cos(theta_m) - triele.Q[i] * sin(theta_m));
//					}
//					else {
//						//直角或钝角三角形的情况
//						if (Se >= 0) {
//							//cout << "i_tri: " << i_tri << ", mut: " << mut << ", domain: " << triele.domain << endl;
//							////cout << "Se: " << endl;
//							////for (int ii = 0; ii < 3; ++ii) {
//							////	for (int jj = 0; jj < 3; ++jj) {
//							////		cout << triele.C[ii][jj] << " ";
//							////	}
//							////	cout << endl;
//							////}
//							////cout << "Coor: " << endl;
//							////for (int ii = 0; ii < 3; ++ii) {
//							////	int n1 = triele.n[ii];
//							////	printf("x[%d]: %.6f, y[%d]: %.6f, At_old[%d]: %.6f\n", ii, mp_node[n1].x, ii, mp_node[n1].y, ii, At_old[n1]);
//							////}
//							//cout << "Se * At_old[triele.n[i]]: " << Se * At_old[triele.n[i]] << endl;
//							//cout << endl;
//							//Se = 0;
//						}
//						F -= Se * At_old[triele.n[i]];
//					}
//				}
//			}
//			At[n] = F / S;
//			A[n] = At[n] / mp_node[n].x;
//		}
//		updateB();
//
//		//判断收敛性
//		double error = 0, a = 0, b = 0;
//		for (int i = 0; i < m_num_nodes; ++i) {
//			a += (At[i] - At_old[i]) * (At[i] - At_old[i]);
//			b += At[i] * At[i];
//		}
//		error = sqrt(a) / sqrt(b);
//		cout << "Relative error: " << error << endl;
// 		if (error > maxerror) {
//			At_old = At;
//		}
//		else {
//			cout << "Nonlinear NDDR iteration finish.\n";
//			return;
//		}
//	}

	//子域内部采用牛顿迭代
	clock_t start, end;
	start = clock();
	if (dimension == FEMModel::DIMENSION::D2AXISM) {
		//solve2DAxim();
		//solve2DAxim1();
		solve2DAxim2();
	}
	else if (dimension == FEMModel::DIMENSION::D2PLANE) {
		//solve2DPlane();
		solve2DPlane1();
		//solve2DPlane2();
	}
	end = clock();
	cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
}

void FEM2DNDDRSolver::solve2DAxim()
{
//	//vector<double> At_old(m_num_nodes, 0);
//	for (int iter = 0; iter < maxitersteps; ++iter) {
//		//cout << "Iteration step " << iter + 1 << " start." << endl;
//#pragma omp parallel for num_threads(8)
//		for (int n = 0; n < m_num_nodes; ++n) {
//			if (mp_node[n].bdr == 1) {
//				continue;
//			}
//			//节点内部迭代过程
//			double A, dA, Jac, ku, J0, k0, NRerror, Js;
//			double A1, A2, A3;
//			double B2, B, V, VB2, B2A;
//			int I, J, K, ID;
//			int i, RelaxCount, count = 0;
//			double maxNRitersteps = 6;
//			double h_c, theta_m;
//			for (int NRiter = 0; NRiter < maxNRitersteps; ++NRiter) {
//				Jac = 0, ku = 0;
//				//装配过程
//				for (int k = 0; k < mp_node[n].NumberofNeighbourElement; ++k) {
//					//vector<vector<double>> J(3, vector<double>(3, 0));	//单元Jacobi矩阵
//					//vector<double> Fj(3, 0);	//新增的右侧项
//					int i_tri = mp_node[n].NeighbourElementId[k];
//					CTriElement triele = mp_triele[i_tri];
//					int nodenumber = mp_node[n].NeighbourElementNumber[k];
//					double mut = triele.material->getMu(triele.B) * triele.xdot;
//					//printf("mu: %f\n", mut);
//					//处理线性单元
//					if (triele.material->getLinearFlag() == true) {
//						for (int i = 0; i < 3; ++i) {
//							double Se = triele.C[nodenumber][i] / mut;
//							if (nodenumber == i) {
//								J0 = Se;
//								k0 = triele.J * triele.area / 3;
//								//永磁部分
//								h_c = triele.material->getH_c();
//								theta_m = triele.material->getTheta_m();
//								k0 += h_c / 2 * (triele.R[i] * cos(theta_m) - triele.Q[i] * sin(theta_m));
//							}
//							else {
//								k0 -= Se * mp_node[triele.n[i]].At_old;
//							}
//						}
//					}
//					//处理非线性单元
//					else {
//						double mu, mut, dvdb, dvdbt, Bt, sigmai[3]{ 0, 0, 0 };
//						mu = triele.material->getMu(triele.B);
//						mut = mu * triele.xdot;
//						dvdb = triele.material->getdvdB(triele.B);
//						dvdbt = dvdb / triele.xdot / triele.xdot;
//						Bt = triele.B * triele.xdot;
//						for (int i = 0; i < 3; ++i) {
//							for (int m = 0; m < 3; ++m) {
//								if (m == nodenumber) {
//									sigmai[i] += triele.C[i][m] * Ati;
//								}
//								else {
//									sigmai[i] += triele.C[i][m] * mp_node[triele.n[m]].At_old;
//								}
//							}
//						}
//						for (int i = 0; i < 3; ++i) {
//							if (Bt != 0) {
//								J = triele.C[nodenumber][i] / mut + dvdbt * sigmai[nodenumber] * sigmai[i] / Bt / triele.area;
//							}
//							else {
//								J = triele.C[nodenumber][i] / mut;
//							}
//							if (nodenumber == i) {
//								S += J;
//								F += (J - triele.C[nodenumber][i] / mut) * Ati;
//							}
//							else {
//								F += (J - triele.C[nodenumber][i] / mut) * mp_node[triele.n[i]].At_old;
//								F -= J * mp_node[triele.n[i]].At_old;
//							}
//						}
//						//if (sigmai != 0 && sigmaj != 0) 
//						//	cout << "sigmai: " << sigmai << ", sigmaj: " << sigmaj << endl;
//					}
//
//				}
//				Ati = F / S;
//				//NR迭代收敛性判断 
//				double a = (Ati - mp_node[n].At) * (Ati - mp_node[n].At);
//				double b = Ati * Ati;
//				double NRerror = sqrt(a) / sqrt(b);
//				if (Ati == 0) {
//					continue;
//				}
//				if (NRerror > maxerror) {
//					mp_node[n].At = Ati;
//					mp_node[n].A = mp_node[n].At / mp_node[n].x;
//					for (int i = 0; i < mp_node[n].NumberofNeighbourElement; ++i) {
//						updateB(mp_node[n].NeighbourElementId[i]);
//					}
//				}
//				else {
//					//printf("NRerror: %f\n", NRerror);
//					//cout << "NRerror: " << NRerror << endl;
//					//printf("n: %d, NRiter: %d\n", n, NRiter);
//					break;
//				}
//			}
//		}
//
//		//判断全局收敛性
//		//double error = 0, a = 0, b = 0;
//		//for (int i = 0; i < m_num_nodes; ++i) {
//		//	a += (mp_node[i].At - mp_node[i].At_old) * (mp_node[i].At - mp_node[i].At_old);
//		//	b += mp_node[i].At * mp_node[i].At;
//		//}
//		//error = sqrt(a) / sqrt(b);
//		//cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
//
//		//if ((iter + 1) % 100 == 0) {
//		//	cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
//		//}
//		//if (error > maxerror) {
//		//	for (int i = 0; i < m_num_nodes; ++i) {
//		//		At_old[i] = mp_node[i].At;
//		//	}
//		//}
//		//else {
//		//	cout << "Iteration step: " << iter + 1 << endl;
//		//	cout << "Nonlinear NDDR iteration finish.\n";
//		//	return;
//		//}
//#pragma omp parallel for num_threads(8)
//		for (int i = 0; i < m_num_nodes; ++i) {
//			mp_node[i].At_old = mp_node[i].At;
//		}
//
//	}
}

void FEM2DNDDRSolver::solve2DAxim1()
{
	//vector<double> At_old(m_num_nodes, 0);
	for (int iter = 0; iter < maxitersteps; ++iter) {
		//cout << "Iteration step " << iter + 1 << " start." << endl;
#pragma omp parallel for num_threads(8)
		for (int n = 0; n < m_num_nodes; ++n) {
			if (mp_node[n].bdr == 1) {
				continue;
			}
			//节点内部迭代过程
			int maxNRitersteps = 6;
			double Ati = 0;
			for (int NRiter = 0; NRiter < maxNRitersteps; ++NRiter) {
				double S = 0, F = 0;
				double J = 0, Fj = 0;
				//装配过程
				for (int k = 0; k < mp_node[n].NumberofNeighbourElement; ++k) {
					//vector<vector<double>> J(3, vector<double>(3, 0));	//单元Jacobi矩阵
					//vector<double> Fj(3, 0);	//新增的右侧项
					int i_tri = mp_node[n].NeighbourElementId[k];
					CTriElement triele = mp_triele[i_tri];
					int nodenumber = mp_node[n].NeighbourElementNumber[k];
					double mut = triele.material->getMu(triele.B) * triele.xdot;
					//printf("mu: %f\n", mut);
					//处理线性单元
					if (triele.material->getLinearFlag() == true) {
						for (int i = 0; i < 3; ++i) {
							double Se = triele.C[nodenumber][i] / mut;
							if (nodenumber == i) {
								S += Se;
								F += triele.J * triele.area / 3;
								//永磁部分
								double h_c = triele.material->getH_c();
								double theta_m = triele.material->getTheta_m();
								F += h_c / 2 * (triele.R[i] * cos(theta_m) - triele.Q[i] * sin(theta_m));
							}
							else {
								F -= Se * mp_node[triele.n[i]].At_old;
							}
						}
					}
					//处理非线性单元
					else {
						double mu, mut, dvdb, dvdbt, Bt, sigmai[3]{0, 0, 0};
						mu = triele.material->getMu(triele.B);
						mut = mu * triele.xdot;
						dvdb = triele.material->getdvdB(triele.B);
						dvdbt = dvdb / triele.xdot / triele.xdot;
						Bt = triele.B * triele.xdot;
						for (int i = 0; i < 3; ++i) {
							for (int m = 0; m < 3; ++m) {
								if (m == nodenumber) {
									sigmai[i] += triele.C[i][m] * Ati;
								}
								else {
									sigmai[i] += triele.C[i][m] * mp_node[triele.n[m]].At_old;
								}
							}
						}
						for (int i = 0; i < 3; ++i) {
							if (Bt != 0) {
								J = triele.C[nodenumber][i] / mut + dvdbt * sigmai[nodenumber] * sigmai[i] / Bt / triele.area;
							}
							else {
								J = triele.C[nodenumber][i] / mut;
							}
							if (nodenumber == i) {
								S += J;
								F += (J - triele.C[nodenumber][i] / mut) * Ati;
							}
							else {
								F += (J - triele.C[nodenumber][i] / mut) * mp_node[triele.n[i]].At_old;
								F -= J * mp_node[triele.n[i]].At_old;
							}
						}
						//if (sigmai != 0 && sigmaj != 0) 
						//	cout << "sigmai: " << sigmai << ", sigmaj: " << sigmaj << endl;
					}

				}
				Ati = F / S;
				//NR迭代收敛性判断 
				double a = (Ati - mp_node[n].At) * (Ati - mp_node[n].At);
				double b = Ati * Ati;
				double NRerror = sqrt(a) / sqrt(b);
				if (Ati == 0) {
					continue;
				}
				if (NRerror > maxerror) {
					mp_node[n].At = Ati;
					mp_node[n].A = mp_node[n].At / mp_node[n].x;
					for (int i = 0; i < mp_node[n].NumberofNeighbourElement; ++i) {
						updateB(mp_node[n].NeighbourElementId[i]);
					}
				}
				else {
					//printf("NRerror: %f\n", NRerror);
					//cout << "NRerror: " << NRerror << endl;
					//printf("n: %d, NRiter: %d\n", n, NRiter);
					break;
				}
			}
		}

		//判断全局收敛性
		double error = 0, a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (mp_node[i].At - mp_node[i].At_old) * (mp_node[i].At - mp_node[i].At_old);
			b += mp_node[i].At * mp_node[i].At;
		}
		error = sqrt(a) / sqrt(b);
		cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;

		if ((iter + 1) % 100 == 0) {
			cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
		}
		if (error > maxerror) {
			for (int i = 0; i < m_num_nodes; ++i) {
				mp_node[i].At_old = mp_node[i].At;
			}
		}
		else {
			cout << "Iteration step: " << iter + 1 << endl;
			cout << "Nonlinear NDDR iteration finish.\n";
			return;
		}
//#pragma omp parallel for num_threads(8)
//		for (int i = 0; i < m_num_nodes; ++i) {
//			mp_node[i].At_old = mp_node[i].At;
//		}

	}
}

void FEM2DNDDRSolver::solve2DAxim2()
{
	for (int iter = 0; iter < maxitersteps; ++iter) {
#pragma omp parallel for num_threads(8)
		for (int i = 0; i < m_num_nodes; ++i) {
			double A, At, dAt, Jac, k, J0, k0, NRerror, Js;
			double B2, B, V, Vt, VB2, B2A;
			int ID;
			int RelaxCount, count = 0;
			int NRCount = 10;
			double RelaxFactor = 1;
			double AtLocal[3]{ 0 ,0, 0 }, ALocal[3]{0, 0, 0};
			if (mp_node[i].bdr != 1)
			{
				At = 0; dAt = 0; NRerror = 1; count = 0;
				while (count < NRCount)
				{
					Jac = 0; k = 0;
					for (int j = 0; j < mp_node[i].NumberofNeighbourElement; j++)
					{
						ID = mp_node[i].NeighbourElementId[j];
						Js = mp_triele[ID].J * mp_triele[ID].area / 3;
						int nodenumber = mp_node[i].NeighbourElementNumber[j];
						for (int m = 0; m < 3; ++m) {
							if (m == nodenumber) {
								AtLocal[m] = At;
							}
							else {
								AtLocal[m] = mp_node[mp_triele[ID].n[m]].At_old;
							}
							ALocal[m] = AtLocal[m] / mp_triele[ID].xdot;
							A = At / mp_triele[ID].xdot;
						}
						double RHSContri = 0;
						for (int m = 0; m < 3; ++m) {

							RHSContri += mp_triele[ID].C[nodenumber][m] * AtLocal[m];
						}

						//线性单元，计算右侧列向量时，也要考虑A？？？
						if (mp_triele[ID].material->getLinearFlag() == true)	//线性空气单元
						{
							Vt = 1.0 / mp_triele[ID].material->getMu(0) / mp_triele[ID].xdot;
							J0 = Vt * mp_triele[ID].C[nodenumber][nodenumber];
							k0 = Js - Vt * RHSContri;

							//V = 1.0 / mp_triele[ID].material->getMu(0);
							//J0 = V * mp_triele[ID].C[nodenumber][nodenumber];
							//k0 = Js - V * RHSContri;
						}

						else	//非线性单元
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (ALocal[0] - ALocal[1]) * (ALocal[0] - ALocal[1]) + mp_triele[ID].C[1][2] * (ALocal[1] - ALocal[2]) * (ALocal[1] - ALocal[2]) + mp_triele[ID].C[0][2] * (ALocal[0] - ALocal[2]) * (ALocal[0] - ALocal[2]));
							B = sqrt(B2);	//计算B，用于处理非线性
							Vt = 1.0 / mp_triele[ID].material->getMu(B) / mp_triele[ID].xdot;
							VB2 = mp_triele[ID].material->getdvdB2(B);
							B2A = -1 / mp_triele[ID].area * (2 * mp_triele[ID].C[nodenumber][0] * (A - ALocal[0]) + 2 * mp_triele[ID].C[nodenumber][1] * (A - ALocal[1]) + 2 * mp_triele[ID].C[nodenumber][2] * (A - ALocal[2]));
							J0 = B2A * VB2 * RHSContri / mp_triele[ID].xdot / mp_triele[ID].xdot + Vt * mp_triele[ID].C[nodenumber][nodenumber];
							k0 = Js - Vt * RHSContri;

							//B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (AtLocal[0] - AtLocal[1]) * (AtLocal[0] - AtLocal[1]) + mp_triele[ID].C[1][2] * (AtLocal[1] - AtLocal[2]) * (AtLocal[1] - AtLocal[2]) + mp_triele[ID].C[0][2] * (AtLocal[0] - AtLocal[2]) * (AtLocal[0] - AtLocal[2]));
							//B = sqrt(B2);	//计算B，用于处理非线性
							//V = 1.0 / mp_triele[ID].material->getMu(B);
							//VB2 = mp_triele[ID].material->getdvdB2(B);
							//B2A = -1 / mp_triele[ID].area * (2 * mp_triele[ID].C[nodenumber][0] * (A - AtLocal[0]) + 2 * mp_triele[ID].C[nodenumber][1] * (A - AtLocal[1]) + 2 * mp_triele[ID].C[nodenumber][2] * (A - AtLocal[2]));
							//J0 = B2A * VB2 * RHSContri + V * mp_triele[ID].C[nodenumber][nodenumber];
							//k0 = Js - V * RHSContri;
						}

						Jac = Jac + J0;
						k = k + k0;
					}
					dAt = k / Jac;
					At = At + RelaxFactor * dAt;
					//printf("k: %f, Jac: %f, dAt: %.20f, At: %.20f\n", k, Jac, dAt, At);
					double a = dAt * dAt;
					double b = At * At;
					double NRerror = sqrt(a) / sqrt(b);
					if (At == 0) {
						break;
					}
					if (NRerror > maxerror) {
						continue;
					}
					else {
						//if (NRiter != 1) {
						//	printf("NRerror: %.20f\n", NRerror);
						//	cout << "NRerror: " << NRerror << endl;
						//	printf("n: %d, NRiter: %d\n", n, NRiter);
						//}
						break;
					}
					count++;
				}
				mp_node[i].At = At;
			}
		}

		//判断全局收敛性
		double error = 0, a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (mp_node[i].At - mp_node[i].At_old) * (mp_node[i].At - mp_node[i].At_old);
			b += mp_node[i].At * mp_node[i].At;
		}
		error = sqrt(a) / sqrt(b);
		if ((iter + 1) % 100 == 0) {
			cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
		}
		cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
		if (error > maxerror) {
			//这一步影响效率
			for (int i = 0; i < m_num_nodes; ++i) {
				mp_node[i].At_old = mp_node[i].At;
			}
		}
		else {
			cout << "Iteration step: " << iter + 1 << endl;
			cout << "Nonlinear NDDR iteration finish.\n";
			return;
		}

		for (int i = 0; i < m_num_nodes; ++i) {
			if (mp_node[i].x != 0)
				mp_node[i].A = mp_node[i].At / mp_node[i].x;
		}
//#pragma omp parallel for num_threads(8)
//		for (int i = 0; i < m_num_nodes; ++i) {
//			mp_node[i].At_old = mp_node[i].At;
//		}
		updateB();
	}
}

void FEM2DNDDRSolver::solve2DPlane()
{
	vector<double> A_old(m_num_nodes, 0);
	for (int iter = 0; iter < maxitersteps; ++iter) {
		//cout << "Iteration step " << iter + 1 << " start." << endl;
#pragma omp parallel for num_threads(8)
		for (int n = 0; n < m_num_nodes; ++n) {
			if (mp_node[n].bdr == 1) {
				continue;
			}
			//节点内部迭代过程
			int maxNRitersteps = 6;
			double Ai = 0;
			for (int NRiter = 0; NRiter < maxNRitersteps; ++NRiter) {
				//cout << "n: " << n <<", NRiter: " << NRiter << endl;
				double S = 0, F = 0;
				double J = 0, Fj = 0;
				//装配过程
				for (int k = 0; k < mp_node[n].NumberofNeighbourElement; ++k) {
					//vector<vector<double>> J(3, vector<double>(3, 0));	//单元Jacobi矩阵
					//vector<double> Fj(3, 0);	//新增的右侧项
					int i_tri = mp_node[n].NeighbourElementId[k];
					CTriElement triele = mp_triele[i_tri];
					int nodenumber = mp_node[n].NeighbourElementNumber[k];
					double mu = triele.material->getMu(triele.B);
					//printf("mu: %f\n", mut);
					//处理线性单元
					if (triele.material->getLinearFlag() == true) {
						for (int i = 0; i < 3; ++i) {
							double Se = triele.C[nodenumber][i] / mu;
							if (nodenumber == i) {
								S += Se;
								F += triele.J * triele.area / 3;
								////永磁部分
								//double h_c = triele.material->getH_c();
								//double theta_m = triele.material->getTheta_m();
								//F += h_c / 2 * (triele.R[i] * cos(theta_m) - triele.Q[i] * sin(theta_m));
							}
							else {
								F -= Se * A_old[triele.n[i]];
							}
						}
					}
					//处理非线性单元
					else {
						double mu, dvdb, B, sigmai[3]{0, 0, 0};
						mu = triele.material->getMu(triele.B);
						dvdb = triele.material->getdvdB(triele.B);
						B = triele.B;
						for (int i = 0; i < 3; ++i) {
							for (int m = 0; m < 3; ++m) {
								if (m == nodenumber) {
									sigmai[i] += triele.C[i][m] * Ai;
								}
								else {
									sigmai[i] += triele.C[i][m] * A_old[triele.n[m]];
								}
							}
						}

						for (int i = 0; i < 3; ++i) {
							if (B != 0) {
								J = triele.C[nodenumber][i] / mu + dvdb * sigmai[nodenumber] * sigmai[i] / B / triele.area;
							}
							else {
								J = triele.C[nodenumber][i] / mu;
							}
							if (nodenumber == i) {
								S += J;
								F += (J - triele.C[nodenumber][i] / mu) * Ai;
							}
							else {
								F += (J - triele.C[nodenumber][i] / mu) * A_old[triele.n[i]];
								F -= J * A_old[triele.n[i]];
							}
						}
						//if (sigmai != 0 && sigmaj != 0) 
						//	cout << "sigmai: " << sigmai << ", sigmaj: " << sigmaj << endl;
					}

				}
				Ai = F / S;
				//cout << "F: " << F << ", S: " << S << endl;

				//NR迭代收敛性判断 
				double a = (Ai - mp_node[n].A) * (Ai - mp_node[n].A);
				double b = Ai * Ai;
				double NRerror = sqrt(a) / sqrt(b);
				if (Ai == 0) {
					break;
				}
				if (NRerror > maxerror) {
					mp_node[n].A = Ai;
					for (int i = 0; i < mp_node[n].NumberofNeighbourElement; ++i) {
						updateB(mp_node[n].NeighbourElementId[i]);
					}
				}
				else {
					//if (NRiter != 1) {
					//	printf("NRerror: %.20f\n", NRerror);
					//	cout << "NRerror: " << NRerror << endl;
					//	printf("n: %d, NRiter: %d\n", n, NRiter);
					//}
					break;
				}
			}
		}

		//判断全局收敛性
		double error = 0, a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (mp_node[i].A - A_old[i]) * (mp_node[i].A - A_old[i]);
			b += mp_node[i].A * mp_node[i].A;
		}
		error = sqrt(a) / sqrt(b);
		//if ((iter + 1) % 100 == 0) {
		//	cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
		//}
		cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
		if (error > maxerror) {
			for (int i = 0; i < m_num_nodes; ++i) {
				A_old[i] = mp_node[i].A;
			}
		}
		else {
			cout << "Iteration step: " << iter + 1 << endl;
			cout << "Nonlinear NDDR iteration finish.\n";
			return;
		}
//#pragma omp parallel for num_threads(8)
//		for (int i = 0; i < m_num_nodes; ++i) {
//			A_old[i] = mp_node[i].A;
//		}
	}
}

void FEM2DNDDRSolver::solve2DPlane1()
{
	for (int iter = 0; iter < maxitersteps; ++iter) {
#pragma omp parallel for num_threads(8)
		for (int i = 0; i < m_num_nodes; ++i) {
			double A, dA, Jac, k, J0, k0, NRerror, Js;
			double B2, B, V, VB2, B2A;
			int ID;
			int RelaxCount, count = 0;
			int NRCount = 2;
			double RelaxFactor = 1;
			if (mp_node[i].bdr != 1)
			{
				A = 0; dA = 0; NRerror = 1; count = 0;
				while (count < NRCount)
				{
					Jac = 0; k = 0;
					for (int j = 0; j < mp_node[i].NumberofNeighbourElement; j++)
					{
						ID = mp_node[i].NeighbourElementId[j];
						Js = mp_triele[ID].J * mp_triele[ID].area / 3;
						double ALocal[3]{ 0, 0, 0 };
						int nodenumber = mp_node[i].NeighbourElementNumber[j];
						for (int m = 0; m < 3; ++m) {
							if (m == nodenumber) {
								ALocal[m] = A;
							}
							else {
								ALocal[m] = mp_node[mp_triele[ID].n[m]].A_old;
							}
						}

						double RHSContri = 0;
						for (int m = 0; m < 3; ++m) {

							RHSContri += mp_triele[ID].C[nodenumber][m] * ALocal[m];
						}

						//线性单元，计算右侧列向量时，也要考虑A？？？
						if (mp_triele[ID].material->getLinearFlag() == true)	//线性空气单元
						{
							V = 1.0 / mp_triele[ID].material->getMu(0);
							J0 = V * mp_triele[ID].C[nodenumber][nodenumber];
							k0 = Js - V * RHSContri;

						}

						else	//非线性单元
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (ALocal[0] - ALocal[1]) * (ALocal[0] - ALocal[1]) + mp_triele[ID].C[1][2] * (ALocal[1] - ALocal[2]) * (ALocal[1] - ALocal[2]) + mp_triele[ID].C[0][2] * (ALocal[0] - ALocal[2]) * (ALocal[0] - ALocal[2]));
							B = sqrt(B2);	//计算B，用于处理非线性
							V = 1.0 / mp_triele[ID].material->getMu(B);
							VB2 = mp_triele[ID].material->getdvdB2(B);
							B2A = -1 / mp_triele[ID].area * (2 * mp_triele[ID].C[nodenumber][0] * (A - ALocal[0]) + 2 * mp_triele[ID].C[nodenumber][1] * (A - ALocal[1]) + 2 * mp_triele[ID].C[nodenumber][2] * (A - ALocal[2]));
							J0 = B2A * VB2 * RHSContri + V * mp_triele[ID].C[nodenumber][nodenumber];
							k0 = Js - V * RHSContri;
						}

						Jac = Jac + J0;
						k = k + k0;
					}
					dA = k / Jac;
					A = A + RelaxFactor * dA;

					double a = dA * dA;
					double b = A * A;
					double NRerror = sqrt(a) / sqrt(b);
					if (A == 0) {
						break;
					}
					if (NRerror > maxerror) {

						continue;
					}
					else {
						//if (NRiter != 1) {
						//	printf("NRerror: %.20f\n", NRerror);
						//	cout << "NRerror: " << NRerror << endl;
						//	printf("n: %d, NRiter: %d\n", n, NRiter);
						//}
						break;
					}
					count++;
				}

				mp_node[i].A = A;
			}
		}

		//判断全局收敛性
		double error = 0, a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (mp_node[i].A - mp_node[i].A_old) * (mp_node[i].A - mp_node[i].A_old);
			b += mp_node[i].A * mp_node[i].A;
		}
		error = sqrt(a) / sqrt(b);
		//if ((iter + 1) % 100 == 0) {
		//	cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
		//}
		cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
		if (error > maxerror) {
			//这一步影响效率
			for (int i = 0; i < m_num_nodes; ++i) {
				mp_node[i].A_old = mp_node[i].A;
			}
		}
		else {
			cout << "Iteration step: " << iter + 1 << endl;
			cout << "Nonlinear NDDR iteration finish.\n";
			return;
		}
	}
}

void FEM2DNDDRSolver::solve2DPlane2()
{
	for (int RelaxCount = 0; RelaxCount < maxitersteps/2; RelaxCount++) {
		Update_Magnetic_Node_A();


		double a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (mp_node[i].A - mp_node[i].A_old) * (mp_node[i].A - mp_node[i].A_old);
			b += mp_node[i].A * mp_node[i].A;
		}
		double error = sqrt(a) / sqrt(b);
		cout << "step: " << 2 * RelaxCount << ", error: " << error << endl;

		Update_Magnetic_Node_A_old();

		a = 0; b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (mp_node[i].A - mp_node[i].A_old) * (mp_node[i].A - mp_node[i].A_old);
			b += mp_node[i].A * mp_node[i].A;
		}
		error = sqrt(a) / sqrt(b);
		cout << "step: " << 2 * RelaxCount  + 1 << ", error: " << error << endl;
	}

	//updateB();
}

void FEM2DNDDRSolver::Update_Magnetic_Node_A()
{
#pragma omp parallel for num_threads(8)
	for (int i = 0; i < m_num_nodes; ++i) {
		double A, dA, Jac, k, J0, k0, NRerror, Js;
		double A1, A2, A3;
		double B2, B, V, VB2, B2A;
		int I, J, K, ID;
		int RelaxCount, count = 0;
		int NRCount = 6;
		double RelaxFactor = 1;
		if (mp_node[i].bdr != 1)
		{
			A = 0; dA = 0; NRerror = 1; count = 0;
			while (count < NRCount)
			{
				Jac = 0; k = 0;
				for (int j = 0; j < mp_node[i].NumberofNeighbourElement; j++)
				{
					ID = mp_node[i].NeighbourElementId[j];
					Js = mp_triele[ID].J * mp_triele[ID].area / 3;
					I = mp_triele[ID].n[0];
					J = mp_triele[ID].n[1];
					K = mp_triele[ID].n[2];
					A1 = mp_node[I].A;
					A2 = mp_node[J].A;
					A3 = mp_node[K].A;
					//线性单元，计算右侧列向量时，也要考虑A？？？
					if (mp_triele[ID].material->getLinearFlag() == true)	//线性空气单元
					{
						double Ve = 1.0 / mp_triele[ID].material->getMu(0);
						if (mp_node[i].NeighbourElementNumber[j] == 0)
						{
							J0 = Ve * mp_triele[ID].C[0][0];
							k0 = Js - A * Ve * mp_triele[ID].C[0][0] - A2 * Ve * mp_triele[ID].C[0][1] - A3 * Ve * mp_triele[ID].C[0][2];
						}
						if (mp_node[i].NeighbourElementNumber[j] == 1)
						{
							J0 = Ve * mp_triele[ID].C[1][1];
							k0 = Js - A1 * Ve * mp_triele[ID].C[1][0] - A * Ve * mp_triele[ID].C[1][1] - A3 * Ve * mp_triele[ID].C[1][2];
						}
						if (mp_node[i].NeighbourElementNumber[j] == 2)
						{
							J0 = Ve * mp_triele[ID].C[2][2];
							k0 = Js - A1 * Ve * mp_triele[ID].C[2][0] - A2 * Ve * mp_triele[ID].C[2][1] - A * Ve * mp_triele[ID].C[2][2];
						}
					}
					else	//非线性单元
					{
						double Mu0 = 0.002;
						if (mp_node[i].NeighbourElementNumber[j] == 0)
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (A - A2) * (A - A2) + mp_triele[ID].C[1][2] * (A2 - A3) * (A2 - A3) + mp_triele[ID].C[0][2] * (A - A3) * (A - A3));
							B = sqrt(B2);	//计算B，用于处理非线性
							if (B < 0.6)
							{
								V = 1 / Mu0;
								VB2 = 0;
								J0 = V * mp_triele[ID].C[0][0];
								k0 = Js - A * V * mp_triele[ID].C[0][0] - A2 * V * mp_triele[ID].C[0][1] - A3 * V * mp_triele[ID].C[0][2];

							}
							else
							{
								V = 1 / Mu0 + 3000 * (B - 0.6) * (B - 0.6)* (B - 0.6) / B;
								VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6)* (B - 0.6)) / B / B / 2 / B;
								B2A = -1 / mp_triele[ID].area * (2 * mp_triele[ID].C[0][1] * (A - A2) + 2 * mp_triele[ID].C[0][2] * (A - A3));//dfrac{dB^2}{dA}?
								J0 = B2A * VB2 * (mp_triele[ID].C[0][0] * A + mp_triele[ID].C[0][1] * A2 + mp_triele[ID].C[0][2] * A3) + V * mp_triele[ID].C[0][0];
								k0 = Js - V * (mp_triele[ID].C[0][0] * A + mp_triele[ID].C[0][1] * A2 + mp_triele[ID].C[0][2] * A3);
							}

						}
						else if (mp_node[i].NeighbourElementNumber[j] == 1)
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (A1 - A) * (A1 - A) + mp_triele[ID].C[1][2] * (A - A3) * (A - A3) + mp_triele[ID].C[0][2] * (A1 - A3) * (A1 - A3));
							B = sqrt(B2);
							if (B < 0.6)
							{
								V = 1 / Mu0;
								VB2 = 0;
								J0 = V * mp_triele[ID].C[1][1];
								k0 = Js - A1 * V * mp_triele[ID].C[1][0] - A * V * mp_triele[ID].C[1][1] - A3 * V * mp_triele[ID].C[1][2];

							}
							else
							{
								V = 1 / Mu0 + 3000 * (B - 0.6) * (B - 0.6)* (B - 0.6) / B;
								VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6)* (B - 0.6)) / B / B / 2 / B;
								B2A = -1 / mp_triele[ID].area * (2 * mp_triele[ID].C[0][1] * (A - A1) + 2 * mp_triele[ID].C[1][2] * (A - A3));
								J0 = B2A * VB2 * (mp_triele[ID].C[1][0] * A1 + mp_triele[ID].C[1][1] * A + mp_triele[ID].C[1][2] * A3) + V * mp_triele[ID].C[1][1];
								k0 = Js - V * (mp_triele[ID].C[1][0] * A1 + mp_triele[ID].C[1][1] * A + mp_triele[ID].C[1][2] * A3);
							}
						}

						else if (mp_node[i].NeighbourElementNumber[j] == 2)
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (A1 - A2) * (A1 - A2) + mp_triele[ID].C[1][2] * (A - A2) * (A - A2) + mp_triele[ID].C[0][2] * (A - A1) * (A - A1));
							B = sqrt(B2);
							if (B < 0.6)
							{
								V = 1 / Mu0;
								VB2 = 0;
								J0 = V * mp_triele[ID].C[2][2];
								k0 = Js - A1 * V * mp_triele[ID].C[2][0] - A2 * V * mp_triele[ID].C[2][1] - A * V * mp_triele[ID].C[2][2];

							}
							else
							{
								V = 1 / Mu0 + 3000 * (B - 0.6) * (B - 0.6)* (B - 0.6) / B;
								VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6)* (B - 0.6)) / B / B / 2 / B;
								B2A = -1 / mp_triele[ID].area * (2 * mp_triele[ID].C[1][2] * (A - A2) + 2 * mp_triele[ID].C[2][0] * (A - A1));
								J0 = B2A * VB2 * (mp_triele[ID].C[2][0] * A1 + mp_triele[ID].C[2][1] * A2 + mp_triele[ID].C[2][2] * A) + V * mp_triele[ID].C[2][2];
								k0 = Js - V * (mp_triele[ID].C[2][0] * A1 + mp_triele[ID].C[2][1] * A2 + mp_triele[ID].C[2][2] * A);
							}
						}
					}

					Jac = Jac + J0;
					k = k + k0;
				}
				dA = k / Jac;
				A = A + RelaxFactor * dA;

				double a = dA * dA;
				double b = A * A;
				double NRerror = sqrt(a) / sqrt(b);
				if (A == 0) {
					break;
				}
				if (NRerror > maxerror) {
					continue;
				}
				else {
					//if (NRiter != 1) {
					//	printf("NRerror: %.20f\n", NRerror);
					//	cout << "NRerror: " << NRerror << endl;
					//	printf("n: %d, NRiter: %d\n", n, NRiter);
					//}
					break;
				}
				count++;
			}

			mp_node[i].A_old = A;
		}
	}

}

void FEM2DNDDRSolver::Update_Magnetic_Node_A_old()
{
#pragma omp parallel for num_threads(8)
	for (int i = 0; i < m_num_nodes; ++i) {
		double A, dA, Jac, k, J0, k0, NRerror, Js;
		double A1, A2, A3;
		double B2, B, V, VB2, B2A;
		int I, J, K, ID;
		int RelaxCount, count = 0;
		int NRCount = 6;
		double RelaxFactor = 1;
		if (mp_node[i].bdr != 1)
		{
			A = 0; dA = 0; NRerror = 1; count = 0;
			while (count < NRCount)
			{
				Jac = 0; k = 0;
				for (int j = 0; j < mp_node[i].NumberofNeighbourElement; j++)
				{
					ID = mp_node[i].NeighbourElementId[j];
					Js = mp_triele[ID].J * mp_triele[ID].area / 3;
					I = mp_triele[ID].n[0];
					J = mp_triele[ID].n[1];
					K = mp_triele[ID].n[2];
					A1 = mp_node[I].A_old;
					A2 = mp_node[J].A_old;
					A3 = mp_node[K].A_old;
					//线性单元，计算右侧列向量时，也要考虑A？？？
					if (mp_triele[ID].material->getLinearFlag() == true)	//线性空气单元
					{
						double Ve = 1.0 / mp_triele[ID].material->getMu(0);
						if (mp_node[i].NeighbourElementNumber[j] == 0)
						{
							J0 = Ve * mp_triele[ID].C[0][0];
							k0 = Js - A * Ve * mp_triele[ID].C[0][0] - A2 * Ve * mp_triele[ID].C[0][1] - A3 * Ve * mp_triele[ID].C[0][2];
						}
						if (mp_node[i].NeighbourElementNumber[j] == 1)
						{
							J0 = Ve * mp_triele[ID].C[1][1];
							k0 = Js - A1 * Ve * mp_triele[ID].C[1][0] - A * Ve * mp_triele[ID].C[1][1] - A3 * Ve * mp_triele[ID].C[1][2];
						}
						if (mp_node[i].NeighbourElementNumber[j] == 2)
						{
							J0 = Ve * mp_triele[ID].C[2][2];
							k0 = Js - A1 * Ve * mp_triele[ID].C[2][0] - A2 * Ve * mp_triele[ID].C[2][1] - A * Ve * mp_triele[ID].C[2][2];
						}
					}
					else	//非线性单元
					{
						double Mu0 = 0.002;
						if (mp_node[i].NeighbourElementNumber[j] == 0)
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (A - A2) * (A - A2) + mp_triele[ID].C[1][2] * (A2 - A3) * (A2 - A3) + mp_triele[ID].C[0][2] * (A - A3) * (A - A3));
							B = sqrt(B2);	//计算B，用于处理非线性
							if (B < 0.6)
							{
								V = 1 / Mu0;
								VB2 = 0;
								J0 = V * mp_triele[ID].C[0][0];
								k0 = Js - A * V * mp_triele[ID].C[0][0] - A2 * V * mp_triele[ID].C[0][1] - A3 * V * mp_triele[ID].C[0][2];

							}
							else
							{
								V = 1 / Mu0 + 3000 * (B - 0.6) * (B - 0.6) * (B - 0.6) / B;
								VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6) * (B - 0.6)) / B / B / 2 / B;
								B2A = -1 / mp_triele[ID].area * (2 * mp_triele[ID].C[0][1] * (A - A2) + 2 * mp_triele[ID].C[0][2] * (A - A3));//dfrac{dB^2}{dA}?
								J0 = B2A * VB2 * (mp_triele[ID].C[0][0] * A + mp_triele[ID].C[0][1] * A2 + mp_triele[ID].C[0][2] * A3) + V * mp_triele[ID].C[0][0];
								k0 = Js - V * (mp_triele[ID].C[0][0] * A + mp_triele[ID].C[0][1] * A2 + mp_triele[ID].C[0][2] * A3);

								//V = 1 / Mu0 + 1e5 * powf(B - 1.3, 2) / B;
								//VB2 = (B * 2e5 * powf(B - 1.3, 1) - 1e5 * powf(B - 1.3, 2)) / B / B / 2 / B;	//dfrac{dv}{d(B^2)}?
								//B2A = -1 / MyElem_d[ID].Area * (2 * MyElem_d[ID].k01 * (A - A2) + 2 * MyElem_d[ID].k02 * (A - A3));//dfrac{dB^2}{dA}?
								//J0 = B2A * VB2 * (MyElem_d[ID].k00 * A + MyElem_d[ID].k01 * A2 + MyElem_d[ID].k02 * A3) + V * MyElem_d[ID].k00;
								//k0 = Js - V * (MyElem_d[ID].k00 * A + MyElem_d[ID].k01 * A2 + MyElem_d[ID].k02 * A3);
							}

						}
						else if (mp_node[i].NeighbourElementNumber[j] == 1)
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (A1 - A) * (A1 - A) + mp_triele[ID].C[1][2] * (A - A3) * (A - A3) + mp_triele[ID].C[0][2] * (A1 - A3) * (A1 - A3));
							B = sqrt(B2);
							if (B < 0.6)
							{
								V = 1 / Mu0;
								VB2 = 0;
								J0 = V * mp_triele[ID].C[1][1];
								k0 = Js - A1 * V * mp_triele[ID].C[1][0] - A * V * mp_triele[ID].C[1][1] - A3 * V * mp_triele[ID].C[1][2];

							}
							else
							{
								V = 1 / Mu0 + 3000 * (B - 0.6) * (B - 0.6) * (B - 0.6) / B;
								VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6) * (B - 0.6)) / B / B / 2 / B;
								B2A = -1 / mp_triele[ID].area * (2 * mp_triele[ID].C[0][1] * (A - A1) + 2 * mp_triele[ID].C[1][2] * (A - A3));
								J0 = B2A * VB2 * (mp_triele[ID].C[1][0] * A1 + mp_triele[ID].C[1][1] * A + mp_triele[ID].C[1][2] * A3) + V * mp_triele[ID].C[1][1];
								k0 = Js - V * (mp_triele[ID].C[1][0] * A1 + mp_triele[ID].C[1][1] * A + mp_triele[ID].C[1][2] * A3);

								//V = 1 / Mu0 + 1e5 * powf(B - 1.3, 2) / B;
								//VB2 = (B * 2e5 * powf(B - 1.3, 1) - 1e5 * powf(B - 1.3, 2)) / B / B / 2 / B;
								//B2A = -1 / MyElem_d[ID].Area * (2 * MyElem_d[ID].k01 * (A - A1) + 2 * MyElem_d[ID].k12 * (A - A3));
								//J0 = B2A * VB2 * (MyElem_d[ID].k10 * A1 + MyElem_d[ID].k11 * A + MyElem_d[ID].k12 * A3) + V * MyElem_d[ID].k11;
								//k0 = Js - V * (MyElem_d[ID].k10 * A1 + MyElem_d[ID].k11 * A + MyElem_d[ID].k12 * A3);
							}
						}

						else if (mp_node[i].NeighbourElementNumber[j] == 2)
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (A1 - A2) * (A1 - A2) + mp_triele[ID].C[1][2] * (A - A2) * (A - A2) + mp_triele[ID].C[0][2] * (A - A1) * (A - A1));
							B = sqrt(B2);
							if (B < 0.6)
							{
								V = 1 / Mu0;
								VB2 = 0;
								J0 = V * mp_triele[ID].C[2][2];
								k0 = Js - A1 * V * mp_triele[ID].C[2][0] - A2 * V * mp_triele[ID].C[2][1] - A * V * mp_triele[ID].C[2][2];

							}
							else
							{
								V = 1 / Mu0 + 3000 * (B - 0.6) * (B - 0.6) * (B - 0.6) / B;
								VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6) * (B - 0.6)) / B / B / 2 / B;
								B2A = -1 / mp_triele[ID].area * (2 * mp_triele[ID].C[1][2] * (A - A2) + 2 * mp_triele[ID].C[2][0] * (A - A1));
								J0 = B2A * VB2 * (mp_triele[ID].C[2][0] * A1 + mp_triele[ID].C[2][1] * A2 + mp_triele[ID].C[2][2] * A) + V * mp_triele[ID].C[2][2];
								k0 = Js - V * (mp_triele[ID].C[2][0] * A1 + mp_triele[ID].C[2][1] * A2 + mp_triele[ID].C[2][2] * A);

								//V = 1 / Mu0 + 1e5 * powf(B - 1.3, 2) / B;
								//VB2 = (B * 2e5 * powf(B - 1.3, 1) - 1e5 * powf(B - 1.3, 2)) / B / B / 2 / B;
								//B2A = -1 / MyElem_d[ID].Area * (2 * MyElem_d[ID].k12 * (A - A2) + 2 * MyElem_d[ID].k20 * (A - A1));
								//J0 = B2A * VB2 * (MyElem_d[ID].k20 * A1 + MyElem_d[ID].k21 * A2 + MyElem_d[ID].k22 * A) + V * MyElem_d[ID].k22;
								//k0 = Js - V * (MyElem_d[ID].k20 * A1 + MyElem_d[ID].k21 * A2 + MyElem_d[ID].k22 * A);
							}
						}
					}

					Jac = Jac + J0;
					k = k + k0;
				}
				dA = k / Jac;
				A = A + RelaxFactor * dA;

				double a = dA * dA;
				double b = A * A;
				double NRerror = sqrt(a) / sqrt(b);
				if (A == 0) {
					break;
				}
				if (NRerror > maxerror) {
					continue;
				}
				else {
					//if (NRiter != 1) {
					//	printf("NRerror: %.20f\n", NRerror);
					//	cout << "NRerror: " << NRerror << endl;
					//	printf("n: %d, NRiter: %d\n", n, NRiter);
					//}
					break;
				}

				count++;
			}

			mp_node[i].A = A;
		}
	}
}
