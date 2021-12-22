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
	//���������ε�Ԫ���β���
	makeTrangle();
	//����߽����������Ϻ͸���
	processBoundaryCondition();
	processMaterial();
	processLoad();
	//����NDDR�Ľڵ���Ϣ
	for (int i = 0; i < m_num_nodes; ++i) {
		mp_node[i].NumberofNeighbourElement = 0;
	}
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		for (int j = 0; j < 3; ++j) {
			int n = mp_triele[i_tri].n[j];
			int id = mp_node[n].NumberofNeighbourElement;
			mp_node[n].NeighbourElementId[id] = i_tri;
			mp_node[n].NeighbourElementNumber[id] = j;
			mp_node[n].NumberofNeighbourElement++;
		}
	}

	//�����ڲ�����ţ�ٵ���
	clock_t start, end;
	start = clock();
	if (dimension == FEMModel::DIMENSION::D2AXISM) {
		//solve2DAxim1();	//����汾��NDDR�㷨
		//solve2DAxim2();	//���ò�ֵNR������汾NDDR�㷨
		//solve2DAximOpt();	//��һ���Ż����NDDR�㷨������㷨���ڵ����⣺Jacobi����ȡ�õĽⲢ����ȷ������NR�����Ĵ��������ֱ�ӷ�Ҫ�಻��
		//solve2DAximOpt1();	//�ڶ����Ż�NDDR�㷨���������Է������ϵ������С�
		//solve2DAximPrecondition();	//Ԥ�����Ż�NDDR����Ȼ��Ч��Ż��
		solve2DAximRobin();
	}
	else if (dimension == FEMModel::DIMENSION::D2PLANE) {
		//solve2DPlane();
		//solve2DPlane1();
		solve2DPlane2();
	}
	end = clock();
	cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
}

  void FEM2DNDDRSolver::solveDynamic()
{
	//COMSOL��̬����
	string name = "RelayDynamic";

	const int n = 101;
	current = new double[n];
	dis = new double[n];
	velocity = new double[n];
	acc = new double[n];
	magneticforce = new double[n];
	springforce = new double[n];
	flux = new double[n];
	mass = 0.024;
	h = 5e-4;
	U = 24;
	R = 40;
	int dynamicstepsarray[n];
	double time[n];

	dis[0] = 0;
	velocity[0] = 0;
	acc[0] = 0;
	springforce[0] = solveSpringForce(4, 0);
	magneticforce[0] = 0;
	current[0] = 0;
	flux[0] = 0;
	dynamicstepsarray[0] = 0;
	time[0] = 0;

	bool stopflag = false;


	clock_t start, end;
	start = clock();
	for (int i = 1; i < n; ++i) {
		clock_t start, end;
		start = clock();
		dynamicsteps = 0;
		cout << "solve step " << i << "...\n";

		acc[i] = (magneticforce[i - 1] + springforce[i - 1]) / mass;
		if (acc[i] < 0) acc[i] = 0;
		velocity[i] = velocity[i - 1] + h * acc[i - 1];
		dis[i] = dis[i - 1] + h * velocity[i - 1];

		if (dis[i] >= 0.0024) {
			dis[i] = 0.0024;
			if (stopflag == false) {
				stopflag = true;
			}
			else {
				velocity[i] = 0;
				acc[i] = 0;
			}
		}



		//���η�����
		//springforce[i] = solveSpringForce(1, dis[i]);
		//acc[i] = (magneticforce[i - 1] + springforce[i - 1]) / mass;
		//if (acc[i] < 0) acc[i] = 0;
		//velocity[i] = velocity[i - 1] + 0.5 * h * (acc[i] + acc[i - 1]);
		//dis[i] = dis[i - 1] + 0.5 * h * (velocity[i] + velocity[i - 1]);

		//����ŷ��������
		springforce[i] = solveSpringForce(4, dis[i]);
		acc[i] = (magneticforce[i - 1] + springforce[i - 1]) / mass;
		if (acc[i] < 0) acc[i] = 0;
		velocity[i] = velocity[i - 1] + h * acc[i];
		dis[i] = dis[i - 1] + h * velocity[i];

		if (dis[i] >= 0.0024) {
			dis[i] = 0.0024;
			if (stopflag == false) {
				stopflag = true;
			}
			else {
				velocity[i] = 0;
				acc[i] = 0;
			}
		}

		//��ǰλ��ʱ�Ĵų�-��·���
		//meshmanager->remesh(name, i, 0, dis[i] - dis[i - 1]);
		if (/*dis[i] - dis[i - 1] != 0*/i >= 30 && i <= 46) {
			string meshfile = "D:/femplatform/model/geo/modelcomsol_dynamic_NDDR/modelwithband_";
			meshfile += to_string(i) + ".mphtxt";
			meshmanager->readMeshFile(meshfile);
			setNodes(meshmanager->getNumofNodes(), meshmanager->getNodes());
			setVtxElements(meshmanager->getNumofVtxEle(), meshmanager->getVtxElements());
			setEdgElements(meshmanager->getNumofEdgEle(), meshmanager->getEdgElements());
			setTriElements(meshmanager->getNumofTriEle(), meshmanager->getTriElements());
		}
		//updateLoadmap(3, current[i]);
		//solveStatic();
		solveWeakCouple(i);
		solveMagneticForce();
		magneticforce[i] = Fy;
		springforce[i] = solveSpringForce(4, dis[i]);

		dynamicstepsarray[i] = dynamicsteps;
		end = clock();
		time[i] = double(end - start) / CLOCKS_PER_SEC;
		printf("dynamicsteps: %d, time: %f\n", dynamicstepsarray[i], time[i]);
		printf("step: %d, dis: %f, velocity: %f, acc: %f, springforce: %f, magneticforce: %f\n\n", i, dis[i], velocity[i], acc[i], springforce[i], magneticforce[i]);


	}

	for (int i = 0; i < n; ++i) {
		printf("time: %f, dis: %f, velocity: %f, acc: %f, springforce: %f, magneticforce: %f, current: %f, flux: %f\n", i * h, dis[i], velocity[i], acc[i], springforce[i], magneticforce[i], current[i], flux[i]);
	}

	//д�����ļ�
	char fn[256];
	sprintf(fn, "%s.m", "RelayDynamic");
	FILE* fp;
	fp = fopen(fn, "w");
	fprintf(fp, "%%output by FEEM\n");
	fprintf(fp, "%%timesteps, displacements, velocities, accelerations, magneticforce, current, flux, steps, time\n");

	fprintf(fp, "results = [\n");
	for (int i = 0; i < n; ++i) {
		fprintf(fp, "%10.8e,", i * h);
		fprintf(fp, "%10.8e,", dis[i]);
		fprintf(fp, "%10.8e,", velocity[i]);
		fprintf(fp, "%10.8e,", acc[i]);
		fprintf(fp, "%10.8e,", magneticforce[i]);
		fprintf(fp, "%10.8e,", current[i]);
		fprintf(fp, "%10.8e,", flux[i]);
		fprintf(fp, "%d,", dynamicstepsarray[i]);
		fprintf(fp, "%10.8e,", time[i]);
		fprintf(fp, "; \n");
	}
	fprintf(fp, "];\n");

	fprintf(fp, "subplot(2,3,1);hold on;\n");
	fprintf(fp, "plot(results(:,1),results(:,2),'*-');\n");
	//fprintf(fp, "plot(results(:,1),results(:,3),'*-');\n");
	fprintf(fp, "title(\"%s\");\n\n", "displacement");

	fprintf(fp, "subplot(2,3,2);hold on;\n");
	fprintf(fp, "plot(results(:,1),results(:,3),'*-');\n");
	fprintf(fp, "title(\"%s\");\n\n", "velocities");

	fprintf(fp, "subplot(2,3,3);hold on;\n");
	fprintf(fp, "plot(results(:,1),results(:,4),'*-');\n");
	fprintf(fp, "title(\"%s\");\n\n", "accelerations");

	fprintf(fp, "subplot(2,3,4);hold on;\n");
	fprintf(fp, "plot(results(:,1),results(:,5),'*-');\n");
	fprintf(fp, "title(\"%s\");\n\n", "magforces");

	fprintf(fp, "subplot(2,3,5);hold on;\n");
	fprintf(fp, "plot(results(:,1),results(:,6),'*-');\n");
	//fprintf(fp, "plot(results(:,1),results(:,9),'*-');\n");
	fprintf(fp, "title(\"%s\");\n\n", "ICoil");

	fprintf(fp, "subplot(2,3,6);hold on;\n");
	fprintf(fp, "plot(results(:,1),results(:,7),'*-');\n");
	fprintf(fp, "title(\"%s\");\n", "PhiCoil");

	fclose(fp);

	delete[] flux;
	delete[] springforce;
	delete[] magneticforce;
	delete[] acc;
	delete[] velocity;
	delete[] dis;
	delete[] current;

	//string name = "RelayDynamic";
	//cout << "void FEM2DNRSolver::solveDynamic()" << endl;
	//const int n = 11;
	//double dis[n] = { 0, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00015 };
	//double force[n + 1];
	//vector<int> air_domain{ 4, 5, 6, 7 };
	//for (int i = 0; i < n; ++i) {
	//	cout << "solve step " << i + 1 << "...\n";
	//	//meshmanager->remesh(name, i, 0, dis[i]);
	//	string meshfilename = "D:/femplatform/model/geo/modelcomsol_static/modelwithband_";
	//	meshfilename += to_string(i) + ".mphtxt";
	//	cout << "meshfilename: " << meshfilename << endl;
	//	meshmanager->readMeshFile(meshfilename);
	//	setNodes(meshmanager->getNumofNodes(), meshmanager->getNodes());
	//	setVtxElements(meshmanager->getNumofVtxEle(), meshmanager->getVtxElements());
	//	setEdgElements(meshmanager->getNumofEdgEle(), meshmanager->getEdgElements());
	//	setTriElements(meshmanager->getNumofTriEle(), meshmanager->getTriElements());
	//	solveStatic();
	//	solveMagneticForce();
	//	//solveMagneticForce1();
	//	force[i] = Fy;
	//	writeVtkFile(name + "_" + to_string(i));
	//	writeVtkFileNoAir(name + "_" + to_string(i), air_domain);
	//	cout << "step " << i + 1 << " solve finish.\n\n";
	//}

	//for (int i = 0; i < n; ++i) {
	//	cout << "i: " << i << ", Fy: " << force[i] << endl;
	//}
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
//			//�ڵ��ڲ���������
//			double A, dA, Jac, ku, J0, k0, NRerror, Js;
//			double A1, A2, A3;
//			double B2, B, V, VB2, B2A;
//			int I, J, K, ID;
//			int i, RelaxCount, count = 0;
//			double maxNRitersteps = 6;
//			double h_c, theta_m;
//			for (int NRiter = 0; NRiter < maxNRitersteps; ++NRiter) {
//				Jac = 0, ku = 0;
//				//װ�����
//				for (int k = 0; k < mp_node[n].NumberofNeighbourElement; ++k) {
//					//vector<vector<double>> J(3, vector<double>(3, 0));	//��ԪJacobi����
//					//vector<double> Fj(3, 0);	//�������Ҳ���
//					int i_tri = mp_node[n].NeighbourElementId[k];
//					CTriElement triele = mp_triele[i_tri];
//					int nodenumber = mp_node[n].NeighbourElementNumber[k];
//					double mut = triele.material->getMu(triele.B) * triele.xdot;
//					//printf("mu: %f\n", mut);
//					//�������Ե�Ԫ
//					if (triele.material->getLinearFlag() == true) {
//						for (int i = 0; i < 3; ++i) {
//							double Se = triele.C[nodenumber][i] / mut;
//							if (nodenumber == i) {
//								J0 = Se;
//								k0 = triele.J * triele.area / 3;
//								//���Ų���
//								h_c = triele.material->getH_c();
//								theta_m = triele.material->getTheta_m();
//								k0 += h_c / 2 * (triele.R[i] * cos(theta_m) - triele.Q[i] * sin(theta_m));
//							}
//							else {
//								k0 -= Se * mp_node[triele.n[i]].At_old;
//							}
//						}
//					}
//					//��������Ե�Ԫ
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
//				//NR�����������ж� 
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
//		//�ж�ȫ��������
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
			//�ڵ��ڲ���������
			int maxNRitersteps = 100 ;
			double Ati = 0;
			for (int NRiter = 0; NRiter < maxNRitersteps; ++NRiter) {
				double S = 0, F = 0;
				double J = 0, Fj = 0;
				//װ�����
				for (int k = 0; k < mp_node[n].NumberofNeighbourElement; ++k) {
					//vector<vector<double>> J(3, vector<double>(3, 0));	//��ԪJacobi����
					//vector<double> Fj(3, 0);	//�������Ҳ���
					int i_tri = mp_node[n].NeighbourElementId[k];
					CTriElement triele = mp_triele[i_tri];
					int nodenumber = mp_node[n].NeighbourElementNumber[k];
					double mut = triele.material->getMu(triele.B) * triele.xdot;
					//printf("mu: %f\n", mut);
					//�������Ե�Ԫ
					if (triele.material->getLinearFlag() == true) {
						for (int i = 0; i < 3; ++i) {
							double Se = triele.C[nodenumber][i] / mut;
							if (nodenumber == i) {
								S += Se;
								F += triele.J * triele.area / 3;
								//���Ų���
								double h_c = triele.material->getH_c();
								double theta_m = triele.material->getTheta_m();
								F += h_c / 2 * (triele.R[i] * cos(theta_m) - triele.Q[i] * sin(theta_m));
							}
							else {
								F -= Se * mp_node[triele.n[i]].At_old;
							}
						}
					}
					//��������Ե�Ԫ
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
				//NR�����������ж� 
				double a = (Ati - mp_node[n].At) * (Ati - mp_node[n].At);
				double b = Ati * Ati;
				double NRerror = sqrt(a) / sqrt(b);
				if (Ati == 0) {
					continue;
				}
				if (NRerror > 1e-6) {
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

		//�ж�ȫ��������
		double error = 0, a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (mp_node[i].At - mp_node[i].At_old) * (mp_node[i].At - mp_node[i].At_old);
			b += mp_node[i].At * mp_node[i].At;
		}
		error = sqrt(a) / sqrt(b);
		cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;

		//if ((iter + 1) % 100 == 0) {
		//	cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
		//}
		if (error > maxerror) {
			for (int i = 0; i < m_num_nodes; ++i) {
				mp_node[i].At_old = mp_node[i].At;
			}
		}
		else {
			staticsteps = iter;
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
//#pragma omp parallel for num_threads(8)
		for (int i = 0; i < m_num_nodes; ++i) {
			double A, At, dAt, Jac, k, J0, k0, NRerror, Js;
			double B2, B, V, Vt, VB2, B2A;
			int ID;
			int RelaxCount, count = 0;
			int NRCount = 1000;
			double RelaxFactor = 1;
			double AtLocal[3]{ 0 ,0, 0 }, ALocal[3]{0, 0, 0};
			if (mp_node[i].bdr != 1)
			{
				At = 0; dAt = 0; NRerror = 10; count = 0;
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

						//���Ե�Ԫ�������Ҳ�������ʱ��ҲҪ����A������
						if (mp_triele[ID].material->getLinearFlag() == true)	//���Կ�����Ԫ
						{
							Vt = 1.0 / mp_triele[ID].material->getMu(0) / mp_triele[ID].xdot;
							J0 = Vt * mp_triele[ID].C[nodenumber][nodenumber];
							k0 = Js - Vt * RHSContri;
						}

						else	//�����Ե�Ԫ
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (ALocal[0] - ALocal[1]) * (ALocal[0] - ALocal[1]) + mp_triele[ID].C[1][2] * (ALocal[1] - ALocal[2]) * (ALocal[1] - ALocal[2]) + mp_triele[ID].C[0][2] * (ALocal[0] - ALocal[2]) * (ALocal[0] - ALocal[2]));
							B = sqrt(B2);	//����B�����ڴ��������
							Vt = 1.0 / mp_triele[ID].material->getMu(B) / mp_triele[ID].xdot;
							VB2 = mp_triele[ID].material->getdvdB2(B);
							B2A = -1 / mp_triele[ID].area * (2 * mp_triele[ID].C[nodenumber][0] * (A - ALocal[0]) + 2 * mp_triele[ID].C[nodenumber][1] * (A - ALocal[1]) + 2 * mp_triele[ID].C[nodenumber][2] * (A - ALocal[2]));
							J0 = B2A * VB2 * RHSContri / mp_triele[ID].xdot / mp_triele[ID].xdot + Vt * mp_triele[ID].C[nodenumber][nodenumber];
							k0 = Js - Vt * RHSContri;
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
						//cout << "node: " << i << ", NRerror: " << NRerror << endl;
						count++;
						continue;
					}
					else {
						//cout << "node: " << i << ", NRerror: " << NRerror << endl;
						//if (NRiter != 1) {
						//	printf("NRerror: %.20f\n", NRerror);
						//	cout << "NRerror: " << NRerror << endl;
						//	printf("n: %d, NRiter: %d\n", n, NRiter);
						//}
						break;
					}

				}
				mp_node[i].At = At;
			}
		}

		//�ж�ȫ��������
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
			//��һ��Ӱ��Ч��
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

void FEM2DNDDRSolver::solve2DAximOpt()
{
	for (int iter = 0; iter < maxitersteps; ++iter) {
		//#pragma omp parallel for num_threads(8)
		for (int n = 0; n < m_num_nodes; ++n) {
			if (mp_node[n].bdr == 1) {
				continue;
			}

			double S = 0, F = 0, NA = 0;
			double M = 0;
			double Se[3] = { 0, 0, 0 };
			for (int k = 0; k < mp_node[n].NumberofNeighbourElement; ++k) {
				int i_tri = mp_node[n].NeighbourElementId[k];	//ͨ����Ԫ�������ڵı�ż�����Ԫ��ȫ�ֱ��
				CTriElement triele = mp_triele[i_tri];
				int nodenumber = mp_node[n].NeighbourElementNumber[k];	//���������Ľڵ��Ӧ�ĵ�Ԫ�ֲ����
				double mu = triele.material->getMu(triele.B);
				double mut = mu * triele.xdot;
				//--------------------------------------��Ԫ����
				if (triele.material->getLinearFlag() == true) {
					for (int i = 0; i < 3; ++i) {
						Se[i] = triele.C[nodenumber][i] / mut;
						if (i != nodenumber) {
							NA -= Se[i] * mp_node[triele.n[i]].delta_At_old;
						}
					}
					M += Se[nodenumber];
					F += (Se[0] * mp_node[triele.n[0]].At_old + Se[1] * mp_node[triele.n[1]].At_old + Se[2] * mp_node[triele.n[2]].At_old - triele.J * triele.area / 3);
				}
				else
				{
					//-------------------------------------���㵥Ԫ��Jacobi����
					double dvdb, dvdbt, Bt, sigmai[3]{ 0, 0, 0 }, J[3]{ 0, 0, 0 };
					dvdb = triele.material->getdvdB(triele.B);
					dvdbt = dvdb / triele.xdot / triele.xdot;
					Bt = triele.B * triele.xdot;
					for (int i = 0; i < 3; ++i) {
						Se[i] = triele.C[nodenumber][i] / mut;
						for (int m = 0; m < 3; ++m) {
							sigmai[i] += triele.C[i][m] * mp_node[triele.n[m]].At;
						}
						J[i] = triele.C[nodenumber][i] / mut;
						if (Bt != 0) {
							J[i] += dvdbt * sigmai[nodenumber] * sigmai[i] / Bt / triele.area;
						}
						if (i != nodenumber) {
							NA -= J[i] * mp_node[triele.n[i]].delta_At_old;
						}
					}
					//-------------------------------------------------------
					M += J[nodenumber];
					F += (Se[0] * mp_node[triele.n[0]].At_old + Se[1] * mp_node[triele.n[1]].At_old + Se[2] * mp_node[triele.n[2]].At_old - triele.J * triele.area / 3);
				}
				//----------------------------------------------
			}
			mp_node[n].delta_At = (NA - F) / M;
		}

		//-----------------------------------���Ե����������ж���������ж���ʽ���ǲ��ǿ�����һЩ�޸��أ�
		double errorJacobi = 0, a = 0, b = 0;
		for (int n = 0; n < m_num_nodes; ++n) {
			a += (mp_node[n].delta_At - mp_node[n].delta_At_old) * (mp_node[n].delta_At - mp_node[n].delta_At_old);
			b += mp_node[n].delta_At * mp_node[n].delta_At;
		}
		errorJacobi = sqrt(a) / sqrt(b);
		if ((iter + 1) % 100 == 0) {
			cout << "Jacobi Iteration step: " << iter + 1 << ", Relative error: " << errorJacobi << endl;
		}
		//------------------------------------------------

		//--------������Ե���������Ҫ�󣬽�At�滻��At + delta_At
		if (errorJacobi > maxerror) {
			for (int n = 0; n < m_num_nodes; ++n) {
				mp_node[n].delta_At_old = mp_node[n].delta_At;
			}
		}
		else {
			double errorNR = 0, aNR = 0, bNR = 0;
			for (int n = 0; n < m_num_nodes; ++n) {
				mp_node[n].At += mp_node[n].delta_At;
				mp_node[n].delta_At_old = 0;
				if (mp_node[n].x != 0) {
					mp_node[n].A = mp_node[n].At / mp_node[n].x;
				}
				aNR += mp_node[n].delta_At * mp_node[n].delta_At;
				bNR += mp_node[n].At * mp_node[n].At;
			}
			updateB();
			errorNR = sqrt(aNR) / sqrt(bNR);
			cout << "errorNR: " << errorNR << endl << endl;
			if (errorNR < maxerror) {
				return;
			}
			else
			{
				for (int n = 0; n < m_num_nodes; ++n) {
					mp_node[n].At_old = mp_node[n].At;
				}
			}
		}
		//------------------------------------------------------
	}
}

void FEM2DNDDRSolver::solve2DAximOpt1()
{
	for (int iter = 0; iter < maxitersteps; ++iter) {
#pragma omp parallel for num_threads(8)
		for (int n = 0; n < m_num_nodes; ++n) {
			if (mp_node[n].bdr == 1) {
				continue;
			}

			//�ڵ��ڲ���������
			int maxNRitersteps = 1;
			double Ati = 0, dAt = 0, NRerror = 1, count = 0;
			//double Ati = 0;
			double AtLocal[3]{ 0 ,0, 0 }, ALocal[3]{ 0, 0, 0 };
			for (int NRiter = 0; NRiter < maxNRitersteps; ++NRiter) {
				double S = 0, F = 0, NA = 0;
				double M = 0;
				double Se[3] = { 0, 0, 0 };
				for (int k = 0; k < mp_node[n].NumberofNeighbourElement; ++k) {
					int i_tri = mp_node[n].NeighbourElementId[k];	//ͨ����Ԫ�������ڵı�ż�����Ԫ��ȫ�ֱ��
					CTriElement triele = mp_triele[i_tri];
					int nodenumber = mp_node[n].NeighbourElementNumber[k];	//���������Ľڵ��Ӧ�ĵ�Ԫ�ֲ����
					double mu = triele.material->getMu(triele.B);
					double mut = mu * triele.xdot;
					//--------------------------------------���������ڵ�A
					for (int m = 0; m < 3; ++m) {
						if (m == nodenumber) {
							AtLocal[m] = Ati;
						}
						else {
							AtLocal[m] = mp_node[mp_triele[i_tri].n[m]].At_old;
						}
					}
					//-------------------------------------------------

					//--------------------------------------��Ԫ����
					if (triele.material->getLinearFlag() == true) {
						for (int i = 0; i < 3; ++i) {
							Se[i] = triele.C[nodenumber][i] / mut;
							if (i != nodenumber) {
								//NA -= Se[i] * mp_node[triele.n[i]].delta_At_old;
							}
						}
						M += Se[nodenumber];
						F += (Se[0] * AtLocal[0] + Se[1] * AtLocal[1] + Se[2] * AtLocal[2] - triele.J * triele.area / 3);
					}
					else
					{
						//-------------------------------------���㵥Ԫ��Jacobi����
						double dvdb, dvdbt, Bt, sigmai[3]{ 0, 0, 0 }, J[3]{ 0, 0, 0 };
						dvdb = triele.material->getdvdB(triele.B);
						dvdbt = dvdb / triele.xdot / triele.xdot;
						Bt = triele.B * triele.xdot;
						for (int i = 0; i < 3; ++i) {
							Se[i] = triele.C[nodenumber][i] / mut;
							for (int m = 0; m < 3; ++m) {
								sigmai[i] += triele.C[i][m] * AtLocal[m];
							}
							J[i] = triele.C[nodenumber][i] / mut;
							if (Bt != 0) {
								J[i] += dvdbt * sigmai[nodenumber] * sigmai[i] / Bt / triele.area;
							}
							if (i != nodenumber) {
								//NA -= J[i] * mp_node[triele.n[i]].delta_At_old;
							}
						}
						//-------------------------------------------------------
						M += J[nodenumber];
						F += (Se[0] * AtLocal[0] + Se[1] * AtLocal[1] + Se[2] * AtLocal[2] - triele.J * triele.area / 3);
					}
					//----------------------------------------------
				}
				dAt = (NA - F) / M;
				Ati = Ati + dAt;
				double a = dAt * dAt;
				double b = Ati * Ati;
				double NRerror = sqrt(a) / sqrt(b);
				if (Ati == 0) {
					break;
				}
				if (NRerror > 1e-6) {
					//cout << "n:" << n << ", NRiter: " << NRiter << ", deltaAt: " << mp_node[n].delta_At << ", Ati: " << Ati << endl;
					//cout << "n:" << n << ", NRerror: " << NRerror << endl;
					mp_node[n].delta_At = dAt;
					mp_node[n].At = Ati;
					mp_node[n].A = Ati / mp_node[n].x;
					for (int i = 0; i < mp_node[n].NumberofNeighbourElement; ++i) {
						updateB(mp_node[n].NeighbourElementId[i]);	//���updateB��λ���Ƿ��������һ�����޸��أ�
					}
				}
				else {
					//cout << "n:" << n << ", NRerror: " << NRerror << endl;
					break;
				}
			}
		}

		//-------------------------------------------ȫ���������ж�
		double error = 0, a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			/*a += mp_node[i].delta_At * mp_node[i].delta_At;*/
			a += (mp_node[i].At - mp_node[i].At_old) * (mp_node[i].At - mp_node[i].At_old);
			b += mp_node[i].At * mp_node[i].At;
		}
		error = sqrt(a) / sqrt(b);
		cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;

		//if ((iter + 1) % 100 == 0) {
		//	cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
		//}
		if (error > maxerror) {
			for (int i = 0; i < m_num_nodes; ++i) {
				mp_node[i].delta_At_old = mp_node[i].At - mp_node[i].At_old;
				mp_node[i].At_old = mp_node[i].At;
			}
		}
		else {
			staticsteps = iter;
			cout << "Iteration step: " << iter + 1 << endl;
			cout << "Nonlinear NDDR iteration finish.\n";
			return;
		}
		//------------------------------------------------
	}
}

void FEM2DNDDRSolver::solve2DAximPrecondition()
{
	//-----------1.NR������ȡ��һ���ֲڵĳ�ֵ-------------------------------------------
	//���ǵ�һ��߽�������װ��
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
						double mu = triele.material->getMu();	//����mu=0�����
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
				//���ż���
				double h_c = triele.material->getH_c();
				double theta_m = triele.material->getTheta_m();
				F[node_pos[n1]] += h_c / 2 * (triele.R[i] * cos(theta_m) - triele.Q[i] * sin(theta_m));
			}
		}
	}

	//�����Բ��ֵ���
	int pos1 = pos;
	std::vector<double> F1 = F;
	//std::vector<double> A_old(m_num_nodes, 0);
	for (int step = 0; step < maxitersteps; ++step) {
		cout << "Iteration step " << step + 1 << " start." << endl;
		pos = pos1;
		F = F1;
		for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
			CTriElement triele = mp_triele[i_tri];
			if (triele.material->getLinearFlag() == false) {
				//���㵥ԪJacobi����
				vector<vector<double>> J(3, vector<double>(3, 0));	//��ԪJacobi����
				vector<double> Fj(3, 0);	//�������Ҳ���
				double mu, mut, dvdb, dvdbt, Bt, sigmai, sigmaj;
				mu = triele.material->getMu(mp_triele[i_tri].B);
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
				//װ��
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
		//���
		locs[0].resize(pos);
		locs[1].resize(pos);
		cout << "pos: " << pos << endl;
		vector<double> res1 = matsolver->solveMatrix(locs, vals, F, pos, num_freenodes);
		for (int i = 0; i < num_freenodes; ++i) {
			int index = node_reorder[i];
			mp_node[index].At = res1[i];
			mp_node[index].At_old = res1[i];
			if (mp_node[index].x != 0) {
				mp_node[index].A = mp_node[index].At / mp_node[index].x;
			}
		}
		//���´ų����
		updateB();

		if (nonlinearelesize == 0) {
			cout << "Linear problem solved!\n";
			return;
		}
		//�ж�������
		double error = 0, a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (mp_node[i].A - mp_node[i].A_old) * (mp_node[i].A - mp_node[i].A_old);
			b += mp_node[i].A * mp_node[i].A;
		}
		error = sqrt(a) / sqrt(b);
		cout << "Relative error: " << error << endl;
		if (error > 1e-2) {
			for (int i = 0; i < m_num_nodes; ++i) {
				mp_node[i].A_old = mp_node[i].A;
			}
		}
		else {
			staticsteps = step + 1;
			cout << "Nonlinear iteration finish.\n";
			//return;
			break;
		}
	}
	//staticsteps = maxitersteps;
	//cout << "Warning: Number of iterations out of limit.\n";

	//-----------2.NDDR��������ȷ����-------------------------------------------
	solve2DAxim1();
//	for (int iter = 0; iter < maxitersteps; ++iter) {
//#pragma omp parallel for num_threads(8)
//		for (int n = 0; n < m_num_nodes; ++n) {
//			if (mp_node[n].bdr == 1) {
//				continue;
//			}
//
//			//�ڵ��ڲ���������
//			int maxNRitersteps = 10;
//			double Ati = 0, dAt = 0, NRerror = 1, count = 0;
//			//double Ati = 0;
//			double AtLocal[3]{ 0 ,0, 0 }, ALocal[3]{ 0, 0, 0 };
//			for (int NRiter = 0; NRiter < maxNRitersteps; ++NRiter) {
//				double S = 0, F = 0, NA = 0;
//				double M = 0;
//				double Se[3] = { 0, 0, 0 };
//				for (int k = 0; k < mp_node[n].NumberofNeighbourElement; ++k) {
//					int i_tri = mp_node[n].NeighbourElementId[k];	//ͨ����Ԫ�������ڵı�ż�����Ԫ��ȫ�ֱ��
//					CTriElement triele = mp_triele[i_tri];
//					int nodenumber = mp_node[n].NeighbourElementNumber[k];	//���������Ľڵ��Ӧ�ĵ�Ԫ�ֲ����
//					double mu = triele.material->getMu(triele.B);
//					double mut = mu * triele.xdot;
//					//--------------------------------------���������ڵ�A
//					for (int m = 0; m < 3; ++m) {
//						if (m == nodenumber) {
//							AtLocal[m] = Ati;
//						}
//						else {
//							AtLocal[m] = mp_node[mp_triele[i_tri].n[m]].At_old;
//						}
//					}
//					//-------------------------------------------------
//
//					//--------------------------------------��Ԫ����
//					if (triele.material->getLinearFlag() == true) {
//						for (int i = 0; i < 3; ++i) {
//							Se[i] = triele.C[nodenumber][i] / mut;
//							if (i != nodenumber) {
//								//NA -= Se[i] * mp_node[triele.n[i]].delta_At_old;
//							}
//						}
//						M += Se[nodenumber];
//						F += (Se[0] * AtLocal[0] + Se[1] * AtLocal[1] + Se[2] * AtLocal[2] - triele.J * triele.area / 3);
//					}
//					else
//					{
//						//-------------------------------------���㵥Ԫ��Jacobi����
//						double dvdb, dvdbt, Bt, sigmai[3]{ 0, 0, 0 }, J[3]{ 0, 0, 0 };
//						dvdb = triele.material->getdvdB(triele.B);
//						dvdbt = dvdb / triele.xdot / triele.xdot;
//						Bt = triele.B * triele.xdot;
//						for (int i = 0; i < 3; ++i) {
//							Se[i] = triele.C[nodenumber][i] / mut;
//							for (int m = 0; m < 3; ++m) {
//								sigmai[i] += triele.C[i][m] * AtLocal[m];
//							}
//							J[i] = triele.C[nodenumber][i] / mut;
//							if (Bt != 0) {
//								J[i] += dvdbt * sigmai[nodenumber] * sigmai[i] / Bt / triele.area;
//							}
//							if (i != nodenumber) {
//								//NA -= J[i] * mp_node[triele.n[i]].delta_At_old;
//							}
//						}
//						//-------------------------------------------------------
//						M += J[nodenumber];
//						F += (Se[0] * AtLocal[0] + Se[1] * AtLocal[1] + Se[2] * AtLocal[2] - triele.J * triele.area / 3);
//					}
//					//----------------------------------------------
//				}
//				dAt = (NA - F) / M;
//				Ati = Ati + dAt;
//				double a = dAt * dAt;
//				double b = Ati * Ati;
//				double NRerror = sqrt(a) / sqrt(b);
//				if (Ati == 0) {
//					break;
//				}
//				if (NRerror > 1e-6) {
//					//cout << "n:" << n << ", NRiter: " << NRiter << ", deltaAt: " << mp_node[n].delta_At << ", Ati: " << Ati << endl;
//					//cout << "n:" << n << ", NRerror: " << NRerror << endl;
//					mp_node[n].delta_At = dAt;
//					mp_node[n].At = Ati;
//					mp_node[n].A = Ati / mp_node[n].x;
//					for (int i = 0; i < mp_node[n].NumberofNeighbourElement; ++i) {
//						updateB(mp_node[n].NeighbourElementId[i]);	//���updateB��λ���Ƿ��������һ�����޸��أ�
//					}
//				}
//				else {
//					//cout << "n:" << n << ", NRerror: " << NRerror << endl;
//					break;
//				}
//			}
//		}
//
//		//-------------------------------------------ȫ���������ж�
//		double error = 0, a = 0, b = 0;
//		for (int i = 0; i < m_num_nodes; ++i) {
//			/*a += mp_node[i].delta_At * mp_node[i].delta_At;*/
//			a += (mp_node[i].At - mp_node[i].At_old) * (mp_node[i].At - mp_node[i].At_old);
//			b += mp_node[i].At * mp_node[i].At;
//		}
//		error = sqrt(a) / sqrt(b);
//		cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
//
//		//if ((iter + 1) % 100 == 0) {
//		//	cout << "Iteration step: " << iter + 1 << ", Relative error: " << error << endl;
//		//}
//		if (error > maxerror) {
//			for (int i = 0; i < m_num_nodes; ++i) {
//				mp_node[i].delta_At_old = mp_node[i].At - mp_node[i].At_old;
//				mp_node[i].At_old = mp_node[i].At;
//			}
//		}
//		else {
//			staticsteps = iter;
//			cout << "Iteration step: " << iter + 1 << endl;
//			cout << "Nonlinear NDDR iteration finish.\n";
//			return;
//		}
//		//------------------------------------------------
//	}
}

void FEM2DNDDRSolver::solve2DAximRobin()
{
	DataPrepare();
	JsSumCalculate();
	SumNeiborJsSumCalculate();

	for (int iter = 0; iter < maxitersteps; ++iter) {
		ElmRHSContriCalculate();
		SumNodeRHSCalculate();
		UpdateSolutiontoA1();
		CopyA1toA0();
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
			//�ڵ��ڲ���������
			int maxNRitersteps = 6;
			double Ai = 0;
			for (int NRiter = 0; NRiter < maxNRitersteps; ++NRiter) {
				//cout << "n: " << n <<", NRiter: " << NRiter << endl;
				double S = 0, F = 0;
				double J = 0, Fj = 0;
				//װ�����
				for (int k = 0; k < mp_node[n].NumberofNeighbourElement; ++k) {
					//vector<vector<double>> J(3, vector<double>(3, 0));	//��ԪJacobi����
					//vector<double> Fj(3, 0);	//�������Ҳ���
					int i_tri = mp_node[n].NeighbourElementId[k];
					CTriElement triele = mp_triele[i_tri];
					int nodenumber = mp_node[n].NeighbourElementNumber[k];
					double mu = triele.material->getMu(triele.B);
					//printf("mu: %f\n", mut);
					//�������Ե�Ԫ
					if (triele.material->getLinearFlag() == true) {
						for (int i = 0; i < 3; ++i) {
							double Se = triele.C[nodenumber][i] / mu;
							if (nodenumber == i) {
								S += Se;
								F += triele.J * triele.area / 3;
								////���Ų���
								//double h_c = triele.material->getH_c();
								//double theta_m = triele.material->getTheta_m();
								//F += h_c / 2 * (triele.R[i] * cos(theta_m) - triele.Q[i] * sin(theta_m));
							}
							else {
								F -= Se * A_old[triele.n[i]];
							}
						}
					}
					//��������Ե�Ԫ
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

				//NR�����������ж� 
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

		//�ж�ȫ��������
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
			int NRCount = 100;
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
						
						//���Ե�Ԫ�������Ҳ�������ʱ��ҲҪ����A������
						if (mp_triele[ID].material->getLinearFlag() == true)	//���Կ�����Ԫ
						{
							V = 1.0 / mp_triele[ID].material->getMu(0);
							J0 = V * mp_triele[ID].C[nodenumber][nodenumber];
							k0 = Js - V * RHSContri;

						}

						else	//�����Ե�Ԫ
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (ALocal[0] - ALocal[1]) * (ALocal[0] - ALocal[1]) + mp_triele[ID].C[1][2] * (ALocal[1] - ALocal[2]) * (ALocal[1] - ALocal[2]) + mp_triele[ID].C[0][2] * (ALocal[0] - ALocal[2]) * (ALocal[0] - ALocal[2]));
							B = sqrt(B2);	//����B�����ڴ��������
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
						count++;
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
				}

				mp_node[i].A = A;
			}
		}

		//�ж�ȫ��������
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
			//��һ��Ӱ��Ч��
			for (int i = 0; i < m_num_nodes; ++i) {
				mp_node[i].A_old = mp_node[i].A;
			}
		}
		else {
			cout << "Iteration step: " << iter + 1 << endl;
			cout << "Nonlinear NDDR iteration finish.\n";
			break;
		}
	}

	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		double bx = 0, by = 0;
		for (int i = 0; i < 3; ++i) {
			int n = mp_triele[i_tri].n[i];
			bx += mp_triele[i_tri].R[i] * mp_node[n].A;
			by += mp_triele[i_tri].Q[i] * mp_node[n].A;
		}
		bx = bx / 2 / mp_triele[i_tri].area;
		mp_triele[i_tri].Bx = bx;
		by = -by / 2 / mp_triele[i_tri].area;
		mp_triele[i_tri].By = by;
		mp_triele[i_tri].B = sqrt(bx * bx + by * by);
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
					//���Ե�Ԫ�������Ҳ�������ʱ��ҲҪ����A������
					if (mp_triele[ID].material->getLinearFlag() == true)	//���Կ�����Ԫ
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
					else	//�����Ե�Ԫ
					{
						double Mu0 = 0.002;
						if (mp_node[i].NeighbourElementNumber[j] == 0)
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (A - A2) * (A - A2) + mp_triele[ID].C[1][2] * (A2 - A3) * (A2 - A3) + mp_triele[ID].C[0][2] * (A - A3) * (A - A3));
							B = sqrt(B2);	//����B�����ڴ��������
							V = 1 / mp_triele[ID].material->getMu(B);
							VB2 = mp_triele[ID].material->getdvdB2(B);
							if (B < 0.6)
							{
								//V = 1 / Mu0;
								//VB2 = 0;
								J0 = V * mp_triele[ID].C[0][0];
								k0 = Js - A * V * mp_triele[ID].C[0][0] - A2 * V * mp_triele[ID].C[0][1] - A3 * V * mp_triele[ID].C[0][2];

							}
							else
							{
								//V = 1 / Mu0 + 3000 * (B - 0.6) * (B - 0.6)* (B - 0.6) / B;
								//VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6)* (B - 0.6)) / B / B / 2 / B;
								B2A = -1 / mp_triele[ID].area * (2 * mp_triele[ID].C[0][1] * (A - A2) + 2 * mp_triele[ID].C[0][2] * (A - A3));//dfrac{dB^2}{dA}?
								J0 = B2A * VB2 * (mp_triele[ID].C[0][0] * A + mp_triele[ID].C[0][1] * A2 + mp_triele[ID].C[0][2] * A3) + V * mp_triele[ID].C[0][0];
								k0 = Js - V * (mp_triele[ID].C[0][0] * A + mp_triele[ID].C[0][1] * A2 + mp_triele[ID].C[0][2] * A3);
							}

						}
						else if (mp_node[i].NeighbourElementNumber[j] == 1)
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (A1 - A) * (A1 - A) + mp_triele[ID].C[1][2] * (A - A3) * (A - A3) + mp_triele[ID].C[0][2] * (A1 - A3) * (A1 - A3));
							B = sqrt(B2);
							V = 1 / mp_triele[ID].material->getMu(B);
							VB2 = mp_triele[ID].material->getdvdB2(B);
							if (B < 0.6)
							{
								//V = 1 / Mu0;
								//VB2 = 0;
								J0 = V * mp_triele[ID].C[1][1];
								k0 = Js - A1 * V * mp_triele[ID].C[1][0] - A * V * mp_triele[ID].C[1][1] - A3 * V * mp_triele[ID].C[1][2];

							}
							else
							{
								//V = 1 / Mu0 + 3000 * (B - 0.6) * (B - 0.6)* (B - 0.6) / B;
								//VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6)* (B - 0.6)) / B / B / 2 / B;
								B2A = -1 / mp_triele[ID].area * (2 * mp_triele[ID].C[0][1] * (A - A1) + 2 * mp_triele[ID].C[1][2] * (A - A3));
								J0 = B2A * VB2 * (mp_triele[ID].C[1][0] * A1 + mp_triele[ID].C[1][1] * A + mp_triele[ID].C[1][2] * A3) + V * mp_triele[ID].C[1][1];
								k0 = Js - V * (mp_triele[ID].C[1][0] * A1 + mp_triele[ID].C[1][1] * A + mp_triele[ID].C[1][2] * A3);
							}
						}

						else if (mp_node[i].NeighbourElementNumber[j] == 2)
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (A1 - A2) * (A1 - A2) + mp_triele[ID].C[1][2] * (A - A2) * (A - A2) + mp_triele[ID].C[0][2] * (A - A1) * (A - A1));
							B = sqrt(B2);
							V = 1 / mp_triele[ID].material->getMu(B);
							VB2 = mp_triele[ID].material->getdvdB2(B);
							if (B < 0.6)
							{
								//V = 1 / Mu0;
								//VB2 = 0;
								J0 = V * mp_triele[ID].C[2][2];
								k0 = Js - A1 * V * mp_triele[ID].C[2][0] - A2 * V * mp_triele[ID].C[2][1] - A * V * mp_triele[ID].C[2][2];

							}
							else
							{
								//V = 1 / Mu0 + 3000 * (B - 0.6) * (B - 0.6)* (B - 0.6) / B;
								//VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6)* (B - 0.6)) / B / B / 2 / B;
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
		double B2, B, V, VB2, B2A, Vtest, VB2test;
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
					//���Ե�Ԫ�������Ҳ�������ʱ��ҲҪ����A������
					if (mp_triele[ID].material->getLinearFlag() == true)	//���Կ�����Ԫ
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
					else	//�����Ե�Ԫ
					{
						double Mu0 = 0.002;
						if (mp_node[i].NeighbourElementNumber[j] == 0)
						{
							B2 = -1 / mp_triele[ID].area * (mp_triele[ID].C[0][1] * (A - A2) * (A - A2) + mp_triele[ID].C[1][2] * (A2 - A3) * (A2 - A3) + mp_triele[ID].C[0][2] * (A - A3) * (A - A3));
							B = sqrt(B2);	//����B�����ڴ��������
							V = 1 / mp_triele[ID].material->getMu(B);
							VB2 = mp_triele[ID].material->getdvdB2(B);
							if (B < 0.6)
							{
								//V = 1 / Mu0;
								//VB2 = 0;
								J0 = V * mp_triele[ID].C[0][0];
								k0 = Js - A * V * mp_triele[ID].C[0][0] - A2 * V * mp_triele[ID].C[0][1] - A3 * V * mp_triele[ID].C[0][2];

							}
							else
							{
								//V = 1 / Mu0 + 3000 * (B - 0.6) * (B - 0.6) * (B - 0.6) / B;
								//VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6) * (B - 0.6)) / B / B / 2 / B;
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
							V = 1 / mp_triele[ID].material->getMu(B);
							VB2 = mp_triele[ID].material->getdvdB2(B);
							if (B < 0.6)
							{
								//V = 1 / Mu0;
								//VB2 = 0;
								J0 = V * mp_triele[ID].C[1][1];
								k0 = Js - A1 * V * mp_triele[ID].C[1][0] - A * V * mp_triele[ID].C[1][1] - A3 * V * mp_triele[ID].C[1][2];

							}
							else
							{
								//V = 1 / Mu0 + 3000 * (B - 0.6) * (B - 0.6) * (B - 0.6) / B;
								//VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6) * (B - 0.6)) / B / B / 2 / B;
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
							V = 1 / mp_triele[ID].material->getMu(B);
							VB2 = mp_triele[ID].material->getdvdB2(B);
							if (B < 0.6)
							{
								//V = 1 / Mu0;
								//VB2 = 0;
								J0 = V * mp_triele[ID].C[2][2];
								k0 = Js - A1 * V * mp_triele[ID].C[2][0] - A2 * V * mp_triele[ID].C[2][1] - A * V * mp_triele[ID].C[2][2];

							}
							else
							{
								//V = 1 / Mu0 + 3000 * (B - 0.6) * (B - 0.6) * (B - 0.6) / B;
								//VB2 = (B * 9000.0 * (B - 0.6) * (B - 0.6) - 3000.0 * (B - 0.6) * (B - 0.6) * (B - 0.6)) / B / B / 2 / B;
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

void FEM2DNDDRSolver::DataPrepare()
{
	int* dummy = new int[m_num_nodes]();

	for (int i = 0; i < m_num_nodes; ++i) {
		for (int j = 0; j < m_num_nodes; ++j) {
			dummy[j] = 0;
		}
		for (int j = 0; j < mp_node[i].NumberofNeighbourElement; ++j) {
			for (int k = 0; k < 3; ++k) {
				dummy[mp_triele[mp_node[i].NeighbourElementId[j]].n[k]] = 1;
			}
		}
		for (int j = 0; j < m_num_nodes; ++j) {
			if (dummy[j] == 1 && j != i) {
				mp_node[i].NeiborNode[mp_node[i].NumNeiborNodes] = j;
				mp_node[i].NumNeiborNodes++;
			}
		}
	}

	delete[] dummy;

	for (int e = 0; e < m_num_triele; ++e) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				mp_triele[e].ElmRowSum[j][k] = 0;
			}
		}

		double temp;
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				if (mp_node[mp_triele[e].n[k]].bdr != 1) {
					if (k == j)
						temp = 1;
					else
						temp = 1 / Gamma;
					mp_triele[e].ElmRowSum[j][0] += temp * mp_triele[e].C[k][0];
					mp_triele[e].ElmRowSum[j][1] += temp * mp_triele[e].C[k][1];
					mp_triele[e].ElmRowSum[j][2] += temp * mp_triele[e].C[k][2];
				}
			}
		}
	}
}

void FEM2DNDDRSolver::JsSumCalculate()
{
	for (int node = 0; node < m_num_nodes; ++node) {
		mp_node[node].JsSum = 0;
		for (int e = 0; e < mp_node[node].NumberofNeighbourElement; ++e) {
			mp_node[node].JsSum += (mp_triele[mp_node[node].NeighbourElementId[e]].J) * mp_triele[mp_node[node].NeighbourElementId[e]].area / 3;
		}
	}
}

void FEM2DNDDRSolver::SumNeiborJsSumCalculate()
{
	for (int node = 0; node < m_num_nodes; ++node) {
		for (int n = 0; n < mp_node[node].NumNeiborNodes; ++n) {
			mp_node[node].SumNeiborJsSum += mp_node[mp_node[node].NeiborNode[n]].JsSum;
		}
	}
}

void FEM2DNDDRSolver::ElmRHSContriCalculate()
{
	for (int e = 0; e < m_num_triele; ++e) {
		int row, col;
		for (row = 0; row < 3; ++row) {
			mp_triele[e].RHSContri[row] = 0;
		}
		for (row = 0; row < 3; ++row) {
			for (col = 0; col < 3; ++col) {
				mp_triele[e].RHSContri[row] += mp_triele[e].C[row][col] * mp_node[mp_triele[e].n[col]].At_old;
			}
			mp_triele[e].RHSContri[row] = mp_triele[e].RHSContri[row] / mp_triele[e].mut;
		}
	}
}

void FEM2DNDDRSolver::SumNodeRHSCalculate()
{
	for (int node = 0; node < m_num_nodes; ++node) {
		int e;
		mp_node[node].SumRHSContri = 0;
		for (e = 0; e < mp_node[node].NumberofNeighbourElement; e++) {
			mp_node[node].SumRHSContri += mp_triele[mp_node[node].NeighbourElementId[e]].RHSContri[mp_node[node].NeighbourElementNumber[e]];
		}
	}
}

void FEM2DNDDRSolver::UpdateSolutiontoA1()
{
#pragma omp parallel for num_threads(8)
	for (int n = 0; n < m_num_nodes; ++n) {
		double RHS, LHS, dF_dA, NRsolution, temp;
		int i, j, NeiborID, e, LocalPos;
		double ALocal[3], B2, B, Bt, V, dvdb, dvdbt, dbda, sigma = 0;
		int NRct;


		if (mp_node[n].bdr != 1) {
			//-------------Get RHS
			RHS = 0;

			//�߽��ϵ�ȫ���ڵ��SumRHSContri֮��
			for (i = 0; i < mp_node[n].NumNeiborNodes; ++i) {
				NeiborID = mp_node[n].NeiborNode[i];
				if (mp_node[NeiborID].bdr != 1) {
					RHS -= mp_node[NeiborID].SumRHSContri;
				}
			}

			//������ȫ����Ԫ��Ϊ�߽繱�׵�RHSContri֮��
			for (i = 0; i < mp_node[n].NumberofNeighbourElement; ++i) {
				e = mp_node[n].NeighbourElementId[i];
				for (j = 0; j < 3; j++) {
					if (n != mp_triele[e].n[j] && mp_node[mp_triele[e].n[j]].bdr != 1) {
						RHS += mp_triele[e].RHSContri[j];
					}
				}
			}

			RHS = (RHS + mp_node[n].SumNeiborJsSum) / Gamma + mp_node[n].JsSum;
			//---------------Get LHS and dF_dA
			NRsolution = mp_node[n].At_old;
			//---------------NRiteration
			for (NRct = 0; NRct < 1; ++NRct) {
				LHS = 0;
				dF_dA = 0;
				for (i = 0; i < mp_node[n].NumberofNeighbourElement; ++i) {
					e = mp_node[n].NeighbourElementId[i];
					LocalPos = mp_node[n].NeighbourElementNumber[i];
					ALocal[0] = mp_node[mp_triele[e].n[0]].At_old;
					ALocal[1] = mp_node[mp_triele[e].n[1]].At_old;
					ALocal[2] = mp_node[mp_triele[e].n[2]].At_old;
					ALocal[LocalPos] = NRsolution;

					//��Ҫ����B�ļ���
					//updateB(e);
					//B = mp_triele[e].B;

					//����B
					double bx = 0, by = 0;
					for (int bi = 0; bi < 3; ++bi) {
						int n = mp_triele[e].n[bi];
						bx += mp_triele[e].R[bi] * ALocal[bi];
						by += mp_triele[e].Q[bi] * ALocal[bi];
					}
					bx = bx / 2 / mp_triele[e].area / mp_triele[e].xdot;
					mp_triele[e].Bx = bx;
					by = -by / 2 / mp_triele[e].area / mp_triele[e].xdot;
					mp_triele[e].By = by;
					mp_triele[e].B = sqrt(bx * bx + by * by);
					B = mp_triele[e].B;
					Bt = B * mp_triele[e].xdot;
					dvdb = mp_triele[e].material->getdvdB(mp_triele[e].B);
					dvdbt = dvdb / mp_triele[e].xdot / mp_triele[e].xdot;
					//cout << "n: " << n << ", i: " << i << ", Bt: " << Bt << ", dvdb: " << dvdb << ", dvdbt:" << dvdbt << endl;

					for (int m = 0; m < 3; ++m) {
						sigma += mp_triele[e].C[LocalPos][m] * ALocal[m];
					}
					if (Bt != 0) {
						dbda = sigma / Bt / mp_triele[e].area;
					}

					temp = 0;
					for (j = 0; j < 3; j++) {
						dF_dA += dvdb * dbda * ALocal[j] * mp_triele[e].ElmRowSum[LocalPos][j];
						temp += ALocal[j] * mp_triele[e].ElmRowSum[LocalPos][j];
					}
					dF_dA += mp_triele[e].ElmRowSum[LocalPos][LocalPos] / (mp_triele[e].material->getMu(mp_triele[e].B) * mp_triele[e].xdot);
					LHS += temp / (mp_triele[e].material->getMu(mp_triele[e].B) * mp_triele[e].xdot);

					////�����������
					//for (j = 0; j < 3; ++j) {
					//	LHS += ALocal[j] * mp_triele[e].ElmRowSum[LocalPos][j] / (4 * PI * 1e-3);
					//}
					//dF_dA += mp_triele[e].ElmRowSum[LocalPos][LocalPos] / (4 * PI * 1e-3);
				}

				NRsolution += (RHS - LHS) / dF_dA;
			}
			mp_node[n].At = NRsolution;
		}
	}
}

void FEM2DNDDRSolver::CopyA1toA0()
{
	double a, b, error;
	for (int i = 0; i < m_num_nodes; ++i) {
		/*a += mp_node[i].delta_At * mp_node[i].delta_At;*/
		a += (mp_node[i].At - mp_node[i].At_old) * (mp_node[i].At - mp_node[i].At_old);
		b += mp_node[i].At * mp_node[i].At;
	}
	error = sqrt(a) / sqrt(b);
	cout << "Relative error: " << error << endl;

	for (int i = 0; i < m_num_nodes; ++i) {

		mp_node[i].At_old = mp_node[i].At;
		//mp_node[i].A = mp_node[i].At / mp_node[i].x;
	}
}

void FEM2DNDDRSolver::UpdateVe()
{

}
