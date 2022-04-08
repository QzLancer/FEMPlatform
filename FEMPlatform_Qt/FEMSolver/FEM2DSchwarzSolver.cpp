#include "FEM2DSchwarzSolver.h"
#include "MatrixOutput.h"

FEM2DSchwarzSolver::FEM2DSchwarzSolver(int _numofdomain):
	d_node_pos(NULL),
	d_node_reorder(NULL),
	d_nodeid(NULL),
	nparttable(NULL),
	crossbdrtable(NULL),
	eparttable(NULL),
	d_eleid(NULL),
	d_elesize(NULL),
	d_dof(NULL),
	d_nodesize(NULL),
	metisfile()
{
	numofdomain = _numofdomain;
	
}

FEM2DSchwarzSolver::~FEM2DSchwarzSolver()
{
	for (int i = 0; i < numofdomain; ++i) {
		if (d_eleid[i] != NULL) delete[] d_eleid[i];
		if (eparttable[i] != NULL) delete[] eparttable[i];
		if (crossbdrtable[i] != NULL) delete[] crossbdrtable[i];
		if (nparttable[i] != NULL) delete[] nparttable[i];
		if (d_nodeid[i] != NULL) delete[] d_nodeid[i];
		if (d_node_reorder[i] != NULL) delete[] d_node_reorder[i];
		if (d_node_pos[i] != NULL) delete[] d_node_pos[i];
	}
	if (d_node_pos != NULL) delete[] d_node_pos;
	if (d_node_reorder != NULL) delete[] d_node_reorder;
	if (d_nodeid != NULL) delete[] d_nodeid;
	if (nparttable != NULL) delete[] nparttable;
	if (d_nodesize != NULL) delete[] d_nodesize;
	if (d_dof != NULL) delete[] d_dof;
	if (d_elesize != NULL) delete[] d_elesize;
	if (d_eleid != NULL) delete[] d_eleid;
	if (eparttable != NULL) delete[] eparttable;
	if (crossbdrtable != NULL) delete[] crossbdrtable;
}

void FEM2DSchwarzSolver::solveStatic()
{
	generateMetisMesh();
	makeTrangle();
	processBoundaryCondition();
	processMaterial();
	processLoad();
	if (dimension == FEMModel::DIMENSION::D2AXISM) {
		solve2DAxim();
	}
	else {
		cout << "Error: can't solve this dimension" << endl;
	}
}

void FEM2DSchwarzSolver::solveDynamic()
{
	generateMetisMesh();
}

void FEM2DSchwarzSolver::processBoundaryCondition()
{
	//schwarz����ֽ��µı߽���������
	//���������ֱ߽磺1��ԭ��ģ���еĵ�һ��߽�������2�����򽻽紦�ĵ�һ��߽�����
	//����ȫ����Ȼ�߽�����
	num_freenodes = m_num_nodes;
	std::vector<int> boundarynodes;
	//ԭ�а汾�ı߽����
	for (auto iter = boundarymap.begin(); iter != boundarymap.end(); ++iter) {
		if (iter->second->getBoundaryType() == 1) {
			int boundarynodedomain = iter->first;
			for (int i = 0; i < m_num_edgele; ++i) {
				if (mp_edgele[i].domain == boundarynodedomain) {
					boundarynodes.push_back(mp_edgele[i].n[0]);
					boundarynodes.push_back(mp_edgele[i].n[1]);
				}
			}
		}
		else {
			cout << "Invalid boundarytype!\n";
			exit(0);
		}
	}

	//��x=0��ȫ���ڵ���Ϊ�߽�ڵ�
	for (int i = 0; i < m_num_nodes; ++i) {
		if (mp_node[i].x == 0) {
			boundarynodes.push_back(i);
		}
	}

	std::sort(boundarynodes.begin(), boundarynodes.end());
	boundarynodes.erase(unique(boundarynodes.begin(), boundarynodes.end()), boundarynodes.end());

	////����߽�㣬����ͨ��matlab��ͼ�ж��Ƿ������ȫ���ı߽�
	//for (auto a :boundarynodes) {
	//	cout << mp_node[a].x << " " << mp_node[a].y << endl;
	//}

	num_freenodes -= boundarynodes.size();
	for (auto a : boundarynodes) {
		mp_node[a].bdr = 1;
	}

	cout << "num_bdr: " << boundarynodes.size() << ", num_freenodes: " << num_freenodes << endl;

	//����ÿ����������ɽڵ���Ŀ���������ڵ����ɽڵ�ͱ߽�ڵ��������
	d_dof = new int[numofdomain];
	d_node_reorder = new int* [numofdomain];
	d_node_pos = new int* [numofdomain];
	for (int domain = 0; domain < numofdomain; ++domain) {
		d_node_reorder[domain] = new int[d_nodesize[domain]]();
		d_node_pos[domain] = new int[d_nodesize[domain]]();
		d_dof[domain] = d_nodesize[domain];
		int node_bdr = 0;
		for (int n = 0; n < d_nodesize[domain]; ++n) {
			int globaln = d_nodeid[domain][n];
			if (mp_node[globaln].bdr == 1 || crossbdrtable[domain][globaln] == 1) {
				d_dof[domain]--;
				++node_bdr;
				d_node_reorder[domain][d_nodesize[domain] - node_bdr] = n;
				d_node_pos[domain][n] = d_nodesize[domain] - node_bdr;
			}
			else {
				d_node_reorder[domain][n - node_bdr] = n;
				d_node_pos[domain][n] = n - node_bdr;
			}
		}
		printf("domain:%d, nodesize:%d, dof:%d\n", domain, d_nodesize[domain], d_dof[domain]);
	}
	printf("***********************************************************\n");

}

void FEM2DSchwarzSolver::generateMetisMesh()
{
	//����metis����������ļ�
	sprintf(metisfile, "../model/metis/%s_%d.metis", modelname.c_str(), m_num_triele);
	FILE* fp = fopen(metisfile, "w+");
	if (!fp) {
		std::cout << "FEM2DSchwarzSolver::FEM2DSchwarzSolver Error: Can't write metisfile.\n";
		return;
	}
	//д�뵥Ԫ��Ŀ�ͽڵ��ţ���1��ʼ��,�˴�ֻ���������ε�Ԫ���Ƿ���ԶԱ߽絥Ԫͬʱ���в�����
	fprintf(fp, "%d\n", m_num_triele);
	for (int i = 0; i < m_num_triele; ++i) {
		fprintf(fp, "%d %d %d\n", mp_triele[i].n[0] + 1, mp_triele[i].n[1] + 1, mp_triele[i].n[2] + 1);
	}
	fclose(fp);
	//����metis API���з���
	char part[3];
	sprintf(part, "%d", numofdomain);
	char a[] = "Administrator";
	char* myargv[] = { a, metisfile, part };
	mpmetis(3, myargv);
	//��ȡ�ڵ�����ֽ��ļ�
	char npartname[256];
	sprintf(npartname, "%s.npart.%d", metisfile, numofdomain);
	fp = fopen(npartname, "r");
	if (!fp) {
		std::cout << "FEM2DSchwarzSolver::FEM2DSchwarzSolver Error: Can't read npart file.\n";
		return;
	}
	int* ndomain = new int[m_num_nodes];
	int b;
	for (int i = 0; i < m_num_nodes; ++i) {
		fscanf(fp, "%d\n", &ndomain[i]);
	}
	//��ȡ�ڵ��ܱߵ�ȫ����Ԫ
	for (int e = 0; e < m_num_triele; ++e) {
		for (int i = 0; i < 3; ++i) {
			int nodeid = mp_triele[e].n[i];
			mp_node[nodeid].NeighbourElementId[mp_node[nodeid].NumberofNeighbourElement] = e;
			mp_node[nodeid].NumberofNeighbourElement += 1;
		}
		//printf("e:%d, %d %d %d\n", e, mp_triele[e].n[0], mp_triele[e].n[1], mp_triele[e].n[2]);
	}
	//�����ڵ����ڵ����򣬽������ڰ����ĵ�Ԫ�ֲ���ű�����eparttable��
	eparttable = new int* [numofdomain];
	for (int domain = 0; domain < numofdomain; ++domain) {
		eparttable[domain] = new int[m_num_triele]();
		for (int e = 0; e < m_num_triele; ++e) {
			eparttable[domain][e] = -1;
		}
	}
	for (int n = 0; n < m_num_nodes; ++n) {
		CNode node = mp_node[n];
		//cout << "n: " << n + 1 << ", numofneibourele: " << node.NumberofNeighbourElement << endl;
		int domain = ndomain[n];
		for (int e = 0; e < node.NumberofNeighbourElement; ++e) {
			int eleid = node.NeighbourElementId[e];
			eparttable[domain][eleid] = 0;
		}
	}

	//���������ڵĵ�Ԫ����������epartTable���±�ţ�ȫ�ֵ�Ԫ���ӳ�䵽����Ԫ��ţ�
	d_elesize = new int[numofdomain]();
	for (int domain = 0; domain < numofdomain; ++domain) {
		for (int e = 0; e < m_num_triele; ++e) {
			if (eparttable[domain][e] != -1) {
				eparttable[domain][e] = d_elesize[domain];
				d_elesize[domain]++;
			}
		}
		//cout << "domainelesize: " << domainelesize[domain] << endl;
	}
	printintMatrix("eparttable.csv", numofdomain, m_num_triele, eparttable);

	//����domaineleid������Ԫ���ӳ�䵽ȫ�ֵ�Ԫ��ţ�
 	d_eleid = new int*[numofdomain];
	for (int domain = 0; domain < numofdomain; ++domain) {
		d_eleid[domain] = new int[d_elesize[domain]];
		for (int e = 0; e < m_num_triele; ++e) {
			if (eparttable[domain][e] != -1) {
				d_eleid[domain][eparttable[domain][e]] = e;
			}
		}
	}

	//����ÿ����Ԫ�ı߽�ڵ㣬������boundaryTable��
	crossbdrtable = new int* [numofdomain];
	for (int domain = 0; domain < numofdomain; ++domain) {
		crossbdrtable[domain] = new int[m_num_nodes]();
	}
	for (int e = 0; e < m_num_triele; ++e) {
		for (int i = 0; i < 3; ++i) {
			int nodeid = mp_triele[e].n[i];
			int nodedomain = ndomain[nodeid];
			for (int eledomain = 0; eledomain < numofdomain; ++eledomain) {
				if (eparttable[eledomain][e] != -1 && nodedomain != eledomain) {
					crossbdrtable[eledomain][nodeid] = 1;
					//cout << "e: " << e << ", n: " << nodeid << ", nodedomain: " << nodedomain <<endl;
				}
			}
		}
	}
	//printintMatrix("crossbdrtable.csv", numofdomain, m_num_nodes, crossbdrtable);

	//�����������Ľڵ������������ֲ���ŵ�ȫ�ֱ�ŵ�ӳ��
	d_nodesize = new int[numofdomain]();
	d_nodeid = new int*[numofdomain];
	nparttable = new int* [numofdomain];
	for (int domain = 0; domain < numofdomain; ++domain) {
		nparttable[domain] = new int[m_num_nodes];
		for (int n = 0; n < m_num_nodes; ++n) {
			if (ndomain[n] == domain || crossbdrtable[domain][n] == 1) {
				d_nodesize[domain]++;
			}
			nparttable[domain][n] = -1;
		}
		d_nodeid[domain] = new int[d_nodesize[domain]];
		int pos = 0;
		for (int n = 0; n < m_num_nodes; ++n) {
			if (ndomain[n] == domain || crossbdrtable[domain][n] == 1) {
				d_nodeid[domain][pos] = n;
				nparttable[domain][n] = pos;
				pos++;
			}
		}
	}
	//printintMatrix("nparttable.csv", numofdomain, m_num_nodes, nparttable);
	delete[] ndomain;
}

void FEM2DSchwarzSolver::solve2DAxim()
{
	//ÿ��������з���
//#pragma omp parallel for num_threads(numofdomain)
	for (int outeriter = 0; outeriter < 100; ++outeriter) {
		vector<vector<double>> res(numofdomain);
		for (int domain = 0; domain < numofdomain; ++domain) {
			int elesize = d_elesize[domain];
			//��ȡ��ǰ����߽���ֵ

			//���ǵ�һ��߽�������װ��
			std::vector<std::vector<int>> locs(2, std::vector<int>(9 * (size_t)(elesize)));
			std::vector<double> vals(9 * (size_t)elesize);
			std::vector<double> F(d_dof[domain]);
			std::vector<double> F_bdr(d_dof[domain]);	//���򽻽紦��Ϊ��һ��߽������������õ����Ҳ�������
			int pos = 0;
			int linearelesize = 0, nonlinearsize = 0;

			for (int i_tri = 0; i_tri < elesize; ++i_tri) {
				int globaleleid = d_eleid[domain][i_tri];
				CTriElement triele = mp_triele[globaleleid];
				if (triele.material->getLinearFlag() == true)
					linearelesize++;
				else
					nonlinearsize++;
				for (int i = 0; i < 3; ++i) {
					int n1 = triele.n[i];	//ȫ��n1���
					int d_n1 = nparttable[domain][n1];	//�ֲ�δ������
					double F_bdr = 0;
					for (int j = 0; j < 3; ++j) {
						int n2 = triele.n[j];	//ȫ��n2���
						int d_n2 = nparttable[domain][n2];	//�ֲ�δ������
						//����װ��
						if (triele.material->getLinearFlag() == true) {
							if (crossbdrtable[domain][n1] != 1 && mp_node[n1].bdr != 1) {
								double mu = triele.material->getMu();
								double mut = mu * triele.xdot;
								double Se = triele.C[i][j] / mut;
								if (crossbdrtable[domain][n2] != 1 && mp_node[n2].bdr != 1) {
									locs[0][pos] = d_node_pos[domain][d_n1];	//ת�����ֲ���������
									locs[1][pos] = d_node_pos[domain][d_n2];
									vals[pos] = Se;
									++pos;
								}
								else {
									//�����һ��߽�����
									F_bdr += Se * mp_node[n2].At_old;
								}

							}
						}
					}
					if (crossbdrtable[domain][n1] != 1 && mp_node[n1].bdr != 1) {
						double Fe = triele.J * triele.area / 3;
						F[d_node_pos[domain][d_n1]] += Fe;
						F[d_node_pos[domain][d_n1]] -= F_bdr;
					}

				}
			}

			//���
			locs[0].resize(pos);
			locs[1].resize(pos);
			res[domain] = matsolver->solveMatrix(locs, vals, F, pos, d_dof[domain]);
		}

		//���������
		for (int domain = 0; domain < numofdomain; ++domain) {
			for (int n = 0; n < d_nodesize[domain]; ++n) {	//n�Ǿֲ�δ������
				int g_nodeid = d_nodeid[domain][n];	//ȫ�ֽڵ���
				if (crossbdrtable[domain][g_nodeid] != 1 && mp_node[g_nodeid].bdr != 1) {
					mp_node[g_nodeid].At = res[domain][d_node_pos[domain][n]];
				}
			}
		}

		//���´ų����
		updateB();

		//�ж����
		double error = 0, a = 0, b = 0;
		for (int i = 0; i < m_num_nodes; ++i) {
			a += (mp_node[i].At - mp_node[i].At_old) * (mp_node[i].At - mp_node[i].At_old);
			b += mp_node[i].At * mp_node[i].At;
		}
		error = sqrt(a) / sqrt(b);
		printf("Outer Iteration step: %d, error = %f\n", outeriter, error);
		if (error > maxerror) {
			for (int i = 0; i < m_num_nodes; ++i) {
				mp_node[i].At_old = mp_node[i].At;
			}
		}
		else {
			printf("Outer Iteration finish!\n");
		}
	}


}
