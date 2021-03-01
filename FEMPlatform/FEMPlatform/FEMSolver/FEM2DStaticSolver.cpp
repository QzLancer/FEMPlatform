#include "FEM2DStaticSolver.h"

void FEM2DStaticSolver::solve()
{
	//���������ε�Ԫ���β���
	for (int i = 0; i < m_num_triele; ++i) {
		makeTrangle(i);
	}

	//����߽����������Ϻ͸���
	processBoundaryCondition();
	processMaterial();
	processLoad();

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
			}
		}
	}
	A.resize(m_num_nodes);
	Bx.resize(m_num_triele);
	By.resize(m_num_triele);
	B.resize(m_num_triele);

	//������������
	if (nonlinearelesize == 0) {
		cout << "Linear Model" << endl;
		locs[0].resize(pos);
		locs[1].resize(pos);
		cout << "pos: " << pos << endl;
		double* res1 = matsolver->solveMatrix(locs, vals, F, pos, num_freenodes);
		for (int i = 0; i < num_freenodes; ++i) {
			int index = node_reorder[i];
			A[index] = res1[i] / mp_node[index].x;
			cout << A[index] << endl;
		}
		updateB();
		return;
	}

	//�����Բ��ֵ���
	int pos1 = pos;
	std::vector<double> F1 = F;
	cout << "NonLinear Model" << endl;
	for (int step = 0; step < maxitersteps; ++step) {
		pos = pos1;
		F = F1;
		for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
			if (mp_triele[i_tri].material->getLinearFlag() == false) {
				double mu, dvdb;
				mu = mp_triele[i_tri].material->getMu(B[i_tri]);
				dvdb = mp_triele[i_tri].material->getdvdB(B[i_tri]);

			}

		}

		//���´ų����
		updateB();
	}

}


