#include "FEM2DSolver.h"

void FEM2DSolver::solveMagneticForce()
{
	vector<bool> deformednode(m_num_nodes, false);
	vector<bool> movenode(m_num_nodes, false);
	vector<int> deformelement;
	vector<int> moveelement;

	//����ȫ����Ԫ����ǽڵ�
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		for (auto a : deformedlist) {
			if (mp_triele[i_tri].domain == a) {
				deformelement.push_back(i_tri);
				//�������ڵ���Ϊ�α�����ڵ�
				for (int i = 0; i < 3; ++i) {
					deformednode[mp_triele[i_tri].n[i]] = true;
				}
			}
		}
	}

	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		for (auto& a : movingmap) {
			if (mp_triele[i_tri].domain == a.first) {
				moveelement.push_back(i_tri);
				//�������ڵ���Ϊ�˶�����ڵ�
				for (int i = 0; i < 3; ++i) {
					movenode[mp_triele[i_tri].n[i]] = true;
				}
			}
		}
	}

	//����ȫ���α�����Ԫ�������Ԫ�ڽڵ�Ҳ�����˶���Ԫ��,
	//������Ӧ�ĵ�Ԫ������ͽڵ�����
	//ʦ�ִ������A��At�𣿣���
	for (auto i_tri : deformelement) {
		CTriElement triele = mp_triele[i_tri];
		double rc = triele.rc;
		double area = triele.area;
		double mu = triele.material->getMu(triele.B);
		double QA = triele.Q[0] * mp_node[triele.n[0]].At + triele.Q[1] * mp_node[triele.n[1]].At + triele.Q[2] * mp_node[triele.n[2]].At;
		double RA = triele.R[0] * mp_node[triele.n[0]].At + triele.R[1] * mp_node[triele.n[1]].At + triele.R[2] * mp_node[triele.n[2]].At;
		double tmp = PI / 4 / mu;	//ΪʲôҪ����4???
		double dvdB2 = triele.material->getdvdB2(triele.B);
		if (dvdB2 != 0) {
			cout << "i_tri: " << i_tri << endl;
		}
		double tmp1 = PI / 4 * (QA * QA + RA * RA) / area / rc * dvdB2 / 4;
		double dQdy[3][3];
		double dRdx[3][3];
		double ac = area * rc;
		double ac2 = ac * ac;
		double ac3 = ac2 * ac;

		dQdy[0][0] = 0; dQdy[0][1] = -1; dQdy[0][2] = 1;
		dQdy[1][0] = 1; dQdy[1][1] = 0; dQdy[1][2] = -1;
		dQdy[2][0] = -1; dQdy[2][1] = 1; dQdy[2][2] = 0;

		dRdx[0][0] = 0; dRdx[0][1] = 1; dRdx[0][2] = -1;
		dRdx[1][0] = -1; dRdx[1][1] = 0; dRdx[1][2] = 1;
		dRdx[2][0] = 1; dRdx[2][1] = -1; dRdx[2][2] = 0;

		double dSdx, dSdy;
		double dQAdy, dRAdx;
		double dBxdx, dBydx, dBxdy, dBydy;
		for (int i = 0; i < 3; ++i) {
			int n = mp_triele[i_tri].n[i];
			if (movenode[n] == true) {
				//x����
				dSdx = 0.5 * (dRdx[i][2] * triele.Q[1] - dRdx[i][1] * triele.Q[2]);
				dRAdx = (dRdx[i][0] * mp_node[triele.n[0]].At + dRdx[i][1] * mp_node[triele.n[1]].At + dRdx[i][2] * mp_node[triele.n[2]].At);
				mp_node[n].NodeForcex -= tmp * (2 * RA / ac * dRAdx);
				mp_node[n].NodeForcex -= tmp * ((QA * QA + RA * RA) * (-1 / ac2 * (area / 3 + rc * dSdx)));
				mp_node[n].NodeForcex -= tmp1 * (2 * RA / ac2 * dRAdx);
				mp_node[n].NodeForcex -= tmp1 * ((QA * QA + RA * RA) * (-2 / ac3 * (area / 3 + rc * dSdx)));
				//y����
				dSdy = 0.5 * (dQdy[i][1] * triele.R[2] - dQdy[i][2] * triele.R[1]);
				dQAdy = (dQdy[i][0] * mp_node[triele.n[0]].At + dQdy[i][1] * mp_node[triele.n[1]].At + dQdy[i][2] * mp_node[triele.n[2]].At);
				mp_node[n].NodeForcey -= tmp * ((QA * QA + RA * RA) * (-1) / ac2 * rc * dSdy);	//Ϊʲô���и��ţ�
				mp_node[n].NodeForcey -= tmp * (2 * QA / ac * dQAdy);
				mp_node[n].NodeForcey -= tmp1 * ((QA * QA + RA * RA) * (-2) / ac3 * rc * dSdy);
				mp_node[n].NodeForcey -= tmp1 * (2 * QA / ac2 * dQAdy);
			}
		}
	}

	////����ȫ���˶�����Ԫ�������Ԫ�ڽڵ�Ҳ�����α䵥Ԫ��,
	////������Ӧ�ĵ�Ԫ������ͽڵ�����
	//double dvdB2_sum = 0;
	//for (auto i_tri : moveelement) {
	//	CTriElement triele = mp_triele[i_tri];
	//	double rc = triele.rc;
	//	double area = triele.area;
	//	double mu = triele.material->getMu(triele.B);
	//	double QA = triele.Q[0] * mp_node[triele.n[0]].At + triele.Q[1] * mp_node[triele.n[1]].At + triele.Q[2] * mp_node[triele.n[2]].At;
	//	double RA = triele.R[0] * mp_node[triele.n[0]].At + triele.R[1] * mp_node[triele.n[1]].At + triele.R[2] * mp_node[triele.n[2]].At;
	//	double tmp = PI / 4 / mu;	//???
	//	double dvdB2 = triele.material->getdvdB2(triele.B);
	//	//cout << "i_tri: " << i_tri << ", dvdB2: " << dvdB2 << endl;
	//	dvdB2_sum += dvdB2;
	//	double tmp1 = PI / 4 * (QA * QA + RA * RA) / area / rc * dvdB2 / 4;
	//	double dQdy[3][3];
	//	double dRdx[3][3];
	//	double ac = area * rc;
	//	double ac2 = ac * ac;
	//	double ac3 = ac2 * ac;

	//	dQdy[0][0] = 0; dQdy[0][1] = -1; dQdy[0][2] = 1;
	//	dQdy[1][0] = 1; dQdy[1][1] = 0; dQdy[1][2] = -1;
	//	dQdy[2][0] = -1; dQdy[2][1] = 1; dQdy[2][2] = 0;

	//	dRdx[0][0] = 0; dRdx[0][1] = 1; dRdx[0][2] = -1;
	//	dRdx[1][0] = -1; dRdx[1][1] = 0; dRdx[1][2] = 1;
	//	dRdx[2][0] = 1; dRdx[2][1] = -1; dRdx[2][2] = 0;

	//	double dSdx, dSdy;
	//	double dQAdy, dRAdx;
	//	double dBxdx, dBydx, dBxdy, dBydy;
	//	for (int i = 0; i < 3; ++i) {
	//		int n = mp_triele[i_tri].n[i];
	//		if (deformednode[n] == true) {
	//			//x����
	//			dSdx = 0.5 * (dRdx[i][2] * triele.Q[1] - dRdx[i][1] * triele.Q[2]);
	//			dRAdx = (dRdx[i][0] * mp_node[triele.n[0]].At + dRdx[i][1] * mp_node[triele.n[1]].At + dRdx[i][2] * mp_node[triele.n[2]].At);
	//			mp_node[n].NodeForcex -= tmp * (2 * RA / ac * dRAdx);
	//			mp_node[n].NodeForcex -= tmp * ((QA * QA + RA * RA) * (-1 / ac2 * (area / 3 + rc * dSdx)));
	//			mp_node[n].NodeForcex -= tmp1 * (2 * RA / ac2 * dRAdx);
	//			mp_node[n].NodeForcex -= tmp1 * ((QA * QA + RA * RA) * (-2 / ac3 * (area / 3 + rc * dSdx)));
	//			//y����
	//			dSdy = 0.5 * (dQdy[i][1] * triele.R[2] - dQdy[i][2] * triele.R[1]);
	//			dQAdy = (dQdy[i][0] * mp_node[triele.n[0]].At + dQdy[i][1] * mp_node[triele.n[1]].At + dQdy[i][2] * mp_node[triele.n[2]].At);
	//			mp_node[n].NodeForcey -= tmp * ((QA * QA + RA * RA) * (-1) / ac2 * rc * dSdy);
	//			mp_node[n].NodeForcey -= tmp * (2 * QA / ac * dQAdy);
	//			mp_node[n].NodeForcey -= tmp1 * ((QA * QA + RA * RA) * (-2) / ac3 * rc * dSdy);
	//			mp_node[n].NodeForcey -= tmp1 * (2 * QA / ac2 * dQAdy);
	//		}
	//	}
	//}

	//�������������֮��
	Fx = 0, Fy = 0;
	for (int n = 0; n < m_num_nodes; ++n) {
		if (movenode[n] == true) {
			Fx += mp_node[n].NodeForcex;
			Fy += mp_node[n].NodeForcey;
		}
	}
	//cout << "dvdB2_sum: " << dvdB2_sum << endl;
	printf("Fx:%10.8e,Fy:%10.8e\n", Fx, Fy);
}

void FEM2DSolver::solveMagneticForce1()
{
	vector<bool> deformednode(m_num_nodes, false);
	vector<bool> movenode(m_num_nodes, false);
	vector<int> deformelement;
	vector<int> moveelement;

	//����ȫ����Ԫ����ǽڵ�
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		for (auto a : deformedlist) {
			if (mp_triele[i_tri].domain == a) {
				deformelement.push_back(i_tri);
				//�������ڵ���Ϊ�α�����ڵ�
				for (int i = 0; i < 3; ++i) {
					deformednode[mp_triele[i_tri].n[i]] = true;
				}
			}
		}
	}

	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		for (auto& a : movingmap) {
			if (mp_triele[i_tri].domain == a.first) {
				moveelement.push_back(i_tri);
				//�������ڵ���Ϊ�˶�����ڵ�
				for (int i = 0; i < 3; ++i) {
					movenode[mp_triele[i_tri].n[i]] = true;
				}
			}
		}
	}

	//����ȫ���α�����Ԫ�������Ԫ�ڽڵ�Ҳ�����˶���Ԫ��,
	//������Ӧ�ĵ�Ԫ������ͽڵ�����
	//ʦ�ִ������A��At�𣿣���
	for (auto i_tri : deformelement) {
		CTriElement triele = mp_triele[i_tri];
		double rc = triele.rc;
		double area = triele.area;
		double mu = triele.material->getMu(triele.B);
		double QAt = triele.Q[0] * mp_node[triele.n[0]].At + triele.Q[1] * mp_node[triele.n[1]].At + triele.Q[2] * mp_node[triele.n[2]].At;
		double RAt = triele.R[0] * mp_node[triele.n[0]].At + triele.R[1] * mp_node[triele.n[1]].At + triele.R[2] * mp_node[triele.n[2]].At;
		double tmp = PI / 4 / mu;	//ΪʲôҪ����4???
		double dvdB2 = triele.material->getdvdB2(triele.B);
		if (dvdB2 != 0) {
			cout << "i_tri: " << i_tri << endl;
		}
		double tmp1 = PI / 4 * (QAt * QAt + RAt * RAt) / area / rc * dvdB2 / 4;
		double dQdy[3][3];
		double dRdx[3][3];
		double ac = area * rc;
		double ac2 = ac * ac;
		double ac3 = ac2 * ac;

		dQdy[0][0] = 0; dQdy[0][1] = -1; dQdy[0][2] = 1;
		dQdy[1][0] = 1; dQdy[1][1] = 0; dQdy[1][2] = -1;
		dQdy[2][0] = -1; dQdy[2][1] = 1; dQdy[2][2] = 0;

		dRdx[0][0] = 0; dRdx[0][1] = 1; dRdx[0][2] = -1;
		dRdx[1][0] = -1; dRdx[1][1] = 0; dRdx[1][2] = 1;
		dRdx[2][0] = 1; dRdx[2][1] = -1; dRdx[2][2] = 0;

		double dSdx, dSdy;
		double dQAdy, dRAdx;
		double dBxdx, dBydx, dBxdy, dBydy;
		for (int i = 0; i < 3; ++i) {
			int n = mp_triele[i_tri].n[i];
			if (movenode[n] == true) {
				//x����
				dSdx = 0.5 * (dRdx[i][2] * triele.Q[1] - dRdx[i][1] * triele.Q[2]);
				dRAdx = (dRdx[i][0] * mp_node[triele.n[0]].At + dRdx[i][1] * mp_node[triele.n[1]].At + dRdx[i][2] * mp_node[triele.n[2]].At);
				mp_node[n].NodeForcex -= tmp * (2 * RAt / ac * dRAdx);
				mp_node[n].NodeForcex -= tmp * ((QAt * QAt + RAt * RAt) * (-1 / ac2 * (area / 3 + rc * dSdx)));
				mp_node[n].NodeForcex -= tmp1 * (2 * RAt / ac2 * dRAdx);
				mp_node[n].NodeForcex -= tmp1 * ((QAt * QAt + RAt * RAt) * (-2 / ac3 * (area / 3 + rc * dSdx)));
				//y����
				dSdy = 0.5 * (dQdy[i][1] * triele.R[2] - dQdy[i][2] * triele.R[1]);
				dQAdy = (dQdy[i][0] * mp_node[triele.n[0]].At + dQdy[i][1] * mp_node[triele.n[1]].At + dQdy[i][2] * mp_node[triele.n[2]].At);
				double dB2dy = -((RAt * RAt + QAt * QAt) * dSdy / 2 / area / area / area - QAt * dQAdy / 2 / area / area) / rc / rc;
				//mp_node[n].NodeForcey -= tmp * ((QAt * QAt + RAt * RAt) /** (-1)*/ / ac2 * rc * dSdy);
				mp_node[n].NodeForcey -= PI * rc * (QAt * QAt + RAt * RAt) / 4 / rc / rc / area / area * dSdy / mu;
				mp_node[n].NodeForcey -= PI * rc * area * dB2dy / mu;
				mp_node[n].NodeForcey -= PI * rc * area * (QAt * QAt + RAt * RAt) / 4 / rc / rc / area / area * dvdB2 * dB2dy;
				//mp_node[n].NodeForcey -= tmp * (2 * QAt / ac * dQAdy);
				//mp_node[n].NodeForcey -= tmp1 * ((QAt * QAt + RAt * RAt) * (-2) / ac3 * rc * dSdy);
				//mp_node[n].NodeForcey -= tmp1 * (2 * QAt / ac2 * dQAdy);
			}
		}
	}

	//����ȫ���˶�����Ԫ�������Ԫ�ڽڵ�Ҳ�����α䵥Ԫ��,
	//������Ӧ�ĵ�Ԫ������ͽڵ�����
	double dvdB2_sum = 0;
	for (auto i_tri : moveelement) {
		//CTriElement triele = mp_triele[i_tri];
		//double rc = triele.rc;
		//double area = triele.area;
		//double mu = triele.material->getMu(triele.B);
		//double QAt = triele.Q[0] * mp_node[triele.n[0]].At + triele.Q[1] * mp_node[triele.n[1]].At + triele.Q[2] * mp_node[triele.n[2]].At;
		//double RAt = triele.R[0] * mp_node[triele.n[0]].At + triele.R[1] * mp_node[triele.n[1]].At + triele.R[2] * mp_node[triele.n[2]].At;
		//double tmp = PI / 4 / mu;	//???
		//double dvdB2 = triele.material->getdvdB2(triele.B);
		////cout << "i_tri: " << i_tri << ", dvdB2: " << dvdB2 << endl;
		//dvdB2_sum += dvdB2;
		//double tmp1 = PI / 4 * (QAt * QAt + RAt * RAt) / area / rc * dvdB2 / 4;
		//double dQdy[3][3];
		//double dRdx[3][3];
		//double ac = area * rc;
		//double ac2 = ac * ac;
		//double ac3 = ac2 * ac;

		//dQdy[0][0] = 0; dQdy[0][1] = -1; dQdy[0][2] = 1;
		//dQdy[1][0] = 1; dQdy[1][1] = 0; dQdy[1][2] = -1;
		//dQdy[2][0] = -1; dQdy[2][1] = 1; dQdy[2][2] = 0;

		//dRdx[0][0] = 0; dRdx[0][1] = 1; dRdx[0][2] = -1;
		//dRdx[1][0] = -1; dRdx[1][1] = 0; dRdx[1][2] = 1;
		//dRdx[2][0] = 1; dRdx[2][1] = -1; dRdx[2][2] = 0;

		//double dSdx, dSdy;
		//double dQAdy, dRAdx;
		//double dBxdx, dBydx, dBxdy, dBydy;
		//for (int i = 0; i < 3; ++i) {
		//	int n = mp_triele[i_tri].n[i];
		//	if (deformednode[n] == true) {
		//		//x����
		//		dSdx = 0.5 * (dRdx[i][2] * triele.Q[1] - dRdx[i][1] * triele.Q[2]);
		//		dRAdx = (dRdx[i][0] * mp_node[triele.n[0]].At + dRdx[i][1] * mp_node[triele.n[1]].At + dRdx[i][2] * mp_node[triele.n[2]].At);
		//		mp_node[n].NodeForcex -= tmp * (2 * RAt / ac * dRAdx);
		//		mp_node[n].NodeForcex -= tmp * ((QAt * QAt + RAt * RAt) * (-1 / ac2 * (area / 3 + rc * dSdx)));
		//		mp_node[n].NodeForcex -= tmp1 * (2 * RAt / ac2 * dRAdx);
		//		mp_node[n].NodeForcex -= tmp1 * ((QAt * QAt + RAt * RAt) * (-2 / ac3 * (area / 3 + rc * dSdx)));
		//		//y����
		//		dSdy = 0.5 * (dQdy[i][1] * triele.R[2] - dQdy[i][2] * triele.R[1]);
		//		dQAdy = (dQdy[i][0] * mp_node[triele.n[0]].At + dQdy[i][1] * mp_node[triele.n[1]].At + dQdy[i][2] * mp_node[triele.n[2]].At);
		//		double dB2dy = -((RAt * RAt + QAt * QAt) * dSdy / 2 / area / area / area - QAt * dQAdy / 2 / area / area) / rc / rc;
		//		//mp_node[n].NodeForcey -= tmp * ((QAt * QAt + RAt * RAt) /** (-1) *// ac2 * rc * dSdy);
		//		mp_node[n].NodeForcey -= PI * rc * (QAt * QAt + RAt * RAt) / 4 / rc / rc / area / area * dSdy / mu;
		//		mp_node[n].NodeForcey -= PI * rc * area * dB2dy / mu;
		//		mp_node[n].NodeForcey -= PI * rc * area * (QAt * QAt + RAt * RAt) / 4 / rc / rc / area / area * dvdB2 * dB2dy;
		//		//mp_node[n].NodeForcey -= tmp * (2 * QAt / ac * dQAdy);
		//		//mp_node[n].NodeForcey -= tmp1 * ((QA * QA + RA * RA) * (-2) / ac3 * rc * dSdy);
		//		//mp_node[n].NodeForcey -= tmp1 * (2 * QA / ac2 * dQAdy);
		//	}
		//}
	}

	//�������������֮��
	Fx = 0, Fy = 0;
	for (int n = 0; n < m_num_nodes; ++n) {
		if (movenode[n] == true) {
			Fx += mp_node[n].NodeForcex;
			Fy += mp_node[n].NodeForcey;
		}
	}
	cout << "dvdB2_sum: " << dvdB2_sum << endl;
	printf("Fx:%10.8e,Fy:%10.8e\n", Fx, Fy);
}

void FEM2DSolver::makeTrangle()
{
	for (int index = 0; index < m_num_triele; ++index) {
		int k, m, n;
		double q0, q1, q2, r0, r1, r2, area;
		k = mp_triele[index].n[0];
		m = mp_triele[index].n[1];
		n = mp_triele[index].n[2];

		q0 = mp_node[m].y - mp_node[n].y;
		mp_triele[index].Q[0] = q0;
		q1 = mp_node[n].y - mp_node[k].y;
		mp_triele[index].Q[1] = q1;
		q2 = mp_node[k].y - mp_node[m].y;
		mp_triele[index].Q[2] = q2;

		r0 = mp_node[n].x - mp_node[m].x;
		mp_triele[index].R[0] = r0;
		r1 = mp_node[k].x - mp_node[n].x;
		mp_triele[index].R[1] = r1;
		r2 = mp_node[m].x - mp_node[k].x;
		mp_triele[index].R[2] = r2;

		area = 0.5 * std::abs(q1 * r2 - r1 * q2);
		mp_triele[index].area = area;

		mp_triele[index].rc = (mp_node[k].x +
			mp_node[m].x +
			mp_node[n].x) / 3;
		mp_triele[index].zc = (mp_node[k].y +
			mp_node[m].y +
			mp_node[n].y) / 3;

		int flag = 0;
		for (int f = 0; f < 3; f++) {
			if (mp_node[mp_triele[index].n[f]].x < 1e-7) {
				flag++;
			}
		}

		//�������������İ뾶
		if (flag == 2) {
			mp_triele[index].xdot = mp_triele[index].rc;
		}
		else {
			mp_triele[index].xdot = 1 / (mp_node[k].x + mp_node[m].x);
			mp_triele[index].xdot += 1 / (mp_node[k].x + mp_node[n].x);
			mp_triele[index].xdot += 1 / (mp_node[m].x + mp_node[n].x);
			mp_triele[index].xdot = 1.5 / mp_triele[index].xdot;
		}

		//����һ����������ԳƵ�Ԫϵ������
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				mp_triele[index].C[i][j] = ((mp_triele[index].R[i] * mp_triele[index].R[j] + mp_triele[index].Q[i] * mp_triele[index].Q[j])) / (4.0 * mp_triele[index].area);
			}
		}
		//printf("i_tri: %d, n0: %d, n1: %d, n2: %d, Y11: %f, Y12: %f, Y13: %f, Y22: %f, Y23: %f, Y33: %f\n", index, mp_triele[index].n[0], mp_triele[index].n[1], mp_triele[index].n[2], mp_triele[index].C[0][0], mp_triele[index].C[0][1], mp_triele[index].C[0][2], mp_triele[index].C[1][1], mp_triele[index].C[1][2], mp_triele[index].C[2][2]);
	}
}

void FEM2DSolver::processBoundaryCondition()
{
	//�����boundarymapֻ���˵�һ��߽�������ֻ�����ߵ�Ԫ��
	//1.��ȡ������ȫ����2.��ȡ���Ӧ�ĵ�Ԫ��ţ�3.��ȡ�ڵ�
	//4.ȥ�� 5.��node���ñ߽����� 6.���߽������ĩβ
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
	//���߽������ĩβ

	node_reorder.resize(m_num_nodes);
	node_pos.resize(m_num_nodes);
	int node_bdr = 0;
	for (int ibdr = 0; ibdr < m_num_nodes; ++ibdr) {
		if (mp_node[ibdr].bdr == 1) {
			++node_bdr;
			node_reorder[m_num_nodes - node_bdr] = ibdr;
			node_pos[ibdr] = m_num_nodes - node_bdr;
			mp_node[ibdr].A = 0;
		}
		else
		{
			node_reorder[ibdr - node_bdr] = ibdr;
			node_pos[ibdr] = ibdr - node_bdr;
		}
	}
}

void FEM2DSolver::processMaterial()
{
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		int domain = mp_triele[i_tri].domain;
		mp_triele[i_tri].material = materialmap[domain];
	}
}

void FEM2DSolver::processLoad()
{
	//���õ�Ԫ����
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		int domain = mp_triele[i_tri].domain;
		if (loadmap.find(domain) != loadmap.end()) {
			mp_triele[i_tri].J = loadmap[domain];
		}
	}
}

void FEM2DSolver::setCoilCurrent(double I)
{
}

void FEM2DSolver::updateB()
{
	for( int i_tri = 0; i_tri < m_num_triele; ++i_tri ) {
		double bx = 0, by = 0;
		for (int i = 0; i < 3; ++i) {
			int n = mp_triele[i_tri].n[i];
			bx += mp_triele[i_tri].R[i] * mp_node[n].At;
			by += mp_triele[i_tri].Q[i] * mp_node[n].At;
		}
		bx = bx / 2 / mp_triele[i_tri].area / mp_triele[i_tri].xdot;
		mp_triele[i_tri].Bx = bx;
		by = -by / 2 / mp_triele[i_tri].area / mp_triele[i_tri].xdot;
		mp_triele[i_tri].By = by;
		mp_triele[i_tri].B = sqrt(bx * bx + by * by);
	}
}

void FEM2DSolver::updateB(int i_tri)
{
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

