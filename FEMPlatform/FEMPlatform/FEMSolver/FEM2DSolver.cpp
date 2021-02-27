#include "FEM2DSolver.h"

void FEM2DSolver::makeTrangle(int index)
{
	int k, m, n;
	double p0, p1, p2, q0, q1, q2, area;
	k = mp_triele[index].n[0];
	m = mp_triele[index].n[1];
	n = mp_triele[index].n[2];

	p0 = mp_node[m].y - mp_node[n].y;
	mp_triele[index].P[0] = p0;
	p1 = mp_node[n].y - mp_node[k].y;
	mp_triele[index].P[1] = p1;
	p2 = mp_node[k].y - mp_node[m].y;
	mp_triele[index].P[2] = p2;

	q0 = mp_node[n].x - mp_node[m].x;
	mp_triele[index].Q[0] = q0;
	q1 = mp_node[k].x - mp_node[n].x;
	mp_triele[index].Q[1] = q1;
	q2 = mp_node[m].x - mp_node[k].x;
	mp_triele[index].Q[2] = q2;

	area = 0.5 * std::abs(p1 * q2 - q1 * p2);
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

	//计算三角形重心半径
	if (flag == 2) {
		mp_triele[index].ydot = mp_triele[index].rc;
	}
	else {
		mp_triele[index].ydot = 1 / (mp_node[k].x + mp_node[m].x);
		mp_triele[index].ydot += 1 / (mp_node[k].x + mp_node[n].x);
		mp_triele[index].ydot += 1 / (mp_node[m].x + mp_node[n].x);
		mp_triele[index].ydot = 1.5 / mp_triele[index].ydot;
	}
}
