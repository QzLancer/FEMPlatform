#include "FEM2DStaticSolver.h"

void FEM2DStaticSolver::solve()
{
	////���������Ƿ���ȷ���䵽�������
	//cout << "Number of Nodes: " << m_num_nodes << endl;
	//cout << "Number of VtxElement: " << m_num_vtxele << endl;
	//cout << "Number of EdgElement: " << m_num_edgele << endl;
	//cout << "Number of TriElement: " << m_num_triele << endl;
	//for (int i = 0; i < m_num_nodes; ++i) {
	//	cout << "Nodes: " << mp_node[i].x << " , " << mp_node[i].y << endl;
	//}
	//for (int i = 0; i < m_num_vtxele; ++i) {
	//	cout << "VtxElement: " << mp_vtxele[i].n << ", Domain: " << mp_vtxele[i].domain << endl;
	//}
	//for (int i = 0; i < m_num_edgele; ++i) {
	//	cout << "EdgElement: " << mp_edgele[i].n[0] << ", " << mp_edgele[i].n[1] << ", Domain: " << mp_edgele[i].domain << endl;
	//}
	//for (int i = 0; i < m_num_triele; ++i) {
	//	cout << "TriElement: " << mp_triele[i].n[0] << ", " << mp_triele[i].n[1] << ", " << mp_triele[i].n[2] << ", Domain: " << mp_triele[i].domain << endl;
	//}
	//cout << "Number of Materialmap: " << materialmap.size() << endl;
	//cout << "Number of Loadmap: " << loadmap.size() << endl;
	//cout << "Number of Boundarymap : " << boundarymap.size() << endl;

	//���������ε�Ԫ���β���
	for (int i = 0; i < m_num_triele; ++i) {
		makeTrangle(i);
	}

	//���Ե�Ԫװ��


	//�����Ե���
	for (int iter = 0; iter < maxitersteps; ++iter) {
		//�����Ե�Ԫװ�����
		int pos = 0;
		

	}
}
