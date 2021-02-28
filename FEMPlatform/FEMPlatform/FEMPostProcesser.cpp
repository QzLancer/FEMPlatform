#include "FEMPostProcesser.h"

void FEMPostProcesser::writeVtkFile(std::string _name, FEMMeshManager* _meshmanager, std::vector<double> _A)
{
	std::string name = std::string("../../result/") + _name + std::string(".vtk");
	FILE* fp = nullptr;
	int err;
	char ch[256];
	err = fopen_s(&fp, name.c_str(), "w");
	if (!fp) {
		std::cout << "Error: openning file!" << endl;
		exit(0);
	}
	/*
		 1: points
		 3: line
		 5: Triangular element
		 9: Quadrilateral element
		10: Tetrahedral element
		12: Hexahedral element
		13: Triangular prism element
		14: Pyramid element
	*/
	/** 数据版本声明 **/
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	/** 标题 **/
	fprintf(fp, "vtk title\n");
	/** 文件格式声明 **/
	fprintf(fp, "ASCII\n");
	/** 几何拓扑结构 **/
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	
	//节点
	int numofnodes = _meshmanager->getNumofNodes();
	CNode* node = _meshmanager->getNodes();
	fprintf(fp, "\nPOINTS %d float\n", numofnodes);
	for (int i = 0; i < numofnodes; ++i) {
		fprintf(fp, "%lf %lf %lf\n", node[i].x, node[i].y, node[i].z);
	}
	//一阶三角形单元
	int numoftriele = _meshmanager->getNumofTriEle();
	CTriElement* triele = _meshmanager->getTriElements();
	fprintf(fp, "\nCELLS %d %d\n", numoftriele, 4 * numoftriele);
	for (int i = 0; i < numoftriele; ++i) {
		fprintf(fp, "3 %d %d %d\n", triele[i].n[0], triele[i].n[1], triele[i].n[2]);
	}
	fprintf(fp, "\nCELL_TYPES %d\n", numoftriele);
	int type = 5;
	for (int i = 0; i < numoftriele; ++i) {
		fprintf(fp, "%d\n", type);
	}
	//节点磁矢位
	fprintf(fp, "\nPOINT_DATA %d\n", numofnodes);
	fprintf(fp, "SCALARS A double 1\n");
	fprintf(fp, "LOOKUP_TABLE %s\n", "Atable");
	for (int i = 0; i < numofnodes; ++i) {
		fprintf(fp, "%f\n", _A[i]);
	}

	fclose(fp);
}
