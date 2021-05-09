#include <QCoreApplication>
#include <iostream>
#include <vector>
#include "meshGFace.h"
#include "meshGRegionNetgen.h"
#include "meshGRegion.h"
#include "Generator.h"
#include "GModel.h"
#include "gmsh.h"
#include "MLine.h"
#include "Field.h"

using namespace std;
void remesh(int i);
void remesh(double dx, double dy);
void moveFace(GFace *f, double dx, double dy, double dz);
void updateField();
void deleteFaceMesh(GFace* f);

char fileName[] = "D:/femplatform/model/modelwithband";
int current_step = 1;
vector<double> dis = { 0, 0.0005, 0.0005, 0.0005, 0.0005 };
GModel* model;
int tag_remesh = 5;
int tag_xiantie = 1;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    int size = dis.size();

    //initialization
    int myargn = 4;
    char *myargv[] = { (char*)"gmsh",(char*)"-format",(char*)"msh2",(char*)"-v",(char*)"0" };
    gmsh::initialize(myargn, myargv);
    gmsh::option::setNumber("General.Terminal", 1);
    char geoName[256];
    sprintf(geoName, "%s.geo", fileName);
    gmsh::open(geoName);
    cout << "Opening model " << geoName << endl << endl;

    model = GModel::current();

    //第一次分网
    cout << "Generate first step mesh..." << endl;
    char nameMsh[256];
    GenerateMesh(model, 3);
    sprintf(nameMsh, "%s_%02d.msh", fileName, 0);
    gmsh::write(nameMsh);
    cout << "First step mesh finish." << endl << endl;

    for(int i = 0; i < size; ++i){
        remesh(i);
        current_step++;
    }
    gmsh::finalize();
        std::cout << "Hello World!\n";
    return a.exec();
}

void remesh(int i){
        printf("mesh step %d\n", i);
    remesh(0, dis[i]);
}

void remesh(double dx, double dy){
    char nameMsh[256];
    sprintf(nameMsh, "%s_%02d.msh", fileName, current_step);
    /** 未发生位移就不要分了 **/
    if (fabs(dx) < 1e-10 && fabs(dy) < 1e-10) {
        printf("Skipping remeshing.\n");
        gmsh::write(nameMsh);
        return;
    }

    /** 删掉空气的分网 **/
    GFace* f_xiantie = nullptr;
    GFace* f_air = nullptr;
    for (GModel::fiter it = model->firstFace(); it != model->lastFace(); ++it) {
        if ((*it)->tag() == tag_remesh) {
            f_air = (*it);
        }
        if ((*it)->tag() == tag_xiantie) {
            f_xiantie = (*it);
        }
    }

    /** 移动衔铁的分网 **/
    if (f_xiantie) {
        printf("moving armature region %d mesh dx=%lf,dy=%lf...\n", tag_xiantie, dx, dy);
        moveFace(f_xiantie, dx, dy, 0);
    }

    if (f_air) {
        /** 删除空气的分网 **/
        printf("deleting air surface %d\n", tag_remesh);
        deleteFaceMesh(f_air);
        /** 对空气进行重分网 **/
        printf("remesh air domain...\n");
        f_air->mesh(true);
    }
    gmsh::write(nameMsh);
    printf("Finish remesh 2d\n");
}

void moveFace(GFace *f, double dx, double dy, double dz){
    /** 修改内部分网节点坐标 **/
    for (std::size_t j = 0; j < f->getNumMeshVertices(); ++j) {
        f->getMeshVertex(j)->setXYZ(f->getMeshVertex(j)->x() + dx,
            f->getMeshVertex(j)->y() + dy,
            f->getMeshVertex(j)->z() + dz);
    }
    /** 修改边界棱上分网节点坐标 **/
    std::set<MVertex *, MVertexLessThanNum> all_vertices;
    for (auto e : f->edges()) {
        for (auto line : e->lines) {
            MVertex *v1 = line->getVertex(0);
            MVertex *v2 = line->getVertex(1);

            all_vertices.insert(v1);
            all_vertices.insert(v2);
        }
    }
    /** all_vertices比e->getMeshVertex多了顶点 **/
    for (std::set<MVertex *, MVertexLessThanNum>::iterator ite = all_vertices.begin(); ite != all_vertices.end(); ite++) {
        (*ite)->setXYZ((*ite)->x() + dx, (*ite)->y() + dy, (*ite)->z() + dz);
    }
    /** 修改几何顶点坐标 **/
    for (auto v : f->vertices()) {
        for (std::size_t j = 0; j < v->getNumMeshVertices(); ++j) {
            GPoint p(v->x() + dx, v->y() + dy, v->z() + dz);
            v->setPosition(p);
        }
    }
    /** 更新distance field **/
    updateField();
}

void updateField(){
    /** 更新distance field **/
    FieldManager* fM = GModel::current()->getFields();
    std::map<int, Field *>::iterator iter;
    for (iter = fM->begin(); iter != fM->end(); ++iter) {
        iter->second->update_needed = true;
        iter->second->update();
    }
}

void deleteFaceMesh(GFace* f){
    deMeshGFace dem;
    dem(f);
}
