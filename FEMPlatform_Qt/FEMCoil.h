#pragma once

enum CoilType {
    CICULARCOIL,/** 圆形线圈 **/
    RECTANGULARCOIL,/** 方形线圈 **/
};

class FEMCoil
{
public:
    CoilType type{CoilType::CICULARCOIL};/** 0: CicularCoil;1:RectangularCoil **/
    int    xyz{0}; /* 0: parallel to x  1: parallel to y  2: parallel to z */
    int direction{1};/** 电流的流向：1或者-1 **/
    double r0{0};/** 内圆半径 **/
    double width{0};/** 线圈厚度 **/
    double height{0};/** 线圈高度 **/
    double Jor{0};
    double center[3]{0,0,0};/** 线圈的中心点坐标 **/
    double Nc{0};/** 线圈匝数 **/
    double Rc{0};/** 线圈电阻 **/
    double Vc{0};/** 线圈电压 **/
    double xWidth{0};/** 内部矩形x方向的宽度 **/
    double yWidth{0};/** 内部矩形y方向的宽度 **/
    double tau{0};/** 单位面积内的线圈匝数 **/
};