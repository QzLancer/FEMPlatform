#pragma once

enum CoilType {
    CICULARCOIL,/** Բ����Ȧ **/
    RECTANGULARCOIL,/** ������Ȧ **/
};

class FEMCoil
{
public:
    CoilType type{CoilType::CICULARCOIL};/** 0: CicularCoil;1:RectangularCoil **/
    int    xyz{0}; /* 0: parallel to x  1: parallel to y  2: parallel to z */
    int direction{1};/** ����������1����-1 **/
    double r0{0};/** ��Բ�뾶 **/
    double width{0};/** ��Ȧ��� **/
    double height{0};/** ��Ȧ�߶� **/
    double Jor{0};
    double center[3]{0,0,0};/** ��Ȧ�����ĵ����� **/
    double Nc{0};/** ��Ȧ���� **/
    double Rc{0};/** ��Ȧ���� **/
    double Vc{0};/** ��Ȧ��ѹ **/
    double xWidth{0};/** �ڲ�����x����Ŀ�� **/
    double yWidth{0};/** �ڲ�����y����Ŀ�� **/
    double tau{0};/** ��λ����ڵ���Ȧ���� **/
};