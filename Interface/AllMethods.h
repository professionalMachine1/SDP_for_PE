#pragma once

#ifndef ALLMETHODS_H
#define ALLMETHODS_H

#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include "TMatrix.h"

using namespace std;

class AllMethods
{
protected:
	TMatrix<double> V, R, F; //������� ������, ������� � �������� ������� � ������ ����� ��������������
	TVector<double> xBorder, yBorder; //������� �� ��� x � y ��������������
	TVector<double> resultTest, resultMain; //������ ������ � ������������ � ��������
	double h, k; //��� �� ��� x � y ��������������
	double hE, kE, A; //��������������� ������
	int TaskNumber, n, m; //����� ������, ����������� �����
	double pi;

public:
	AllMethods();
	AllMethods(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber); //����������� �������������

	double F_Function(double x, double y); //������� � ������ ����� (-f(x,y))
	void FunctionInicialisation();
	double ExactSolution(double x, double y);//������ �������
	double XInicialConditions(double x, int Num); //��������� ������� � ���� ������� ����������� ��� X
	double YInicialConditions(double y, int Num); //��������� ������� � ���� ������� ����������� ��� Y
	void Inicialisation(); //���������� ��������� ������� � ���� �������

	virtual void SetParametrs(); //�������� ������

	virtual void Runner(int index) { } //
	virtual double Runner(); //������������ �����

	TMatrix<double> GetV() { return V; }

	void GetRes(double **mas); //�������� ������� � ������������ ������

	virtual void MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name); //�������� �������, ��� �������� �����
	virtual void MethodError(double eps, int MaxIterations); //����������� �������, ��� �������� �����

	void VectorNevyazki();
	double NevyazkaInf(); //������� �� ����� �������������
	double NevyazkaEvkl(); //��������� ����� �������

	void ChangeGrid(int N, int M); //�������� �����
	void ChangeBorders(TVector<double> XBorder, TVector<double> YBorder); //�������� �������
	void ChangeTask(int TASKNumber); //�������� ������
	void Reset();//��������� �����������

	void XInterpolation();
	void YInterpolation();
	void Average();

	void SaveGrid(string s);
	void SaveData(string s);

	virtual void SetUserParametr(double Parametr) { }
	virtual double GetParametr() { return 0; }
	virtual void Optimal(bool wish) {}
	virtual TVector<double> GetTestResults() { return resultTest; }
	virtual TVector<double> GetMainResults() { return resultMain; }

	void TestPartOne(TMatrix<double>& U);
	void TestPartTwo(TMatrix<double>& U);

	void MainPartOne(AllMethods* Solution, char Name);
	void MainPartTwo(AllMethods* Solution, TMatrix<double>& V2);

	virtual double GetBet() { return 0; }
	
	double CalculateAccurancy(TMatrix<double>& matr);
};

#endif
