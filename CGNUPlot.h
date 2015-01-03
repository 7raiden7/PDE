#ifndef GNUPLOT_H_
#define GNUPLOT_H_

#include <string>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
class CGNUPlot
{
public:
	CGNUPlot();
	~CGNUPlot();
	void operator ()(const string & command);
	void send_command(const string & command);

	void plot_2D(double* y_x, double* x, int M, string title, string xlabel,
			string ylabel, string legend);

	void plot_2D(double** z_y_x, double* x, double* y, int M, int N,
			string title, string xlabel, string ylabel, string legend);

	void plot_2D_animation(double** z_y_x, double* x, double* y, int M, int N,
			string title, string xlabel, string ylabel, string legend, double delay);

	void plot_3D(double** z_y_x, double* x, double* y, int M, int N,
			string title, string xlabel, string ylabel, string zlabel,
			string legend);

	void plot_3D_animation(double** z_xy_t, double* x, double* y,
			double* t, int Mx, int My, int N, string title, string xlabel,
			string ylabel, string zlabel, string legend, double delay);
protected:
	FILE *gnuplotpipe;
private:
	void write_txt_2D(string file_name, double* y_x, double* x, int M);
	void write_txt_2D(string file_name, double** z_y_x, double* x, double* y, int M, int N);
	void write_txt_3D(const string &file_name, double** z_x_y, double* x,
			double* y, int M, int N);
	void write_txt_3D_animation(const string &file_name, double** z_xy_t, double* x,
			double* y, int Mx, int My, int N);
	void delete_txt(string file_name);
};

#endif
