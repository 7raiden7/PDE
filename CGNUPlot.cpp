/*
 * CGNUPlot.cpp
 *
 *  Created on: 21 Dec 2014
 *      Author: raiden
 */

#include "CGNUPlot.h"

CGNUPlot::CGNUPlot()
{
	gnuplotpipe = popen("gnuplot -persistent", "w");
	if (!gnuplotpipe)
	{
		cerr << ("Gnuplot not found !");
	}
}

CGNUPlot::~CGNUPlot()
{
	fprintf(gnuplotpipe, "exit\n");
	pclose(gnuplotpipe);
}

void CGNUPlot::write_txt_2D(string filename, double* y_x, double* x, int M)
{
	string filename_x = filename + ".txt";
	ofstream myfile(filename_x.c_str());
	if (myfile.is_open())
	{
		for (int i = 0; i < M; ++i)
			myfile << x[i] << "\t" << y_x[i] << endl;

		myfile.close();
	}
	else
	{
		cout << "Unable to open file";
		throw;
	}
}

void CGNUPlot::write_txt_2D(string filename, double** z_y_x, double* x,
		double* y, int M, int N)
{
	string filename_x = filename + ".txt";
	ofstream myfile(filename_x.c_str());
	if (myfile.is_open())
	{
		for (int i = 0; i < M; ++i)
		{
			myfile << x[i] << ", ";
			for (int j = 0; j < N; ++j)
			{
				myfile << z_y_x[j][i] << ", ";
			}
			myfile << "\n";
		}

		myfile.close();
	}
	else
	{
		cout << "Unable to open file";
		throw;
	}
}

void CGNUPlot::write_txt_3D(const string &filename, double** z_x_y, double* x,
		double* y, int M, int N)
{
	string filename_x = filename + ".txt";
	ofstream myfile(filename_x.c_str());
	if (myfile.is_open())
	{
		for (int j = 0; j < N; ++j)
		{
			for (int i = 0; i < M; ++i)
			{
				myfile << x[i] << "\t" << y[j] << "\t" << z_x_y[j][i] << endl;
			}
			myfile << "\n";
		}

		myfile.close();
	}
	else
	{
		cout << "Unable to open file";
		throw;
	}
}

void CGNUPlot::write_txt_3D_animation(const string &filename, double** z_xy_t,
		double* x, double* y, int Mx, int My, int N)
{
	string filename_x = filename + ".txt";
	ofstream myfile(filename_x.c_str());

	int xy = 0;

	if (myfile.is_open())
	{
		for (int i = 0; i < Mx; ++i)
		{
			for (int j = 0; j < My; ++j)
			{
				xy = i + Mx * j;
				myfile << x[i] << ", " << y[j] << ", ";
				for (int t = 0; t < N; ++t)
					myfile << z_xy_t[t][xy] << ", ";
				myfile << endl;
			}
			myfile << endl;
		}

		myfile.close();
	}
	else
	{
		cout << "Unable to open file";
		throw;
	}
}

void CGNUPlot::delete_txt(string filename)
{
	string file_x = filename + ".txt";
	if (remove(file_x.c_str()) != 0)
	{
		cout << "Error deleting file" << endl;
		throw;
	}
}

void CGNUPlot::send_command(const string &command)
{
	fprintf(gnuplotpipe, "%s\n", command.c_str());
	fflush(gnuplotpipe);
}

void CGNUPlot::plot_2D_animation(double** z_y_x, double* x, double* y, int M,
		int N, string title, string xlabel, string ylabel, string legend,
		double delay)
{
	string aux_filename = "_aux";

	ostringstream d_str;
	d_str << delay;
	string delay_str = d_str.str();

	write_txt_2D(aux_filename, z_y_x, x, y, M, N);

	send_command("set title '" + title + "'");
	send_command("set xlabel '" + xlabel + "'");
	send_command("set style data linespoints");

	string command = "";
	for (int j = 0; j < N; ++j)
	{
		ostringstream strs;
		strs << y[j];
		string yj_str = strs.str();

		ostringstream str2s;
		str2s << j + 2;
		string y_str = str2s.str();

		send_command("set ylabel 'u(x, " + yj_str + ")'");
		string format_command = " notitle pt 7 ps .5";
		string plot_command = "plot '" + aux_filename + ".txt" + "' using 1:"
				+ y_str;

		send_command(plot_command + format_command);
		send_command("pause " + delay_str);
	}
}

// Plots multiple items on the same 2D plot
void CGNUPlot::plot_2D(double** z_y_x, double* x, double* y, int M, int N,
		string title, string xlabel, string ylabel, string legend)
{
	string aux_filename = "_aux";

	write_txt_2D(aux_filename, z_y_x, x, y, M, N);

	send_command("set title '" + title + "'");
	send_command("set xlabel '" + xlabel + "'");
	send_command("set ylabel '" + ylabel + "'");
	send_command("set style data linespoints");

	string command = "";
	for (int j = 0; j < N; ++j)
	{
		ostringstream strs;
		strs << y[j];
		string yj_str = strs.str();

		ostringstream str2s;
		str2s << j + 2;
		string y_str = str2s.str();

		//string format_command = " title '1:" + y_str + "' with line";
		string format_command = " notitle pt 7 ps .5";
		string aux = (j == 0 ? "plot" : "");
		string plot_command = aux + "'" + aux_filename + ".txt" + "' using 1:"
				+ y_str;

		command += plot_command + format_command + ", ";
	}

	send_command(command);
}

void CGNUPlot::plot_2D(double* y_x, double* x, int M, string title,
		string xlabel, string ylabel, string legend)
{
	string aux_filename = "_aux";

	write_txt_2D(aux_filename, y_x, x, M);

	send_command("set title '" + title + "'");
	send_command("set xlabel '" + xlabel + "'");
	send_command("set ylabel '" + ylabel + "'");
	send_command(
			"plot '" + aux_filename + ".txt" + "' title '" + legend
					+ "' with line lt 1 lc 7");
}

void CGNUPlot::plot_3D(double** z_y_x, double* x, double* y, int M, int N,
		string title, string xlabel, string ylabel, string zlabel,
		string legend)
{
	string aux_filename = "_aux";

	write_txt_3D(aux_filename, z_y_x, x, y, M, N);

	send_command("set title ' " + title + "'");
	send_command("set xlabel '" + xlabel + "'");
	send_command("set ylabel '" + ylabel + "'");
	send_command("set zlabel '" + zlabel + "'");
	send_command("set hidden3d");
	send_command("set ticslevel 0");
	send_command("set key box");
	send_command("set style data lines");
	send_command("set view 50, 45");
	send_command("set isosample 40");

	string plot_command = "splot '" + aux_filename + ".txt" + "' using 1:2:3"
			+ " title '" + legend;

	send_command(plot_command);
}

void CGNUPlot::plot_3D_animation(double** z_xy_t, double* x, double* y,
		double* t, int Mx, int My, int N, string title, string xlabel,
		string ylabel, string zlabel, string legend, double delay)
{
	string aux_filename = "_aux";
	ostringstream d_str;
	d_str << delay;
	string delay_str = d_str.str();

	send_command("set xlabel '" + xlabel + "'");
	send_command("set ylabel '" + ylabel + "'");
	send_command("set zlabel '" + zlabel + "'");
	send_command("set hidden3d");
	send_command("set view 50, 45");

	write_txt_3D_animation(aux_filename, z_xy_t, x, y, Mx, My, N);

	for (int n = 0; n < N; ++n)
	{
		ostringstream strs;
		strs << t[n];
		string t_str = strs.str();
		ostringstream strs2;
		strs2 << n + 3;
		string n_str = strs2.str();

		string plot_command = "splot '" + aux_filename + ".txt" + "' using 1:2:"
				+ n_str + " with lines" + " title 'u(x,y," + t_str + ")'";

		send_command(plot_command);
		send_command("pause " + delay_str);
	}
}

