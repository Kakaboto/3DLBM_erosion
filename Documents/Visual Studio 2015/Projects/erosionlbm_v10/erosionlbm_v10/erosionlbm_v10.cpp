// erosionlbm.cpp : Defines the entry point for the console application.
//
#pragma once

#include "stdafx.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "string.h"
#include "random"
#include "time.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <sstream>
//#include <Eigen\Dense>
//#include "mpi.h"







// DATA LAYOUT FOR f AND ftemp:
// f: 
// (i,j,k):		(0,0,0):(Nx-1,Ny-1,Nz-1) - physical points. (-1,-1,-1), (Nx, Ny, Nz) - Boundary points ("ghost points").

// solid_list:	-1 = fluid point. 0 = surface solid point. 1 = interior solid point.

//Frågor:	Varför kan jag inte printa ut värden från solid_list å f när dom är const?

// Antar att Nx = Ny = Nz (isotropt grid).

// rho = 1 i början, även i solid points. Detta löser sig med att macrovariables sätter rätt saker till 0 i första iterationen. rho = 1
// för att velocity BC behöver initial density. 

// varför är c = 1/sqrt(3)?

// antar alltid att objekt ligger i mitten av lådan (sfären går att välja plats för)

// erosion: varje solid point antas bestå av 100st korn. Each erosion step determines how many grains get eroded away (if  all 100 goes away then we say the point has eroded).




using namespace std;
typedef vector<double> dvec;
//typedef Eigen::Matrix<double, 3, 3> Rotmat;

double sq3inv = 0.577350269189626;
double sq2inv = 0.707106781186547;
double sq3 = 1.732050807568877;
double sq2 = 1.414213562373095;
double pi = 3.1415926535897932384626433;

int Nx = 80;
int Ny = 50;
int Nz = 50;
int edfforcedir = 0; //0 = no gravity. 1 = x. 2 = y. 3 = z.
					 // for the moment, only velocity BC in x direction.
string x0BC = "velocityBC";
string xNBC = "velocityBC";
string y0BC = "periodicBC";
string yNBC = "periodicBC";
string z0BC = "periodicBC";
string zNBC = "periodicBC";

double Bvel[3] = { 0.1, 0, 0 };

double latspace = 1.;

//double erosiontol = 4e-8;
double masspernode = 0.3; // Amount of mass per solid node. units = kg?
double kappa_er = 1.; // material property of solid. depends on toughness and density.
double VDW_0 = 7.e-6; //Wan-Der-Waals force
int updatefreq = 700;
double c = 1.;
double tau = 0.52;
double mu = c*c*(tau - 0.5);
double umax_theo = 0.1;
double Re = 1;
double Lz = (Nz - 1)*latspace;
//double gg = mu*Re / (tau*L);
double gg = 0.001;
double F = mu * 10 / Lz;
int cellist[27] = { 3, 2, 3, 2, 1, 2, 3, 2, 3  ,  2, 1, 2, 1, 0, 1, 2, 1, 2  ,  3, 2, 3, 2, 1, 2, 3, 2, 3 };
double Delta_T = 10.*latspace/Bvel[0]; //Need Delta_T to be bigger than characteristic time. Hence, the factor 10 in front of L/U.
double dt = 1.;



class momentum_direction {
public:
	int **  e;
	momentum_direction() {
		e = new int*[27];
		//		for (auto &it : e) {
		//			it.resize(3);
		//		}
		for (int i = 0; i < 27; i++)
			e[i] = new int[3];
		int ix = 0;
		int iy = 0;
		int iz = 0;
		int a = 0;
		for (iz = -1; iz < 2; iz++) {
			for (iy = -1; iy < 2; iy++) {
				for (ix = -1; ix < 2; ix++) {
					e[a][0] = ix;
					e[a][1] = iy;
					e[a][2] = iz;
					//					cout << " " << e[a][0] << " " << e[a][1] << " " << e[a][2] << "\n";
					a++;
				}
			}
		}
	}

	int& operator()(int a, int i) {
		return e[a][i];
	}

private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;

};


class Grid {

public:
	dvec x;			// 0:Nx-1
	dvec y;			// 0:Nx-1
	dvec z;			// 0:Nx-1
	Grid() {
		x.resize(Nx);
		y.resize(Nx);
		z.resize(Nz);
		int ix = 0;
		int iy = 0;
		int iz = 0;
		//		latspace_x = 1. / (Nx - 1.);
		//		latspace_y = 1. / (Ny - 1.);
		//		latspace_z = 1. / (Nz - 1.);

		for (ix = 0; ix < Nx; ix++)
			x[ix] = xmin + latspace*ix;
		for (iy = 0; iy < Ny; iy++)
			y[iy] = ymin + latspace*iy;
		for (iz = 0; iz < Nz; iz++)
			z[iz] = zmin + latspace*iz;
		cout << "\n Grid created. \n";

	}
	void printgrid(FILE * gridfile) {
		int ix = 0;
		int iy = 0;
		int iz = 0;
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					fprintf(gridfile, "%e %e %e\n", x[ix], y[iy], z[iz]);
				}
			}
		}
		cout << "\n Printed grid to file for figures. \n";
	}
	~Grid();


private:
	double xmin = 0; double xmax = 1;
	double ymin = 0; double ymax = 1;
	double zmin = 0; double zmax = 1;
	double latspace_x = 0;
	double latspace_y = 0;
	double latspace_z = 0;
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;

};

Grid::~Grid() {
	cout << "\n Grid deleted. \n";
}



class Solid_list {

public:
	vector<int> element;		// solid_list:	-1 = fluid point. 0 = surface solid point. 1 = interior solid point.
	Solid_list(int choice, Grid& grid, double rotation, momentum_direction& e, FILE * parfile);
	~Solid_list();
	int& operator()(int ix, int iy, int iz) {
		return element[(ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot];
	}

	double center[3];
	//double rotation = pi*0.5;

	void updatesurface(int ix, int iy, int iz, momentum_direction& e) {
		int a = 0;
		int n = 0;
		element[(ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot] = -1;
		for (a = 0; a < 27; a++) {
			n = (ix + 1 + e(a, 0)) + (iy + 1 + e(a, 1))*Nxtot + (iz + 1 + e(a, 2))*Nxtot*Nytot;
			if (element[n] == 1)
				element[n] = 0;
		}
	}

	void printsolid_list(FILE * solfile) {
		int ix = 0;
		int iy = 0;
		int iz = 0;
		int n = 0;
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
					fprintf(solfile, "%i ", element[n]);
				}
			}
		}
		fprintf(solfile, "\n");
		cout << "\n Printed grid to file. \n";
	}

	void clear_list() {
		int n = 0;
		for (int iz = -1; iz < Nz + 1; iz++) {
			for (int iy = -1; iy < Ny + 1; iy++) {
				for (int ix = -1; ix < Nx + 1; ix++) {
					n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
					element[n] = 0.;
				}
			}
		}
	}

private:
	//	double latspace = 1. / (Nx - 1); //assuming Nx = Ny = Nz
	double xmin = 0; double xmax = 1;
	double ymin = 0; double ymax = 1;
	double zmin = 0; double zmax = 1;
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;

};
// Unnessesary to fill list with 0:s at the start.
Solid_list::Solid_list(int choice, Grid& grid, double rotation_x, momentum_direction& e, FILE * parfile) {
	element.resize(Ncubedtot);
	fill(element.begin(), element.end(), 0);
	double xpos = 0.3;
	double ypos = 0.5;
	fprintf(parfile, "%e %e ", xpos, ypos);
	//for case 7-----
	double funcvalue = 0;
	int i = 0;
	double triangle_like = 16;
	double y_stretch = 0.6;
	double length = 20;
	double sizelimit = 10000;
	//--------------
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;
	int n = 0;
	int nnext;
	double xmid = 0;
	double ymid = 0;
	int surfacecheck = 0;
	string tempstr;
	double radius1;
	double radius2;
	double radiussq;
	double radiussqh;
	double distsq;
	//	dvec center;
	//	center.resize(3);
	//	Rotmat Rx;
	//	Rx(0, 0) = 1; Rx(0, 1) = 0; Rx(0, 2) = 0;
	//	Rx(1, 0) = 0; Rx(1, 1) = cos(rotation_x); Rx(2, 1) = -sin(rotation_x);
	//	Rx(2, 0) = 0; Rx(2, 1) = sin(rotation_x); Rx(2, 2) = sin(rotation_x);

	switch (choice) {
	case 1:	//sphere
		cout << "\n Sphere chosen! \n";
		cout << "\n Input radius: ";
		//cin >> radius1;
		radius1 = 13.;
		//		stringstream(tempstr) >> radius;
		//		if (radius1 * 2 > 1) {
		//			cout << "Radius too large! must be between 0 and 1";
		//			break;
		//		}
		cout << "\n Input center (x,y,z goes from (0,0,0) to (1,1,1)): ";
		//cin >> center[0];
		//cin >> center[1];
		//cin >> center[2];
		center[0] = (Nx - 1)*xpos;
		center[1] = (Ny - 1)*ypos;
		center[2] = (Nz - 1)*0.5;
		/*		getline(cin, tempstr);
		stringstream(tempstr) >> center[0];
		getline(cin, tempstr);
		stringstream(tempstr) >> center[1];
		getline(cin, tempstr);
		stringstream(tempstr) >> center[2];*/
		//		if ((center[0] > 1) || (center[1] > 1) || (center[2] > 1)) {
		//			cout << "Center outside of grid! must be between 0 and 1";
		//			break;
		//		}


		radiussq = radius1*radius1;
		radiussqh = (radius1 - sq3*latspace)*(radius1 - sq3*latspace);
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					distsq = (grid.x[ix] - center[0])*(grid.x[ix] - center[0]) + (grid.y[iy] - center[1])*(grid.y[iy] - center[1]) + (grid.z[iz] - center[2])*(grid.z[iz] - center[2]);
					n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
					if (distsq > radiussq)
						element[n] = -1;
					else
						element[n] = 1;
				}
			}
		}
		cout << "\n Created obstacle: Sphere! \n";
		n = 0;
		break;

	case 2:
		cout << "Cylinder chosen! \n Input radius: ";
		getline(cin, tempstr);
		stringstream(tempstr) >> radius1;
		if (radius1 * 2 > 1) {
			cout << "Radius too large! must be between 0 and 1";
			break;
		}
		radiussq = radius1*radius1;
		radiussqh = (radius1 + sq3*latspace)*(radius1 + sq3*latspace);
		xmid = floor(((double)Nx - 1.)*0.5);
		ymid = floor(((double)Ny - 1.)*0.5);
		n = 0;
		for (iz = 0; iz < Nz; iz++) {	// can optimize by doing this for 1 iz, then just copiyng it to all other places in the vector.
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					distsq = (grid.x[ix] - grid.x[xmid])*(grid.x[ix] - grid.x[xmid]) + (grid.y[iy] - grid.y[ymid])*(grid.y[iy] - grid.y[ymid]);
					n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
					if (distsq < radiussq)
						element[n] = -1;
					else
						element[n] = 1;
				}
			}
		}
		n = 0;
		break;
	case 3:
		cout << "\n Square pipe chosen! \n";
		cout << "flow in z direction \n";
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
					if ((ix == 0) || (ix == Nx - 1) || (iy == 0) || (iy == Ny - 1))
						element[n] = 1;
					else
						element[n] = -1;
				}
			}
		}
		n = 0;
		break;

	case 4:
		cout << "Cylinder chosen! \n Input radius: ";
		getline(cin, tempstr);
		stringstream(tempstr) >> radius1;
		if (radius1 * 2 > 1) {
			cout << "Radius too large! must be between 0 and 1";
			break;
		}
		radiussq = radius1*radius1;
		radiussqh = (radius1 + sq3*latspace)*(radius1 + sq3*latspace);
		xmid = floor(((double)Nx - 1.)*0.5);
		ymid = floor(((double)Ny - 1.)*0.5);
		n = 0;
		for (iz = 0; iz < Nz; iz++) {	// can optimize by doing this for 1 iz, then just copiyng it to all other places in the vector.
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					distsq = (grid.x[ix] - grid.x[xmid])*(grid.x[ix] - grid.x[xmid]) + (grid.y[iy] - grid.y[ymid])*(grid.y[iy] - grid.y[ymid]);
					n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
					if (distsq < radiussq)
						element[n] = -1;
					else
						element[n] = 1;
				}
			}
		}

		cout << "\n Sphere chosen! \n";
		cout << "\n Input radius: ";
		cin >> radius2;
		//		stringstream(tempstr) >> radius;
		if (radius2 * 2 > radius1 * 2) {
			cout << "Radius too large! must be between 0 and " << radius1 * 2 << "\n";
			break;
		}
		cout << "\n Input center (x,y,z goes from (0,0,0) to (1,1,1)): ";
		cin >> center[0];
		cin >> center[1];
		cin >> center[2];
		if ((center[0] > 1) || (center[1] > 1) || (center[2] > 1)) {
			cout << "Center outside of grid! must be between 0 and 1";
			break;
		}

		n = 0;
		radiussq = radius2*radius2;
		radiussqh = (radius2 - sq3*latspace)*(radius2 - sq3*latspace);
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					distsq = (grid.x[ix] - center[0])*(grid.x[ix] - center[0]) + (grid.y[iy] - center[1])*(grid.y[iy] - center[1]) + (grid.z[iz] - center[2])*(grid.z[iz] - center[2]);
					n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
					if (distsq > radiussq)
						if (element[n] != 1)
							element[n] = -1;
						else
							element[n] = 1;
				}
			}
		}
		cout << "\n Created obstacle: Sphere! \n";
		break;

	case 5:
		cout << "Cylinder object chosen! \n Input radius: ";
		getline(cin, tempstr);
		stringstream(tempstr) >> radius1;
		//		if (radius1 * 2 > 1) {
		//			cout << "Radius too large! must be between 0 and 1";
		//			break;
		//		}
		radiussq = radius1*radius1;
		radiussqh = (radius1 + sq3*latspace)*(radius1 + sq3*latspace);
		xmid = ((double)Nx - 1.)*xpos;
		ymid = ((double)Ny - 1.)*ypos;
		n = 0;
		for (iz = 0; iz < Nz; iz++) {	// can optimize by doing this for 1 iz, then just copiyng it to all other places in the vector.
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					distsq = ((double)grid.x[ix] - xmid)*((double)grid.x[ix] - xmid) + ((double)grid.y[iy] - ymid)*((double)grid.y[iy] - ymid);
					n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
					if (distsq > radiussq)
						element[n] = -1;
					else
						element[n] = 1;
				}
			}
		}
		break;
	case 6:
		n = 0;
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
					element[n] = -1;
				}
			}
		}
		n = 0;
		break;
	case 7: //dolphin case
		n = 0;
		double iytemp;
		double iztemp;
		center[0] = (Nx - 1)*0.5;
		center[1] = (Ny - 1)*0.5;
		center[2] = (Nz - 1)*0.5;
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
					if (ix < (center[0] + length*0.5) && ix >(center[0] - length*0.5)) {
						funcvalue = 0;
						iytemp = (double(iy) - center[1])*cos(rotation_x) - (double(iz) - center[2])*sin(rotation_x);
						iztemp = (double(iy) - center[1])*sin(rotation_x) + (double(iz) - center[2])*cos(rotation_x);
						for (i = 1; i < 4; i++) { // changed from i = 0; i < 3 to i = 1; i < 4.
												  // deforms a sphere continuously to a triangle. higher m = more triangle:y.
							funcvalue += 1e-7*pow(abs(y_stretch*(iytemp)*cos(2. * pi*i / 3. + pi / 2) + (iztemp)*sin(2. * pi*i / 3. + pi / 2) - 1. / 3.), triangle_like);
						}
						if (funcvalue > sizelimit)
							element[n] = -1;
						else
							element[n] = 1;
					}
					else
						element[n] = -1;
				}
			}
		}


	default:
		cout << "\n Error. incorrect input for obstacle type. \n";
	}



	// Placing out surface points.
	for (iz = 0; iz < Nz; iz++) {	// can optimize by doing this for 1 iz, then just copiyng it to all other places in the vector.
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				n = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot;
				if (element[n] == 1) {
					for (a = 0; a < 27; a++) {
						if (((ix + e(a, 0)) >= 0) && ((ix + e(a, 0)) < Nx) && ((iy + e(a, 1)) >= 0) && ((iy + e(a, 1)) < Ny) && ((iz + e(a, 2)) >= 0) && ((iz + e(a, 2)) < Nz)) {
							nnext = (ix + 1 + e(a, 0)) + (iy + 1 + e(a, 1))*Nxtot + (iz + 1 + e(a, 2))*Nxtot*Nytot;
							if (element[nnext] == -1) {
								element[n] = 0;
								break;
							}
						}
					}
				}
			}
		}
	}
}

Solid_list::~Solid_list() {
	//	element.clear;
	cout << "Object deleted. \n";
}



class Surface_fluid_list {
public:
	vector<int> element;
	Surface_fluid_list(Solid_list& solid_list, momentum_direction& e);
	~Surface_fluid_list();

	int& operator()(int ix, int iy, int iz) {
		return element[(ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot];
	}

private:
	//	double latspace = 1. / (Nx - 1); //assuming Nx = Ny = Nz
	double xmin = 0; double xmax = 1;
	double ymin = 0; double ymax = 1;
	double zmin = 0; double zmax = 1;
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};

Surface_fluid_list::Surface_fluid_list(Solid_list& solid_list, momentum_direction& e) {
	element.resize(Ncubedtot);
	fill(element.begin(), element.end(), -1);

	int ix = 0;
	int iy = 0;
	int iz = 0;
	int ixshift = 0;
	int iyshift = 0;
	int izshift = 0;
	int a = 0;
	int nn = 0;

	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {

				if (solid_list(ix, iy, iz) == 0) { //if we're at a surface node
					for (a = 0; a < 27; a++) {	   //check all nearest neighbours.
						ixshift = ix + e(a, 0);
						iyshift = iy + e(a, 1);
						izshift = iz + e(a, 2);
						if (solid_list(ixshift, iyshift, izshift) == -1) { //if nn is a fluid node, it surface fluid node.
							nn = (ixshift + 1) + (iyshift + 1)*Nxtot + (izshift + 1)*Nxtot*Nytot;
							element[nn] = 1;
						}

					}
				}
			}
		}
	}
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				nn = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nztot;
				if (element[nn] == 1) { // if we're at a surface fluid node
					for (a = 0; a < 27; a++) {	   //check all nearest neighbours.
						ixshift = ix + e(a, 0);
						iyshift = iy + e(a, 1);
						izshift = iz + e(a, 2);
						if (solid_list(ixshift, iyshift, izshift) == -1 && element[nn] != 1) { //if nn is a fluid and not already a surface fluid point.
							nn = (ixshift + 1) + (iyshift + 1)*Nxtot + (izshift + 1)*Nxtot*Nytot;
							element[nn] = 2;
						}

					}
				}
			}
		}
	}


}

Surface_fluid_list::~Surface_fluid_list() {
	//	element.clear;
	cout << "\n Surface fluid List deleted \n";
}



class direction_density {
public:
	dvec element;		//size: 27*Ncubedtot
	direction_density() {
		element.resize(27 * Ncubedtot);
	}

	int index(int ix, int iy, int iz, int a) {
		int ind = (ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot + a*Ncubedtot;
		return ind;
	}

	double& operator()(int ix, int iy, int iz, int a) {
		return element[(ix + 1) + (iy + 1)*Nxtot + (iz + 1)*Nxtot*Nytot + a*Ncubedtot];
	}

	void clearvector() {
		// Fortsätt här!!
	}
private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;

};

class EDF {
public:
	dvec element;		//size: 27*Ncubed

	EDF() {
		element.resize(27 * Ncubed);
		fill(element.begin(), element.end(), 0);
	};

	int index(int ix, int iy, int iz, int a) {
		int ind = ix + iy*Nx + iz*Nx*Ny + a*Ncubed;
		return ind;
	}

	double& operator()(int ix, int iy, int iz, int a) {
		return element[ix + iy*Nx + iz*Nx*Ny + a*Ncubed];
	}
private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};

class velocity {
public:
	dvec element;		//size: 3*Ncubed


	velocity() {
		//element.resize(3 * Ncubed);
		element.resize(3 * Ncubed);
	};

	int index(int ix, int iy, int iz, int i) {
		int ind = ix + iy*Nx + iz*Nx*Ny + i;
		return ind;
	}

	double& operator()(int ix, int iy, int iz, int i) {
		return element[3 * (ix + iy*Nx + iz*Nx*Ny) + i];
	}

private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};

class density {
public:
	dvec element;		//size: Ncubed

	density() {
		element.resize(Ncubed);
		fill(element.begin(), element.end(), 1);
	};

	void clear() {
		fill(element.begin(), element.end(), 0.);
	}

	int index(int ix, int iy, int iz) {
		int ind = ix + iy*Nx + iz*Nx*Ny;
		return ind;
	}

	double& operator()(int ix, int iy, int iz) {
		return element[ix + iy*Nx + iz*Nx*Ny];
	}
	void operator()(int ix, int iy, int iz, double num) {
		element[ix + iy*Nx + iz*Nx*Ny] = num;
	}
private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};

class Stresstensor {
public:
	dvec element;
	Stresstensor() {
		element.resize(9 * Ncubed);
		fill(element.begin(), element.end(), 0);
	}
	~Stresstensor() {
		//		element.clear;
	}
	void clear() {
		fill(element.begin(), element.end(), 0);
	}

	double& operator()(int ix, int iy, int iz, int i, int j) {
		return element[(ix + iy*Nx + iz*Nx*Ny) * 9 + i * 3 + j];
	}
private:
	int Ncubed = Nx*Ny*Nz;

};

// same structure as velocity class
class Normalvector {
public:
	dvec element;
	Normalvector() {
		element.resize(3 * Ncubed);
		fill(element.begin(), element.end(), 0);
	}
	~Normalvector() {
		//		element.clear;
	}

	double& operator()(int ix, int iy, int iz, int i) {
		return element[(ix + iy*Nx + iz*Nx*Ny) * 3 + i];
	}
private:
	int Ncubed = Nx*Ny*Nz;
};
// same structure as velocity class
class Wall_force {
public:
	dvec element;		//size: Ncubed

	Wall_force() {
		element.resize(3 * Ncubed);
	};

	void clear() {
		fill(element.begin(), element.end(), 0.);
	}

	double& operator()(int ix, int iy, int iz, int i) {
		return element[(ix + iy*Nx + iz*Nx*Ny) * 3 + i];
	}

private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};

class vector3Ncubed {
public:
	dvec element;		//size: 3*Ncubed


	vector3Ncubed() {
		//element.resize(3 * Ncubed);
		element.resize(3 * Ncubed);
		fill(element.begin(), element.end(), 0.);
	};

	void clear() {
		fill(element.begin(), element.end(), 0.);
	}

	int index(int ix, int iy, int iz, int i) {
		int ind = 3 * (ix + iy*Nx + iz*Nx*Ny) + i;
		return ind;
	}

	double& operator()(int ix, int iy, int iz, int i) {
		return element[3 * (ix + iy*Nx + iz*Nx*Ny) + i];
	}

	double maxsq() {
		double temp[3] = { 0 };
		double max[3] = { 0 };
		int n = 0;
		for (int iz = 0; iz < Nz; iz++) {
			for (int iy = 0; iy < Ny; iy++) {
				for (int ix = 0; ix < Nx; ix++) {
					n = 3 * (ix + iy*Nx + iz*Nx*Ny);
					temp[0] = element[n];
					temp[1] = element[n + 1];
					temp[2] = element[n + 2];
					if (abs(temp[0]) > max[0])
						max[0] = temp[0];
					if (abs(temp[1]) > max[1])
						max[1] = temp[1];
					if (abs(temp[2]) > max[2])
						max[2] = temp[2];
				}
			}
		}
		return max[0] * max[0] + max[1] * max[1] + max[2] * max[2];
	}
private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};

class vectorNcubed {
public:
	dvec element;		//size: Ncubed

	vectorNcubed() {
		element.resize(Ncubed);
		fill(element.begin(), element.end(), 0.);
	};

	void clear() {
		fill(element.begin(), element.end(), 0.);
	}

	int index(int ix, int iy, int iz) {
		int ind = ix + iy*Nx + iz*Nx*Ny;
		return ind;
	}

	double& operator()(int ix, int iy, int iz) {
		return element[ix + iy*Nx + iz*Nx*Ny];
	}

private:
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;
};


//... should prob försöka lägga alla vektorer i en å samma klass med massa olika () operatorer bara.


void IC(const momentum_direction& e, direction_density& f, direction_density& ftemp, Solid_list& solid_list);

void macrovariables(velocity& u, density& rho, Solid_list& solid_list, direction_density& f, momentum_direction& e);
double find_umax(velocity& u);
double find_uav(velocity& u);

void edf(Solid_list& solid_list, velocity& u, density& rho, EDF& feq, momentum_direction& e, int forcedirection);
double calcfeq(density& rho, momentum_direction& e, double ueq[3], int ix, int iy, int iz, int a, int forcedirection);
double * calcforce(double ueq[3], density& rho, int ix, int iy, int iz, int forcedirection);
double edfvecdot(momentum_direction& e, double ueq[3], int a);

void stream(Solid_list& solid_list, direction_density& f, direction_density& ftemp, momentum_direction& e);
void bouncebackBC(direction_density& ftemp, direction_density& f, momentum_direction& e, int ix, int iy, int iz, int a, int ashift);

void collision(Solid_list& solid_list, direction_density& f, direction_density& ftemp, EDF& feq, momentum_direction& e);

void computestress(momentum_direction& e, direction_density& ftemp, direction_density& f, EDF& feq, Solid_list& solid_list, Stresstensor& stresstensor, Normalvector& nhat, Wall_force& tau_stress, vector3Ncubed& F_sum, density& masschange, vectorNcubed& F_vdw, int i_er, int i_Fvdw, density& rho);
double VDWforce(int x, int y, int z);

void updatePBC(direction_density& f);
void updateBC(direction_density& f, int t, double Bvel[3], density& rho, momentum_direction& e, velocity& u);
void updatePBC_solid(Solid_list& solid_list);
void printstuff(FILE * velfile, FILE * densfile, FILE * parfile, FILE * reyfile, FILE * stressfile, FILE * forcefile, FILE * nhatfile, FILE * sttensfile, FILE * torfile, FILE * erodefile, int Renum, velocity& u, density& rho, Wall_force& tau_stress, vector3Ncubed& F_D, Normalvector& nhat, Stresstensor& stresstensor, vector3Ncubed& torque, density& erodelist);
vector<int> readBC(string x0BC, string xNBC, string y0BC, string yNBC, string z0BC, string zNBC);

void computetorque(Solid_list& solid_list, Wall_force& tau_stress, vector3Ncubed& torque);
dvec crossproduct(Wall_force& tau_stress, double r[3], int ix, int iy, int iz);

void erosion(Solid_list& solid_list, momentum_direction& e, vector3Ncubed& F_sum, density& rho, direction_density& f, int forcedirection, FILE * solfile, density& erodelist, Normalvector& nhat, vectorNcubed& ero_reso_check, vectorNcubed& F_vdw, FILE * errorfile);
void erodepoint(int ix, int iy, int iz, density& masschange, Solid_list& solid_list, density& rho, Normalvector& nhat, direction_density& f, momentum_direction& e);
vector<int> BCtype = readBC(x0BC, xNBC, y0BC, yNBC, z0BC, zNBC);
int main()
{
	//===========================================================================================================
	// Initializing variables
	FILE * velfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\stresserror\\velocity_v9.txt", "w");
	FILE * densfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\stresserror\\density_v9.txt", "w");
	FILE * gridfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\stresserror\\grid_v9.txt", "w");
	FILE * parfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\stresserror\\parameters_v9.txt", "w");
	FILE * solfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\stresserror\\solid_list_v9.txt", "w");
	FILE * reyfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\stresserror\\rey_umax_v9.txt", "w");
	FILE * stressfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\stresserror\\shearstress_v9.txt", "w");
	FILE * sttensfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\stresserror\\stresstensor_v9.txt", "w");
	FILE * forcefile = fopen("C:\\Users\\lukas\\documents\\MATLAB\\Erosion_LBM\\stresserror\\force_v9.txt", "w");
	FILE * nhatfile = fopen("C:\\Users\\lukas\\documents\\MATLAB\\Erosion_LBM\\stresserror\\nhat_v9.txt", "w");
	FILE * torfile = fopen("C:\\Users\\lukas\\documents\\MATLAB\\Erosion_LBM\\stresserror\\torfile_v9.txt", "w");
	FILE * erodefile = fopen("C:\\Users\\lukas\\documents\\MATLAB\\Erosion_LBM\\stresserror\\erodefile_v9.txt", "w");
	FILE * errorfile = fopen("C:\\Users\\lukas\\documents\\MATLAB\\Erosion_LBM\\stresserror\\errorfile_v9.txt", "w");
	int obchoice;
	int tend;
	int printi;
	double umax = 0;
	double uav = 0;

	Grid grid;
	momentum_direction e;
	EDF feq;
	velocity u;
	density rho; // IC: rho = 1 on all fluid nodes.
	direction_density f;
	direction_density ftemp;
	Stresstensor stresstensor; //nearest neighbour to surface solid nodes.
	Normalvector nhat;
	Wall_force tau_stress;
	vector3Ncubed F_D;
	vector3Ncubed F_sum;
	vector3Ncubed torque;
	density masschange;
	vectorNcubed ero_reso_check; //erosion resolution check. It's used to check if we erode faster than Delta_T or not.
	vectorNcubed F_vdw;
	masschange.clear();
	//==========================================================================================================
	// Choosing object to put in flow and applying IC.
	grid.printgrid(gridfile);
	cout << "xmax = " << (Nx - 1)*latspace << ", ymax = " << (Ny - 1)*latspace << ", zmax = " << (Nz - 1)*latspace << "\n";

	cout << "Choose object. 1: sphere. 2: cylinder pipe. 3: square pipe. 4: cylinder with sphere inside. 5: cylinder. 6: open box. 7: triangle-like. \n";
	//cin >> obchoice;
	obchoice = 1;
	//for (double rotation = -0.25; rotation <  0.26; rotation += 0.05) {
	Solid_list solid_list(obchoice, grid, pi * 0, e, parfile);
	updatePBC_solid(solid_list);
	solid_list.printsolid_list(solfile);

	// IC for rho (rho = 1 everywhere).
	//	for (int ix = 0; ix < Nx; ix++) {
	//		for (int iy = 0; iy < Ny; iy++) {
	//			for (int iz = 0; iz < Nz; iz++) {
	//				rho(ix, iy, iz) = 1.;
	//			}
	//		}
	//	}
	// Kolla att figuren som printades ut i solid_list roterar korrekt. Sen kör en lång simulering med olika vinklar.
	cout << "Grid and object created. Running simulation for t = \n";
	//cin >> tend;
	tend = 2001;
	IC(e, f, ftemp, solid_list);
	updateBC(f, -1, Bvel, rho, e, u);
	updateBC(ftemp, -1, Bvel, rho, e, u);
	printi = 0;
	//==========================================================================================================
	// Main program.
	int i_er = 1;
	int i_Fvdw = 1;

	for (int t = 0; t < tend; t++) {
		cout << t << " ";
		stream(solid_list, f, ftemp, e);
		updateBC(ftemp, t, Bvel, rho, e, u);
		macrovariables(u, rho, solid_list, ftemp, e); 
//		umax = find_umax(u);
//		uav = find_uav(u);
//		fprintf(reyfile, "%e %e %e\n", umax, Bvel[0] * 16. / mu, uav);
		if (t == 100 * printi) {
			printstuff(velfile, densfile, parfile, reyfile, stressfile, forcefile, nhatfile, sttensfile, torfile, erodefile, t, u, rho, tau_stress, F_D, nhat, stresstensor, torque, masschange);
			printi++;
		}
		edf(solid_list, u, rho, feq, e, edfforcedir);
		collision(solid_list, f, ftemp, feq, e);
		updateBC(f, t, Bvel, rho, e, u);
		computestress(e, ftemp, f, feq, solid_list, stresstensor, nhat, tau_stress, F_sum, masschange, F_vdw, i_er, i_Fvdw, rho);
		if (i_Fvdw == i_er)
			i_Fvdw++;
//		computetorque(solid_list, tau_stress, torque);
		if (t == i_er*Delta_T) {
			erosion(solid_list, e, F_sum, rho, f, edfforcedir, solfile, masschange, nhat, ero_reso_check, F_vdw, errorfile);
			i_er++;
		}
	}
	cout << "\n Done! \n";
	solid_list.clear_list();
	//	};
	//==========================================================================================================


	//	getchar();
	fclose(velfile);
	fclose(densfile);
	fclose(gridfile);
	fclose(parfile);
	fclose(solfile);
	fclose(forcefile);
	fclose(stressfile);
	fclose(torfile);

	return 0;
}





void IC(const momentum_direction& e, direction_density& f, direction_density& ftemp, Solid_list& solid_list) {
	// variables
	// -----------------------------
	int a = 0;
	int ix = 0;		// x 
	int iy = 0;		// y
	int iz = 0;		// z
	int cellist[27] = { 3, 2, 3, 2, 1, 2, 3, 2, 3  ,  2, 1, 2, 1, 0, 1, 2, 1, 2  ,  3, 2, 3, 2, 1, 2, 3, 2, 3 };
	double weights[4];
	double w = 0;
	weights[0] = 0.296296296296296;//8. / 27.;	
	weights[1] = 0.074074074074074;//2. / 27.;
	weights[2] = 0.018518518518519;// 1. / 54.;
	weights[3] = 0.004629629629630;//1. / 216.;
								   //------------------------------
								   // function body
								   //		cout << "\n IC: ------------------------------------" << "\n";
	for (a = 0; a < 27; a++) {
		//		cout << "\n a = " << a << "\n";
		for (iz = 0; iz < Nz; iz++) {
			//				cout << "\n";
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					if (solid_list(ix, iy, iz) != 1) {
						//if (a == 13) {
						w = weights[cellist[a]];
						f(ix, iy, iz, a) = w;
						ftemp(ix, iy, iz, a) = w;
						//}
						//else
						//{
						//f(ix, iy, iz, a) = 0.;
						//ftemp(ix, iy, iz, a) = 0.;
						//}
					}
					else {
						f(ix, iy, iz, a) = 0.;
						ftemp(ix, iy, iz, a) = 0.;
					}
					//						cout << " " << f(ix, iy, iz, a) << " ";
				}
			}
		}
	}


	/*	for (a = 0; a < 27; a++) {
	cout << "\n a = " << a << "\n";
	for (iz = -1 ; iz < Nz+1; iz++) {
	cout << "\n";
	for (iy = -1; iy < Ny+1; iy++) {
	for (ix = -1; ix < Nx+1; ix++) {
	cout << " " << f(ix, iy, iz, a) << " ";
	}
	}
	}
	}*/
}


void macrovariables(velocity& u, density& rho, Solid_list& solid_list, direction_density& f, momentum_direction& e) {
	// variables
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;		// momentum direction
	double umax;
	// function body

	//	cout << "\n Macro: ------------------------------------ \n";
	for (a = 0; a < 27; a++) {   //Could maybe swap 27 for 3^3, where 3 is dimension.
								 //		cout << "\n a = " << a << "\n";
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					//					cout << " " << f(ix, iy, iz, a) << " ";
					if (a == 0) {
						u(ix, iy, iz, 0) = 0.;
						u(ix, iy, iz, 1) = 0.;
						u(ix, iy, iz, 2) = 0.;
						rho(ix, iy, iz) = 0.;
					}

					if (solid_list(ix, iy, iz) == -1) {
						rho(ix, iy, iz) += f(ix, iy, iz, a);
						u(ix, iy, iz, 0) += e(a, 0) * f(ix, iy, iz, a);
						u(ix, iy, iz, 1) += e(a, 1) * f(ix, iy, iz, a);
						u(ix, iy, iz, 2) += e(a, 2) * f(ix, iy, iz, a);
					}
					//					cout << "x: " << u(ix,iy,iz,0) << ", y: " << u(ix,iy,iz,1) << ", z: " << u(ix,iy,iz,2) << "\n";
					//					cout << " " << e(a,0) << " " << e(a,1) << " " << e(a,2);
					//					if (ix == 1 && iy == 1 && iz == 0)
					//						cout << " " << f(ix, iy, iz, a) << " ";
					//					if (u(ix, iy, iz, 0) > 5000 || u(ix, iy, iz, 1) > 5000 || u(ix, iy, iz, 2) > 5000)
					//						cout << "Warning: Very high velocity!! \n";
				}
			}
			//			cout << "\n";
		}
	}
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				if (solid_list(ix, iy, iz) == -1) {
					//					cout << " " << u(ix, iy, iz, 0) << " " << u(ix, iy, iz, 1) << " " << u(ix, iy, iz, 2) << " " << rho(ix, iy, iz) << "\n";
					if (rho(ix, iy, iz) == 0.) {
						cout << "Error: rho = 0\n";
						break;
					}
					u(ix, iy, iz, 0) /= rho(ix, iy, iz);
					u(ix, iy, iz, 1) /= rho(ix, iy, iz);
					u(ix, iy, iz, 2) /= rho(ix, iy, iz);
				}
			}
		}
	}


};
double find_umax(velocity& u) {
	int iz = 0;
	int iy = 0;
	int ix = 0;
	double umax[3] = { 0 };
	double umaxtemp[3] = { 0 };
	double ans[4];
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				umaxtemp[0] = u(ix, iy, iz, 0);
				umaxtemp[1] = u(ix, iy, iz, 1);
				umaxtemp[2] = u(ix, iy, iz, 2);
				if (abs(umaxtemp[0]) > umax[0]) {
					umax[0] = umaxtemp[0];

				}
				if (abs(umaxtemp[1]) > umax[1])
					umax[1] = umaxtemp[1];
				if (abs(umaxtemp[2]) > umax[2])
					umax[2] = umaxtemp[2];
			}
		}
	}
	//	ans[0] = sqrt(umax[0] * umax[0] + umax[1] * umax[1] + umax[2] * umax[2]);
	//	ans[1] = 
	return sqrt(umax[0] * umax[0] + umax[1] * umax[1] + umax[2] * umax[2]);

}
double find_uav(velocity& u) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	double uxsum = 0;
	double uysum = 0;
	double uzsum = 0;
	double uav;

	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				uxsum += u(ix, iy, iz, 0);
				uysum += u(ix, iy, iz, 1);
				uzsum += u(ix, iy, iz, 2);
			}
		}
	}
	uxsum = uxsum / ((double)Nx*(double)Ny*(double)Nz);
	uysum = uysum / ((double)Nx*(double)Ny*(double)Nz);
	uzsum = uzsum / ((double)Nx*(double)Ny*(double)Nz);
	uav = sqrt(uxsum*uxsum + uysum*uysum + uzsum*uzsum);
	return uav;
}


void edf(Solid_list& solid_list, velocity& u, density& rho, EDF& feq, momentum_direction& e, int forcedirection) {
	// variables
	double f1 = 3.;
	double f2 = 9.*0.5;
	double f3 = 3.*0.5;
	double w = 0.;
	double usq = 0.;
	double ueq[3];
	double csq = 0.;
	double dotprod = 0.;
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;	// Direction counter. Used for feq and w. 
				// function body
				// cout << "\n EDF: -------------------------------------- \n";
	for (a = 0; a < 27; a++) {
		//		cout << "\n a = " << a << ": \n \n";
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {

					if (solid_list(ix, iy, iz) == -1) {
						/*						ueq[0] = u(ix, iy, iz, 0) + (tau*gg) / rho(ix, iy, iz);	//Add subroutine for this;
						ueq[1] = u(ix, iy, iz, 1);
						ueq[2] = u(ix, iy, iz, 2);
						w = weights[cellist[a]];
						//							cout << "\n " << ux << " " << uy << " " << uz << "\n";
						usq = ueq[0] * ueq[0] + ueq[1] * ueq[1] + ueq[2] * ueq[2];
						csq = c*c;
						dotprod = edfvecdot(e, ueq, a);		// Testa så att beräkna feq i den nya funktionen funkar.
						feq(ix, iy, iz, a) = rho(ix, iy, iz) * w * (1. + f1*dotprod / (csq)+f2*dotprod*dotprod / (csq*csq) - f3*usq / csq);*/ //c, c^2, c^2
						ueq[0] = u(ix, iy, iz, 0);
						ueq[1] = u(ix, iy, iz, 1);
						ueq[2] = u(ix, iy, iz, 2);
						feq(ix, iy, iz, a) = calcfeq(rho, e, ueq, ix, iy, iz, a, forcedirection);
					}
					//					if (ix == 1 && iy == 1 && iz == 0)
					//						cout << " " << feq(ix, iy, iz, a) << " ";
				}
			}
			//			cout << "\n";
		}
		//			cout << "\n";
	}
}
double edfvecdot(momentum_direction& e, double ueq[3], int a) {
	double answer = 0;
	answer = (double)e(a, 0) * ueq[0] + (double)e(a, 1) * ueq[1] + (double)e(a, 2) * ueq[2];
	return answer;
};
double calcfeq(density& rho, momentum_direction& e, double ueq[3], int ix, int iy, int iz, int a, int forcedirection) {
	double weights[4] = { 0 };
	weights[0] = 0.296296296296296;// 8. / 27.;	
	weights[1] = 0.074074074074074;// 2. / 27.;
	weights[2] = 0.018518518518519;// 1. / 54.;
	weights[3] = 0.004629629629630;// 1. / 216.;
	double feqval;
	double f1 = 3.;
	double f2 = 9.*0.5;
	double f3 = 3.*0.5;
	double w = 0.;
	double usq = 0.;
	double csq = 0.;
	double dotprod = 0.;
	ueq = calcforce(ueq, rho, ix, iy, iz, forcedirection);
	w = weights[cellist[a]];
	//							cout << "\n " << ux << " " << uy << " " << uz << "\n";
	usq = ueq[0] * ueq[0] + ueq[1] * ueq[1] + ueq[2] * ueq[2];
	csq = c*c;
	dotprod = (double)e(a, 0) * ueq[0] + (double)e(a, 1) * ueq[1] + (double)e(a, 2) * ueq[2];
	feqval = rho(ix, iy, iz) * w * (1. + f1*dotprod / (csq)+f2*dotprod*dotprod / (csq*csq) - f3*usq / csq);	// c, c^2, c^2.
	return feqval;
}
double * calcforce(double ueq[3], density& rho, int ix, int iy, int iz, int forcedirection) { //forcedirection == "0" betyder ingen kraft verkar på systemet.
	if (forcedirection == 1)
		ueq[0] += F*cos(((double)iz / Lz)*pi);	//tau*gg / rho(ix, iy, iz);
	if (forcedirection == 2)
		ueq[1] += F*cos(((double)iz / Lz)*pi);  //tau*gg / rho(ix, iy, iz);
	if (forcedirection == 3)
		ueq[2] += F*cos(((double)iz / Lz)*pi);  //tau*gg / rho(ix, iy, iz);
	return ueq;
}


void stream(Solid_list& solid_list, direction_density& f, direction_density& ftemp, momentum_direction& e) {
	// variable defs
	int ix = 0;			// x
	int iy = 0;			// y
	int iz = 0;			// z
	int a = 0;			// direction counter
	int ashift = 0;
	// function body
	/*	cout << "\n Stream: ------------------------------------ \n";
	for (a = 0; a < 27; a++) {
	cout << "\n a = " << a << "\n";
	for (iz = -1; iz < Nz + 1; iz++) {
	//			cout << "\n";
	for (iy = -1; iy < Ny + 1; iy++) {
	for (ix = -1; ix < Nx + 1; ix++) {
	if (ix == 1 && iy == 0 && iz == 0)
	cout << " " << f(ix, iy, iz, a) << " ";
	}
	}
	}
	}*/

	for (a = 0; a < 27; a++) {
		ashift = 2 * (13 - a);
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					switch (solid_list(ix, iy, iz)) {
					case 1:	// interior solid node = Do nothing
						break;
					case 0:	// surface solid node = bounceback BC
						ftemp(ix, iy, iz, a) = f(ix + e(26 - a, 0), iy + e(26 - a, 1), iz + e(26 - a, 2), a);
						break;
					case -1: // fluid node = regular stream
						ftemp(ix, iy, iz, a) = f(ix + e(26 - a, 0), iy + e(26 - a, 1), iz + e(26 - a, 2), a);	//the e:s shifts ix,iy,iz to the node streaming into the node we're looking at. 
						break;
					}
				}
			}
		}
	}

	/*	cout << "---------------------------------------------------------------";
	for (a = 0; a < 27; a++) {
	cout << "\n a = " << a << "\n";
	for (iz = -1; iz < Nz + 1; iz++) {
	//			cout << "\n";
	for (iy = -1; iy < Ny + 1; iy++) {
	for (ix = -1; ix < Nx + 1; ix++) {
	if (ix == 1 && iy == 0 && iz == 0)
	cout << " " << ftemp(ix, iy, iz, a) << " ";
	}
	}
	}
	}*/
}
void bouncebackBC(direction_density& ftemp, direction_density& f, momentum_direction& e, int ix, int iy, int iz, int a, int ashift) {

	/*	int ixtemp = 0;
	int iytemp = 0;
	int iztemp = 0;
	ixtemp = ix + e(26 - a, 0);
	iytemp = iy + e(26 - a, 1);
	iztemp = iz + e(26 - a, 2);
	if (ixtemp < 0)
	ixtemp += Nx;
	if (ixtemp >(Nx - 1))
	ixtemp -= Nx;
	if (iytemp < 0)
	iytemp += Ny;
	if (iytemp >(Ny - 1))
	iytemp -= Ny;
	if (iztemp < 0)
	iztemp += Nz;
	if (iztemp >(Nz - 1))
	iztemp -= Nz;*/
	//	cout << ftemp(ixtemp, iytemp, iztemp, a + ashift) << " " << f(ixtemp, iytemp, iztemp, a) << "\n"; 
	ftemp(ix, iy, iz, a + ashift) = f(ix, iy, iz, a);
	//	cout << ftemp(ixtemp, iytemp, iztemp, a + ashift) << " " << f(ixtemp, iytemp, iztemp, a);
}


void collision(Solid_list& solid_list, direction_density& f, direction_density& ftemp, EDF& feq, momentum_direction& e) {
	int ix = 0;
	int iy = 0;
	int iz = 0;

	int ixtemp;
	int iytemp;
	int iztemp;
	int a = 0;
	int ashift = 0;
	//	cout << "\n Collision: ------------------------------------\n \n";
	/*	for (a = 0; a < 27; a++) {
	cout << "\n a = " << a << "\n";
	for (iz = 0; iz < Nz; iz++) {
	cout << "\n";
	for (iy = 0; iy < Ny; iy++) {
	for (ix = 0; ix < Nx; ix++) {
	if (ix == 1 && iy == 1 && iz == 0)
	cout << " " << ftemp(ix, iy, iz, a) << " ";
	}
	}
	}
	}*/

	for (a = 0; a < 27; a++) {
		ashift = 2 * (13 - a);
		//		cout << "\n a = " << a << "\n";
		for (iz = 0; iz < Nz; iz++) {
			//				cout << "\n";
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					switch (solid_list(ix, iy, iz)) {
					case 1:	// Interior node
						break;
					case 0:	// Surface node
						f(ix, iy, iz, a + ashift) = ftemp(ix, iy, iz, a);
						break;
					case -1: //	Fluid node
						f(ix, iy, iz, a) = ftemp(ix, iy, iz, a) - (ftemp(ix, iy, iz, a) - feq(ix, iy, iz, a)) / tau;
						break;
					}
				}
			}
			//			cout << "\n";
		}
	}

	/*	for (a = 0; a < 27; a++) {
	cout << "\n a = " << a << "\n";
	for (iz = 0; iz < Nz; iz++) {
	cout << "\n";
	for (iy = 0; iy < Ny; iy++) {
	for (ix = 0; ix < Nx; ix++) {
	if (ix == 1 && iy == 1 && iz == 0)
	cout << " " << f(ix, iy, iz, a) << " ";
	}
	}
	}
	}*/
};


void updateBC(direction_density& f, int t, double Bvel[3], density& rho, momentum_direction& e, velocity& u) {
	int a = 0;
	int ix = 0;
	int iy = 0;
	int iz = 0;
	// 1 - periodic BC. 2 - velocity BC.
	for (a = 0; a < 27; a++) {
		if (BCtype[5] == 1) {	//periodic BC
								//----------------------------------------------------------------------
								// upper boundary layer 
								// main surface
			iz = 0;
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, iz + Nz, a)] = f(ix, iy, iz, a);
				}
			}
			// edges
			// copying to left upper edge
			//			for (iy = 0; iy < Ny; iy++)
			//				f.element[f.index(-1, iy, iz + Nz, a)] = f(Nx - 1, iy, iz, a);
			// copying to right upper edge
			//			for (iy = 0; iy < Ny; iy++)
			//				f.element[f.index(Nx, iy, iz + Nz, a)] = f(0, iy, iz, a);

			// copying to forward upper edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, -1, iz + Nz, a)] = f(ix, Ny - 1, iz, a);

			// copying to back upper edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, Ny, iz + Nz, a)] = f(ix, 0, iz, a);

			// corners
			//			f.element[f.index(-1, -1, Nz, a)] = f(Nx - 1, Ny - 1, 0, a);		// upper left-forward
			//			f.element[f.index(Nx, -1, Nz, a)] = f(0, Ny - 1, 0, a);			// upper right-forward
			//			f.element[f.index(Nx, Ny, Nz, a)] = f(0, 0, 0, a);	// upper right-back
			//			f.element[f.index(-1, Ny, Nz, a)] = f(Nx - 1, 0, 0, a);		// upper left-back
		}
		if (BCtype[5] == 2 && t == -1) {	// velocity BC, outlet
											//----------------------------------------------------------------------
											// upper boundary layer
											// main surface
			iz = 0;
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, iz + Nz, a)] = f.element[f.index(ix, iy, iz + Nz - 1, a)];
				}
			}
			// edges
			// left upper edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(-1, iy, iz + Nz, a)] = f.element[f.index(0, iy, iz + Nz - 1, a)];

			// right upper edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(Nx, iy, iz + Nz, a)] = f.element[f.index(Nx - 1, iy, iz + Nz - 1, a)];

			// forward upper edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, -1, iz + Nz, a)] = f.element[f.index(ix, 0, iz + Nz - 1, a)];

			// back upper edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, Ny, iz + Nz, a)] = f.element[f.index(ix, Ny - 1, iz + Nz - 1, a)];

			// corners
			f.element[f.index(-1, -1, Nz, a)] = f.element[f.index(0, 0, Nz - 1, a)];			// upper left-forward
			f.element[f.index(Nx, -1, Nz, a)] = f.element[f.index(Nx - 1, 0, Nz - 1, a)];		// upper right-forward
			f.element[f.index(Nx, Ny, Nz, a)] = f.element[f.index(Nx - 1, Ny - 1, Nz - 1, a)];	// upper right-back
			f.element[f.index(-1, Ny, Nz, a)] = f.element[f.index(0, Ny - 1, Nz - 1, a)];		// upper left-back
		}
		//----------------------------------------------------------------------
		// lower boundary layer (copiyng from real upper layer)
		// main surface
		if (BCtype[4] == 1) {	// Periodic BC
			iz = Nz - 1;
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, -1, a)] = f(ix, iy, iz, a);
				}
			}
			// edges
			// copying to left lower edge
			//			for (iy = 0; iy < Ny; iy++)
			//				f.element[f.index(-1, iy, -1, a)] = f(Nx - 1, iy, iz, a);

			// copying to right lower edge
			//			for (iy = 0; iy < Ny; iy++)
			//				f.element[f.index(Nx, iy, -1, a)] = f(0, iy, iz, a);

			// copying to forward lower edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, -1, -1, a)] = f(ix, Ny - 1, iz, a);

			// copying to back lower edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, Ny, -1, a)] = f(ix, 0, iz, a);

			// corners
			//			f.element[f.index(-1, -1, -1, a)] = f(Nx - 1, Ny - 1, Nz - 1, a);		// bottom left-forward
			//			f.element[f.index(Nx, -1, -1, a)] = f(0, Ny - 1, Nz - 1, a);		// bottom right-forward
			//			f.element[f.index(Nx, Ny, -1, a)] = f(0, 0, Nz - 1, a);			// bottom right-back
			//			f.element[f.index(-1, Ny, -1, a)] = f(Nx - 1, 0, Nz - 1, a);		// bottom left-back
		}
		if (BCtype[4] == 2 && t == -1) {	// Velocity BC, inlet
			iz = Nz - 1;
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);
				}
			}
			// edges

			// copying to right lower edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(Nx, iy, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);

			// copying to forward lower edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, -1, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);

			// copying to back lower edge
			for (ix = 0; ix < Nx; ix++)
				f.element[f.index(ix, Ny, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);

			// corners
			f.element[f.index(-1, -1, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);		// bottom left-forward
			f.element[f.index(Nx, -1, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);		// bottom right-forward
			f.element[f.index(Nx, Ny, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);			// bottom right-back
			f.element[f.index(-1, Ny, -1, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);		// bottom left-back
		}
		//-----------------------------------------------------------------------------------------------------------
		// set right side ghost points
		if (BCtype[1] == 1) {	// Periodic BC
			ix = Nx;
			for (iz = 0; iz < Nz; iz++) {
				f(Nx, -1, iz, a) = f(0, Ny - 1, iz, a);			// copy to right-forward corner side
				f(Nx, Ny, iz, a) = f(0, 0, iz, a);				// copy to right-back corner side
				for (iy = 0; iy < Ny; iy++) {
					f.element[f.index(ix, iy, iz, a)] = f(0, iy, iz, a);
				}
			}
		}
		if (BCtype[1] == 2 && t == -1) {	// Velocity BC, outlet
			ix = Nx;
			for (iz = 0; iz < Nz; iz++) {
				f(Nx, -1, iz, a) = f(Nx - 1, 0, iz, a);			// copy to right-forward edge
				f(Nx, Ny, iz, a) = f(Nx - 1, Ny - 1, iz, a);	// copy to right-back edge
				for (iy = 0; iy < Ny; iy++) {
					f.element[f.index(ix, iy, iz, a)] = f.element[f.index(ix - 1, iy, iz, a)];	// Copy f:s in outlet to ghost points from real points.
				}
			}

			// copying to right upper edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(Nx, iy, Nz, a)] = f(Nx - 1, iy, Nz - 1, a);
			// copying to right lower edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(Nx, iy, -1, a)] = f(Nx - 1, iy, 0, a);

			f.element[f.index(Nx, -1, Nz, a)] = f(Nx - 1, 0, Nz - 1, a);			// upper right-forward
			f.element[f.index(Nx, Ny, Nz, a)] = f(Nx - 1, Ny - 1, Nz - 1, a);		// upper right-back
			f.element[f.index(Nx, -1, -1, a)] = f(Nx - 1, 0, 0, a);					// bottom right-forward
			f.element[f.index(Nx, Ny, -1, a)] = f(Nx - 1, Ny - 1, 0, a);			// bottom right-back
		}
		//-----------------------------------------------------------------------------------------------------------
		// set left side ghost points
		if (BCtype[0] == 1) {	// Periodic BC
			ix = -1;
			for (iz = 0; iz < Nz; iz++) {
				f(-1, -1, iz, a) = f(Nx - 1, Ny - 1, iz, a);	// set left-forward corner side
				f(-1, Ny, iz, a) = f(Nx - 1, 0, iz, a);			// set left-back corner side
				for (iy = 0; iy < Ny; iy++) {
					f.element[f.index(ix, iy, iz, a)] = f(ix + Nx, iy, iz, a);
				}
			}
		}
		if (BCtype[0] == 2 && t == -1) {	// Velocity BC, inlet
			ix = -1;
			for (iz = 0; iz < Nz; iz++) {
				f(-1, -1, iz, a) = calcfeq(rho, e, Bvel, 0, 0, 0, a, 0);			// set left-forward edge 
				f(-1, Ny, iz, a) = calcfeq(rho, e, Bvel, 0, 0, 0, a, 0);			// set left-back edge
				for (iy = 0; iy < Ny; iy++) {
					f.element[f.index(ix, iy, iz, a)] = calcfeq(rho, e, Bvel, 0, 0, 0, a, 0);
				}
			}

			// copying to left upper edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(-1, iy, Nz, a)] = calcfeq(rho, e, Bvel, Nx - 1, iy, 0, a, 0);
			// copying to left lower edge
			for (iy = 0; iy < Ny; iy++)
				f.element[f.index(-1, iy, -1, a)] = calcfeq(rho, e, Bvel, 0, 0, 0, a, 0);

			f.element[f.index(-1, -1, Nz, a)] = calcfeq(rho, e, Bvel, Nx - 1, Ny - 1, 0, a, 0);		// upper left-forward
			f.element[f.index(-1, Ny, Nz, a)] = calcfeq(rho, e, Bvel, Nx - 1, 0, 0, a, 0);				// upper left-back
			f.element[f.index(-1, -1, -1, a)] = calcfeq(rho, e, Bvel, Nx - 1, Ny - 1, Nz - 1, a, 0);	// bottom left-forward
			f.element[f.index(-1, Ny, -1, a)] = calcfeq(rho, e, Bvel, Nx - 1, 0, Nz - 1, a, 0);		// bottom left-back

		}
		//-----------------------------------------------------------------------------------------------------------
		// set back side ghost points
		if (BCtype[3] == 1) {	// Periodic BC
			iy = Ny;
			for (iz = 0; iz < Nz; iz++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, iz, a)] = f(ix, 0, iz, a);
				}
			}
		}
		if (BCtype[3] == 2 && t == -1) {	// Velocity BC, outlet
			iy = Ny;
			for (iz = 0; iz < Nz; iz++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, iz, a)] = f.element[f.index(ix, iy - 1, iz, a)];
				}
			}
		}
		// set front side ghost points
		if (BCtype[2] == 1) {	// Periodic BC
			iy = -1;
			for (iz = 0; iz < Nz; iz++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, iz, a)] = f(ix, iy + Ny, iz, a);
				}
			}
		}
		if (BCtype[2] == 2 && t == -1) {	// Velocity BC, inlet
			iy = -1;
			for (iz = 0; iz < Nz; iz++) {
				for (ix = 0; ix < Nx; ix++) {
					f.element[f.index(ix, iy, iz, a)] = calcfeq(rho, e, Bvel, 1, 1, 1, a, 0);
				}
			}
		}

	}
}

void updatePBC(direction_density& f) {
	int a = 0;
	int ix = 0;
	int iy = 0;
	int iz = 0;
	//	switch (choice) {
	//	case 1:
	for (a = 0; a < 27; a++) {
		//----------------------------------------------------------------------
		// upper boundary layer (copying from real bottom layer)
		// main surface
		iz = 0;
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				f.element[f.index(ix, iy, iz + Nz, a)] = f(ix, iy, iz, a);
			}
		}
		// edges
		// copying to left upper edge
		for (iy = 0; iy < Ny; iy++)
			f.element[f.index(-1, iy, iz + Nz, a)] = f(Nx - 1, iy, iz, a);

		// copying to right upper edge
		for (iy = 0; iy < Ny; iy++)
			f.element[f.index(Nx, iy, iz + Nz, a)] = f(0, iy, iz, a);

		// copying to forward upper edge
		for (ix = 0; ix < Nx; ix++)
			f.element[f.index(ix, -1, iz + Nz, a)] = f(ix, Ny - 1, iz, a);

		// copying to back upper edge
		for (ix = 0; ix < Nx; ix++)
			f.element[f.index(ix, Ny, iz + Nz, a)] = f(ix, 0, iz, a);

		// corners
		f.element[f.index(-1, -1, Nz, a)] = f(Nx - 1, Ny - 1, 0, a);		// upper left-forward
		f.element[f.index(Nx, -1, Nz, a)] = f(0, Ny - 1, 0, a);			// upper right-forward
		f.element[f.index(Nx, Ny, Nz, a)] = f(0, 0, 0, a);	// upper right-back
		f.element[f.index(-1, Ny, Nz, a)] = f(Nx - 1, 0, 0, a);		// upper left-back
																	//----------------------------------------------------------------------
																	// lower boundary layer (copiyng from real upper layer)
																	// main surface
		iz = Nz - 1;
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				f.element[f.index(ix, iy, -1, a)] = f(ix, iy, iz, a);
			}
		}
		// edges
		// copying to left lower edge
		for (iy = 0; iy < Ny; iy++)
			f.element[f.index(-1, iy, -1, a)] = f(Nx - 1, iy, iz, a);

		// copying to right lower edge
		for (iy = 0; iy < Ny; iy++)
			f.element[f.index(Nx, iy, -1, a)] = f(0, iy, iz, a);

		// copying to forward lower edge
		for (ix = 0; ix < Nx; ix++)
			f.element[f.index(ix, -1, -1, a)] = f(ix, Ny - 1, iz, a);

		// copying to back lower edge
		for (ix = 0; ix < Nx; ix++)
			f.element[f.index(ix, Ny, -1, a)] = f(ix, 0, iz, a);

		// corners
		f.element[f.index(-1, -1, -1, a)] = f(Nx - 1, Ny - 1, Nz - 1, a);		// bottom left-forward
		f.element[f.index(Nx, -1, -1, a)] = f(0, Ny - 1, Nz - 1, a);		// bottom right-forward
		f.element[f.index(Nx, Ny, -1, a)] = f(0, 0, Nz - 1, a);			// bottom right-back
		f.element[f.index(-1, Ny, -1, a)] = f(Nx - 1, 0, Nz - 1, a);		// bottom left-ba
																			//----------------------------------------------------------------------
																			// copy to right side
		ix = Nx;
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				f.element[f.index(ix, iy, iz, a)] = f(0, iy, iz, a);
			}
		}
		// copy to left side
		ix = -1;
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				f.element[f.index(ix, iy, iz, a)] = f(ix + Nx, iy, iz, a);
			}
		}
		// copy to back side
		iy = Ny;
		for (iz = 0; iz < Nz; iz++) {
			for (ix = 0; ix < Nx; ix++) {
				f.element[f.index(ix, iy, iz, a)] = f(ix, 0, iz, a);
			}
		}
		// copy to forward side
		iy = -1;
		for (iz = 0; iz < Nz; iz++) {
			for (ix = 0; ix < Nx; ix++) {
				f.element[f.index(ix, iy, iz, a)] = f(ix, iy + Ny, iz, a);
			}
		}
		for (iz = 0; iz < Nz; iz++) {
			f(-1, -1, iz, a) = f(Nx - 1, Ny - 1, iz, a);	// copy to left-forward corner side
			f(Nx, -1, iz, a) = f(0, Ny - 1, iz, a);			// copy to right-forward corner side
			f(-1, Ny, iz, a) = f(Nx - 1, 0, iz, a);			// copy to left-back corner side
			f(Nx, Ny, iz, a) = f(0, 0, iz, a);				// copy to right-back corner side
		}
	}
}

void updatePBC_solid(Solid_list& solid_list) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	//----------------------------------------------------------------------
	// upper boundary layer (copying from real bottom layer)
	// main surface
	iz = 0;
	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < Nx; ix++) {
			solid_list(ix, iy, iz + Nz) = solid_list(ix, iy, iz);
		}
	}
	// edges
	// copying to left upper edge
	for (iy = 0; iy < Ny; iy++)
		solid_list(-1, iy, iz + Nz) = solid_list(Nx - 1, iy, iz);

	// copying to right upper edge
	for (iy = 0; iy < Ny; iy++)
		solid_list(Nx, iy, iz + Nz) = solid_list(0, iy, iz);

	// copying to forward upper edge
	for (ix = 0; ix < Nx; ix++)
		solid_list(ix, -1, iz + Nz) = solid_list(ix, Ny - 1, iz);

	// copying to back upper edge
	for (ix = 0; ix < Nx; ix++)
		solid_list(ix, Ny, iz + Nz) = solid_list(ix, 0, iz);

	// corners
	solid_list(-1, -1, Nz) = solid_list(Nx - 1, Ny - 1, 0);		// upper left-forward
	solid_list(Nx, -1, Nz) = solid_list(0, Ny - 1, 0);			// upper right-forward
	solid_list(Nx, Ny, Nz) = solid_list(0, 0, 0);	// upper right-back
	solid_list(-1, Ny, Nz) = solid_list(Nx - 1, 0, 0);		// upper left-back
															//----------------------------------------------------------------------
															// lower boundary layer (copiyng from real upper layer)
															// main surface
	iz = Nz - 1;
	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < Nx; ix++) {
			solid_list(ix, iy, -1) = solid_list(ix, iy, iz);
		}
	}
	// edges
	// copying to left lower edge
	for (iy = 0; iy < Ny; iy++)
		solid_list(-1, iy, -1) = solid_list(Nx - 1, iy, iz);

	// copying to right lower edge
	for (iy = 0; iy < Ny; iy++)
		solid_list(Nx, iy, -1) = solid_list(0, iy, iz);

	// copying to forward lower edge
	for (ix = 0; ix < Nx; ix++)
		solid_list(ix, -1, -1) = solid_list(ix, Ny - 1, iz);

	// copying to back lower edge
	for (ix = 0; ix < Nx; ix++)
		solid_list(ix, Ny, -1) = solid_list(ix, 0, iz);

	// corners
	solid_list(-1, -1, -1) = solid_list(Nx - 1, Ny - 1, Nz - 1);		// bottom left-forward
	solid_list(Nx, -1, -1) = solid_list(0, Ny - 1, Nz - 1);		// bottom right-forward
	solid_list(Nx, Ny, -1) = solid_list(0, 0, Nz - 1);			// bottom right-back
	solid_list(-1, Ny, -1) = solid_list(Nx - 1, 0, Nz - 1);		// bottom left-ba
																//----------------------------------------------------------------------
																// copy left side
	ix = Nx;
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			solid_list(ix, iy, iz) = solid_list(0, iy, iz);
		}
	}
	// copy right side
	ix = -1;
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			solid_list(ix, iy, iz) = solid_list(ix + Nx, iy, iz);
		}
	}
	// copy forward side
	iy = Ny;
	for (iz = 0; iz < Nz; iz++) {
		for (ix = 0; ix < Nx; ix++) {
			solid_list(ix, iy, iz) = solid_list(ix, 0, iz);
		}
	}
	// copy back side
	iy = -1;
	for (iz = 0; iz < Nz; iz++) {
		for (ix = 0; ix < Nx; ix++) {
			solid_list(ix, iy, iz) = solid_list(ix, iy + Ny, iz);
		}
	}
	for (iz = 0; iz < Nz; iz++) {
		solid_list(-1, -1, iz) = solid_list(Nx - 1, Ny - 1, iz);	// copy left-forward corner side
		solid_list(Nx, -1, iz) = solid_list(0, Ny - 1, iz);			// copy right-forward corner side
		solid_list(-1, Ny, iz) = solid_list(Nx - 1, 0, iz);			// copy left-back corner side
		solid_list(Nx, Ny, iz) = solid_list(0, 0, iz);				// copy right-back corner side
	}
}

void printstuff(FILE * velfile, FILE * densfile, FILE * parfile, FILE * reyfile, FILE * stressfile, FILE * forcefile, FILE * nhatfile, FILE * sttensfile, FILE * torfile, FILE * erodefile, int t, velocity& u, density& rho, Wall_force& tau_stress, vector3Ncubed& F_D, Normalvector& nhat, Stresstensor& stresstensor, vector3Ncubed& torque, density& erodelist) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	double F[3] = { 0 };
	if (t == 0) {
		cout << "\n tau = " << tau << "\n";
		cout << "\n Lz = " << Lz << "\n";
		cout << "\n gg = " << gg << "\n";
		cout << "\n F = " << F << "\n";
		fprintf(parfile, "%e %e %e %e %e %i %i %i ", Re, tau, Lz, umax_theo, gg, Nx, Ny, Nz);
	}
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				fprintf(velfile, "%e %e %e ", u(ix, iy, iz, 0), u(ix, iy, iz, 1), u(ix, iy, iz, 2));
				fprintf(densfile, "%e ", rho(ix, iy, iz));
				fprintf(stressfile, "%e %e %e ", tau_stress(ix, iy, iz, 0), tau_stress(ix, iy, iz, 1), tau_stress(ix, iy, iz, 2));
				fprintf(forcefile, "%e %e %e ", F_D(ix, iy, iz, 0), F_D(ix, iy, iz, 1), F_D(ix, iy, iz, 2));
				fprintf(nhatfile, "%e %e %e ", nhat(ix, iy, iz, 0), nhat(ix, iy, iz, 1), nhat(ix, iy, iz, 2));
				fprintf(torfile, "%e %e %e ", torque(ix, iy, iz, 0), torque(ix, iy, iz, 1), torque(ix, iy, iz, 2));
				fprintf(erodefile, "%i ", erodelist(ix, iy, iz));
				//F[0] += F_D(ix, iy, iz, 0);
				//F[1] += F_D(ix, iy, iz, 1);
				//F[2] += F_D(ix, iy, iz, 2);
			}
		}
	}
	if (t == 3000) {
		for (int a = 0; a < 3; a++) {
			for (int b = 0; b < 3; b++) {

				for (iz = 0; iz < Nz; iz++) {
					for (iy = 0; iy < Ny; iy++) {
						for (ix = 0; ix < Nx; ix++) {
							fprintf(sttensfile, "%e ", stresstensor(ix, iy, iz, a, b));
						}
					}
				}
				fprintf(sttensfile, "\n");
			}
		}
	}
	fprintf(forcefile, "\n");
	//fprintf(forcefile, "%e %e %e\n", F[0], F[1], F[2]);
	fprintf(velfile, "\n");
	fprintf(densfile, "\n");
	fprintf(stressfile, "\n");
	fprintf(nhatfile, "\n");
	fprintf(torfile, "\n");
	fprintf(erodefile, "\n");
}

vector<int> readBC(string x0BC, string xNBC, string y0BC, string yNBC, string z0BC, string zNBC) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;
	int i = 0;
	vector<string> BC = { "periodicBC", "velocityBC" };	// if an additional boundary is implemented, loop the for loop below iend + 1 times.
	vector<int> BCtype;
	BCtype.resize(6);
	fill(BCtype.begin(), BCtype.end(), 0);

	for (i = 0; i < 2; i++) {
		if (x0BC.compare(BC[i]) == 0)
			BCtype[0] = 1 + i;
		if (xNBC.compare(BC[i]) == 0)
			BCtype[1] = 1 + i;
		if (y0BC.compare(BC[i]) == 0)
			BCtype[2] = 1 + i;
		if (yNBC.compare(BC[i]) == 0)
			BCtype[3] = 1 + i;
		if (z0BC.compare(BC[i]) == 0)
			BCtype[4] = 1 + i;
		if (zNBC.compare(BC[i]) == 0)
			BCtype[5] = 1 + i;
	}
	for (i = 0; i < 6; i++) {
		if (BCtype[i] == 0)
			cout << "Error in BC. Edge " << i << " has no BC. Assign a BC to this edge.";
	}

	return BCtype;



}

void computeVDWforce(int ix, int iy, int iz, momentum_direction& e, Solid_list& solid_list, vectorNcubed& F_vdw) {
	double WDWsum = 0.;
	int ixshift = 0;
	int iyshift = 0;
	int izshift = 0;
	for (int a = 0; a < 27; a++) {
		ixshift = ix + e(a, 0);
		iyshift = iy + e(a, 1);
		izshift = iz + e(a, 2);
		if (solid_list(ixshift, iyshift, izshift) == 1 || solid_list(ixshift, iyshift, izshift) == 0) // add wdwforce from all solid nodes.
			WDWsum = WDWsum + VDWforce(e(a, 0), e(a, 1), e(a, 2));
	}
	F_vdw(ix, iy, iz) = WDWsum;
}
void computestress(momentum_direction& e, direction_density& ftemp, direction_density& f, EDF& feq, Solid_list& solid_list, Stresstensor& stresstensor, Normalvector& nhat, Wall_force& tau_stress, vector3Ncubed& F_sum, density& masschange, vectorNcubed& F_vdw, int i_er, int i_Fvdw, density& rho) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;
	int i = 0;	//e component index. Row index. i:th component of the stress.
	int j = 0;	//e component index. Column index. i:th component of the stress acting on the j:th surface.
	int ixshift = 0; //temp nhat x
	int iyshift = 0; //temp nhat y
	int izshift = 0; //temp nhat z
	int normcount = 0;
	double normvec[3] = { 0. };
	double normfactvec[3] = { sq3inv, sq2inv, 1. };
	double normfactor = 0;
	int Nf = 0; //number of nnn fluid points
	double WDWsum = 0.;
	double FF = 0;
	
	// Calculating stress tensor, normal vectors, forces and erosion.
	stresstensor.clear();
	for (a = 0; a < 27; a++) {
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {

					if (solid_list(ix, iy, iz) == -1) { // Fluid node
						for (i = 0; i < 3; i++) {
							for (j = 0; j < 3; j++) {
								stresstensor(ix, iy, iz, i, j) += (1. - (1. / (2.*tau)))*(ftemp(ix, iy, iz, a) - feq(ix, iy, iz, a))*e(a, i)*e(a, j);
							}
						}
					}
					if (solid_list(ix, iy, iz) == 0) { // Surface solid node
						if (solid_list(ix + e(a, 0), iy + e(a, 1), iz + e(a, 2)) != -1) {
							nhat(ix, iy, iz, 0) += e(26 - a, 0);
							nhat(ix, iy, iz, 1) += e(26 - a, 1);
							nhat(ix, iy, iz, 2) += e(26 - a, 2);
						}
					}

				}
			}
		}
	}
	// Normalizing normal vector. We only care about direction and that abs(nhat) = 1. 1/sqrt(3)*(1,1,1) = 1/sqrt(12)*(2,2,2).
	// So we set nhats elements = 1, 0, -1 so we can use them as indeces too. nhat essentially becomes a normalized e vector (mom_dir).
	// This works for our appprox surface. For a real sphere ex. then this is wrong. Then we need nhat = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta)).
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				if (solid_list(ix, iy, iz) == 0) { //Surface solid point
					normcount = 0;
					for (i = 0; i < 3; i++) {
						if (nhat(ix, iy, iz, i) > 0)
							nhat(ix, iy, iz, i) = 1;
						if (nhat(ix, iy, iz, i) < 0)
							nhat(ix, iy, iz, i) = -1;
						if (nhat(ix, iy, iz, i) == 0)
							normcount++;
					}
					if (normcount == 3) { // This corresponds to an isolated point. In theory this would just float away with the fluid so we erode away that point. Might be thin ice here.
						cout << "\n Isolated point! normal vector = (0, 0, 0).\n";
					}
					normfactor = normvec[normcount];
					/*					ixtemp = nhat(ix, iy, iz, 0);
					iytemp = nhat(ix, iy, iz, 1);
					iztemp = nhat(ix, iy, iz, 2);
					normfactor = sqrt(ixtemp*ixtemp + iytemp*iytemp + iztemp*iztemp);
					nhat(ix, iy, iz, 0) /= normfactor;
					nhat(ix, iy, iz, 1) /= normfactor;
					nhat(ix, iy, iz, 2) /= normfactor;*/
				}

			}
		}
	}
	// Computing wall shear stress and drag force (actually the total force, but for now its just drag).
	//F_D.clear();
	tau_stress.clear();
	//	for (a = 0; a < 27; a++) {
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {

				if (solid_list(ix, iy, iz) == 0) { // solid surface point
					// Initialize erodelist with 0:s. This won't work otherwise.
					for (i = 0; i < 3; i++) {
						for (j = 0; j < 3; j++) {
							tau_stress(ix, iy, iz, i) += -nhat(ix, iy, iz, j)*stresstensor(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2), i, j); //Should be normfactor * surface area exposed to the fluid. But these 2 cancel out, so no contribution from them.
							//F_sum(ix, iy, iz, i) += tau_stress(ix, iy, iz, i);
							//tau_stress(ix, iy, iz, i) += -e(a,j)*stresstensor(ix + e(a, 0), iy + e(a, 1), iz + e(a, 2), i, j);
						}//
					}
					if (i_Fvdw == i_er)
						computeVDWforce(ix, iy, iz, e, solid_list, F_vdw);

					FF = sqrt(pow(tau_stress(ix, iy, iz, 0), 2) + pow(tau_stress(ix, iy, iz, 1), 2) + pow(tau_stress(ix, iy, iz, 2), 2));
					if (FF > F_vdw(ix,iy,iz)) { //if the fluid force is greater than the WDW force from all solid nodes.
						masschange(ix, iy, iz) += -kappa_er*FF*dt; // first part of mass change eq. Second part is in function "erosion".
					}
					else
						masschange(ix, iy, iz) += 0; //don't erode point
				}



			}
		}
	}


}
double VDWforce(int x, int y, int z) {
	// VDW force scaling with 1/r^2.
	if (x != 0 && y != 0 && z != 0)
		return 0.33333333333333*VDW_0; // r = sqrt(3). 1/3 ~= 0.33333... 
	if ((x != 0 && y != 0) || (x != 0 && z != 0) || (y != 0 && z != 0))
		return 0.5*VDW_0; // r = sqrt(2) => 1/2
	if (x != 0 || y != 0 || z != 0)
		return VDW_0; // r = 1.
	//if none of the ifs trigger, x = y = z = 0.
	return 0; // the solid node itself.
}

void computetorque(Solid_list& solid_list, Wall_force& tau_stress, vector3Ncubed& torque) {
	double r[3];
	int ix = 0;
	int iy = 0;
	int iz = 0;
	dvec result;
	result.resize(3);
	torque.clear();
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				if (solid_list(ix, iy, iz) == 0) { // surface node
					r[0] = double(ix) - solid_list.center[0];
					r[1] = double(iy) - solid_list.center[1];
					r[2] = double(iz) - solid_list.center[2];
					result = crossproduct(tau_stress, r, ix, iy, iz);
					torque(ix, iy, iz, 0) = result[0];
					torque(ix, iy, iz, 1) = result[1];
					torque(ix, iy, iz, 2) = result[2];
				}
			}
		}
	}


}
dvec crossproduct(Wall_force& tau_stress, double r[3], int ix, int iy, int iz) {
	// assumed vector of dimension 3.
	dvec result;
	result.resize(3);

	result[0] = r[1] * tau_stress(ix, iy, iz, 2) - r[2] * tau_stress(ix, iy, iz, 1);
	result[1] = -r[0] * tau_stress(ix, iy, iz, 2) + r[2] * tau_stress(ix, iy, iz, 0);
	result[2] = r[0] * tau_stress(ix, iy, iz, 1) - r[1] * tau_stress(ix, iy, iz, 0);

	return result;
}

// Function should be after collosion.
void erosion(Solid_list& solid_list, momentum_direction& e, vector3Ncubed& F_sum, density& rho, direction_density& f, int forcedirection, FILE * solfile, density& masschange, Normalvector& nhat, vectorNcubed& ero_reso_check, vectorNcubed& F_vdw, FILE * errorfile) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int ixshift = 0;
	int iyshift = 0;
	int izshift = 0;
	int a = 0;
	double FFsq = 0.; //Fluid Force squared.
	double ueq[3] = { 0 };
	double weights[4] = { 0 };
	weights[0] = 0.296296296296296;// 8. / 27.;	
	weights[1] = 0.074074074074074;// 2. / 27.;
	weights[2] = 0.018518518518519;// 1. / 54.;
	weights[3] = 0.004629629629630;// 1. / 216.;
								   //double Fmaxsq = F_sum.maxsq();

/*
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				erodelist(ix, iy, iz) = 0.;
				if (solid_list(ix, iy, iz) == 0) { //if we're standing on a surface node
					WDWsqsum = 0.;
					for (a = 0; a < 27; a++) {
						ixshift = ix + e(a, 0);
						iyshift = iy + e(a, 1);
						izshift = iz + e(a, 2);
						if (solid_list(ixshift, iyshift, izshift) == 1 || solid_list(ixshift, iyshift, izshift) == 0) // add wdwforce from all solid nodes.
							WDWsqsum = WDWsqsum + WDWforce;
					}
					FFsq = (pow(F_sum(ix, iy, iz, 0) / (double)updatefreq, 2) + pow(F_sum(ix, iy, iz, 1) / (double)updatefreq, 2) + pow(F_sum(ix, iy, iz, 2) / (double)updatefreq, 2));
					if (FFsq > WDWsqsum) { //if the fluid force is greater than the WDW force from all solid nodes.
						erodelist(ix, iy, iz) = (-kappa_er*sqrt(FFsq - WDWsqsum)); //erode point. How to choose limit for erosion of a point? (m_star). Need to add code for F_sum. Only things that contribute should go into the sum.
					}
					else
						erodelist(ix, iy, iz) = 0; //don't erode point
				}
			}
		}
	}
	*/

	//erodes all points that's lost all their mass
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {

				masschange(ix, iy, iz) += kappa_er*Delta_T*F_vdw(ix, iy, iz);
				if (solid_list(ix, iy, iz) == 0) { // if we're at a surface point
					ero_reso_check(ix, iy, iz) += 1.; 
					if (-masschange(ix, iy, iz) > masspernode) { // if more mass has eroded away than existing mass per point
						if (ero_reso_check(ix, iy, iz) == 1.) { // if a point eroded away during one Delta_T window => too low resolution.
							cout << "Resolution too low! Increase Delta_T, mass per node or Kappa.\n";
							fprintf(errorfile, "%i", 1);
						}
						ero_reso_check(ix, iy, iz) = 0.;
						erodepoint(ix, iy, iz, masschange, solid_list, rho, nhat, f, e);
						/*masschange(ix, iy, iz) = 0.;
						solid_list(ix, iy, iz) = -1; //surface node becomes fluid node.
						rho(ix, iy, iz) = rho(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2)); //fluid node is initialized with same density as the interface node. This could prob be extrapolated.
																															 // check that rho != 0 at new points.
						for (a = 0; a < 27; a++) { //Initialize new fluid point with same f as interface node. Also, check all nearby nodes. If it's a interior solid node, it becomes a surface.
							f(ix, iy, iz, a) = f(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2), a);
							ixshift = ix + e(a, 0);
							iyshift = iy + e(a, 1);
							izshift = iz + e(a, 2);
							if (solid_list(ixshift, iyshift, izshift) == 1) {
								solid_list(ixshift, iyshift, izshift) = 0;
								f(ixshift, iyshift, izshift, a) = weights[cellist[a]];
								rho(ixshift, iyshift, izshift) = 1.;
							}
						}*/
					}
				}

			}
		}
	}
	
	//erodes all points that become isolated
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				if (solid_list(ix, iy, iz) == 0) {	// surface node
					if (IsForeverAlone(ix, iy, iz, solid_list, e) == 1) // if surface node is isolated
						erodepoint(ix, iy, iz, masschange, solid_list, rho, nhat, f, e);
				}
			}
		}
	}



	solid_list.printsolid_list(solfile);
	//masschange.clear();
	//F_sum.clear();
}
void erodepoint(int ix, int iy, int iz, density& masschange, Solid_list& solid_list, density& rho, Normalvector& nhat, direction_density& f, momentum_direction& e) {
	double weights[4] = { 0 };
	weights[0] = 0.296296296296296;// 8. / 27.;	
	weights[1] = 0.074074074074074;// 2. / 27.;
	weights[2] = 0.018518518518519;// 1. / 54.;
	weights[3] = 0.004629629629630;// 1. / 216.;
	int ixshift = 0;
	int iyshift = 0;
	int izshift = 0;
	masschange(ix, iy, iz) = 0.;
	solid_list(ix, iy, iz) = -1; //surface node becomes fluid node.
	rho(ix, iy, iz) = rho(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2)); // fluid node is initialized with same density as the interface node. This could prob be extrapolated.
																										 // check that rho != 0 at new points.
	for (int a = 0; a < 27; a++) { //Initialize new fluid point with same f as interface node. Also, check all nearby nodes. If it's a interior solid node, it becomes a surface.
		f(ix, iy, iz, a) = f(ix + nhat(ix, iy, iz, 0), iy + nhat(ix, iy, iz, 1), iz + nhat(ix, iy, iz, 2), a); // same f as interface node.
		ixshift = ix + e(a, 0);
		iyshift = iy + e(a, 1);
		izshift = iz + e(a, 2);
		if (solid_list(ixshift, iyshift, izshift) == 1) { // if solid node, turn it to a surface node, initialize it with IC f, set rho = 1.
			solid_list(ixshift, iyshift, izshift) = 0;
			f(ixshift, iyshift, izshift, a) = weights[cellist[a]];
			rho(ixshift, iyshift, izshift) = 1.;
		}
	}
}
int IsForeverAlone(int ix, int iy, int iz, Solid_list& solid_list, momentum_direction& e) {
	//checks if solid (surface) node is isolated. 
	int ixshift;
	int iyshift;
	int izshift;
	int friends = 26;
	for (int a = 0; a < 27; a++) {		// look through all nearby nodes
		ixshift = ix + e(a, 0);
		iyshift = iy + e(a, 1);
		izshift = iz + e(a, 2);
		if (solid_list(ixshift, iyshift, izshift) == -1)
			friends--;
	}
	if (friends == 0)
		return 1;
	else
		return 0;
}

/*
void updateBC_v2(direction_density& f, int t, double Bvel[3], density& rho, momentum_direction& e, velocity& u) {
int ix = 0;
int iy = 0;
int iz = 0;

for (iz = -1; iz < Nz + 1; iz++) {
for (iy = -1; iy < Ny + 1; iy++) {
for (ix = -1; ix < Nx + 1; ix++) {
if (isnodeBC(ix, iy, iz) == 1) {

}
}
}
}
}
int isnodeBC(int ix, int iy, int iz) {
if (ix == -1 || ix == Nx || iy == -1 || iy == Ny || iz == -1 || iz == Nz)
return 1;
else
return 0;
}
*/
