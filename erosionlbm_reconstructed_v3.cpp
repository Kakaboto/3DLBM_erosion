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




// DATA LAYOUT FOR f AND ftemp:
// f: 
// (i,j,k):		(1,1,1):(Nx,Ny,Nz) - physical points. (0,0,0), (Nx + 1, Ny + 1, Nz + 1) - Boundary points ("ghost points").

// solid_list:	1 = fluid point. 0 = surface solid point. -1 = interior solid point.

//Frågor:	Varför kan jag inte printa ut värden från solid_list å f när dom är const?

// Antar att Nx = Ny = Nz (isotropt grid).



using namespace std;
typedef vector<double> dvec;


int Nx;
int Ny;
int Nz;


double c = 1.;
double tau = 6.;
double mu = 1. / 3.*(tau - 0.5);
double umax_theo = 0.1;
double Re = 0.01;
//double L = 1. - 1. / ((double)Nx - 1.);	//Assuming same L in x and y direction.
double L = 0.5;
										//double gg = 48 * umax_theo / (L*L);
double gg = Re*mu / (tau*L);

double sq3 = 1.732050807568877;


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
		latspace_x = 1. / (Nx - 1.);
		latspace_y = 1. / (Ny - 1.);
		latspace_z = 1. / (Nz - 1.);
		for (ix = 0; ix < Nx; ix++)
			x[ix] = xmin + latspace_x*ix;
		for (iy = 0; iy < Ny; iy++)
			y[iy] = ymin + latspace_y*iy;
		for (iz = 0; iz < Nz; iz++)
			z[iz] = zmin + latspace_z*iz;
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
	Solid_list(int choice, Grid& grid, momentum_direction& e);
	~Solid_list();
	int& operator()(int ix, int iy, int iz) {
		return element[ix + iy*Nx + iz*Nx*Ny];
	}

	void printsolid_list(FILE * solfile) {
		int ix = 0;
		int iy = 0;
		int iz = 0;
		int n = 0;
		for (iz = 0; iz < Nz; iz++) {
			cout << "\n";
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					fprintf(solfile, "%i\n", element[n]);
//					cout << element[n] << " ";
					n++;
				}
			}
		}
		cout << "\n Printed grid to file for figures. \n";
	}

private:
	double latspace = 1. / (Nx - 1); //assuming Nx = Ny = Nz
	double xmin = 0; double xmax = 1;
	double ymin = 0; double ymax = 1;
	double zmin = 0; double zmax = 1;
	int Nxtot = Nx + 2;
	int Nytot = Ny + 2;
	int Nztot = Nz + 2;
	int Ncubed = Nx*Ny*Nz;
	int Ncubedtot = Nxtot*Nytot*Nztot;

};

Solid_list::Solid_list(int choice, Grid& grid, momentum_direction& e) {
	element.resize(Ncubed);
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;
	int n = 0;
	int nnext;
	int xmid = 0;
	int ymid = 0;
	int surfacecheck = 0;
	string tempstr;
	double radius;
	double radiussq;
	double radiussqh;
	double distsq;
	dvec center;
	center.resize(3);

	switch (choice) {
	case 1:	//sphere
		cout << "\n Sphere chosen! \n";
		cout << "\n Input radius: ";
		cin >> radius;
		//		stringstream(tempstr) >> radius;
		if (radius*2 > 1) {
			cout << "Radius too large! must be between 0 and 1";
			break;
		}
		cout << "\n Input center (x,y,z goes from (0,0,0) to (1,1,1)): ";
		cin >> center[0];
		cin >> center[1];
		cin >> center[2];
		/*		getline(cin, tempstr);
		stringstream(tempstr) >> center[0];
		getline(cin, tempstr);
		stringstream(tempstr) >> center[1];
		getline(cin, tempstr);
		stringstream(tempstr) >> center[2];*/
		if ((center[0] > 1) || (center[1] > 1) || (center[2] > 1)) {
			cout << "Center outside of grid! must be between 0 and 1";
			break;
		}


		radiussq = radius*radius;
		radiussqh = (radius - sq3*latspace)*(radius - sq3*latspace);
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					distsq = (grid.x[ix] - center[0])*(grid.x[ix] - center[0]) + (grid.y[iy] - center[1])*(grid.y[iy] - center[1]) + (grid.z[iz] - center[2])*(grid.z[iz] - center[2]);
					if (distsq > radiussq)
						element[n] = -1;
					else
						element[n] = 1;
					/*					if (distsq <= radiussq) {
					if (distsq > radiussqh)
					element[n] = 0;
					else
					element[n] = 1;
					}*/
					n++;
				}
			}
		}
		cout << "\n Created obstacle: Sphere! \n";
		n = 0;
		break;

	case 2:
		cout << "Cylinder chosen! \n Input radius: ";
		getline(cin, tempstr);
		stringstream(tempstr) >> radius;
		if (radius * 2 > 1) {
			cout << "Radius too large! must be between 0 and 1";
			break;
		}
		radiussq = radius*radius;
		radiussqh = (radius + sq3*latspace)*(radius + sq3*latspace);
		xmid = floor(((double)Nx - 1.)*0.5);
		ymid = floor(((double)Ny - 1.)*0.5);
		n = 0;
		for (iz = 0; iz < Nz; iz++) {	// can optimize by doing this for 1 iz, then just copiyng it to all other places in the vector.
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					distsq = (grid.x[ix] - grid.x[xmid])*(grid.x[ix] - grid.x[xmid]) + (grid.y[iy] - grid.y[ymid])*(grid.y[iy] - grid.y[ymid]);
					if (distsq < radiussq)
						element[n] = -1;
					else
						element[n] = 1;
					n++;
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
					if ((ix == 0) || (ix == Nx - 1) || (iy == 0) || (iy == Ny - 1))
						element[n] = 1;
					else
						element[n] = -1;
					n++;
				}
			}
		}
		n = 0;
		break;
	default:
		cout << "\n Error. incorrect input for obstacle type. \n";
	}
	// Placing out surface points.
	for (iz = 0; iz < Nz; iz++) {	// can optimize by doing this for 1 iz, then just copiyng it to all other places in the vector.
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				n = ix + iy*Nx + iz*Nx*Ny;
				if (element[n] == 1) {
					for (a = 0; a < 27; a++) {
						if (((ix + e(a, 0)) >= 0) && ((ix + e(a, 0)) < Nx) && ((iy + e(a, 1)) >= 0) && ((iy + e(a, 1)) < Ny) && ((iz + e(a, 2)) >= 0) && ((iz + e(a, 2)) < Nz)) {
							nnext = (ix + e(a, 0)) + (iy + e(a, 1))*Nx + (iz + e(a, 2))*Nx*Ny;
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
	cout << "Object deleted. \n";
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
	double * element;		//size: 3*Ncubed


	velocity() {
		//element.resize(3 * Ncubed);
		element = new double[3 * Ncubed];
	};

	int index(int ix, int iy, int iz, int i) {
		int ind = ix + iy*Nx + iz*Nx*Ny + i;
		return ind;
	}

	double& operator()(int ix, int iy, int iz, int i) {
		return element[3*(ix + iy*Nx + iz*Nx*Ny) + i];
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
	double * element;		//size: Ncubed

	density() {
		element = new double[Ncubed];
	};

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



void IC(const momentum_direction& e, direction_density& f, direction_density& ftemp, Solid_list& solid_list);
void macrovariables(velocity& u, density& rho, Solid_list& solid_list, direction_density& f, momentum_direction& e);
void edf(Solid_list& solid_list, velocity& u, density& rho, EDF& feq, momentum_direction& e);
double edfvecdot(momentum_direction& e, double ueq[3], int a);
void stream(Solid_list& solid_list, direction_density& f, direction_density& ftemp, momentum_direction& e);
void collision(Solid_list& solid_list, direction_density& f, direction_density& ftemp, EDF& feq);
void updatePBC(direction_density& f);
void printstuff(FILE * velfile, FILE * densfile, FILE * parfile, int Renum, velocity& u, density& rho);


int main()
{
	FILE * velfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\velocity_reconstructed.txt", "w");
	FILE * densfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\density_reconstructed.txt", "w");
	FILE * gridfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\grid_reconstructed.txt", "w");
	FILE * parfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\parameters_reconstructed.txt", "w");
	FILE * solfile = fopen("C:\\Users\\lukas\\Documents\\MATLAB\\Erosion_LBM\\solid_list_reconstructed.txt", "w");
	int obchoice;
	int tend;



	cout << "Input Nx, Ny, Nz: \n";
	//cin >> Nx;
	//cin >> Ny;
	//cin >> Nz;
	Nx = 30;
	Ny = 30;
	Nz = 30;
	Grid grid;
	momentum_direction e;
	EDF feq;
	velocity u;
	density rho;
	direction_density f;
	direction_density ftemp;						// Kör med 30x30x30 grid, cylinder med 0.45 i radius och en sfär i mitten av cylindern.
	grid.printgrid(gridfile);

	cout << "Choose object. 1: sphere. 2: cylinder. 3: square pipe. 4: cylinder with sphere inside.\n";
	//cin >> obchoice;
	obchoice = 2;
	Solid_list solid_list(obchoice, grid, e);
	solid_list.printsolid_list(solfile);

	cout << "Grid and object created. Running simulation for t = ";
	//cin >> tend;
	tend = 50;
	IC(e, f, ftemp, solid_list);
	for (int t = 0; t < tend; t++) {
		stream(solid_list, f, ftemp, e);
		macrovariables(u, rho, solid_list, f, e);
		printstuff(velfile, densfile, parfile, t, u, rho);
		edf(solid_list, u, rho, feq, e);
		collision(solid_list, f, ftemp, feq);
		cout << t << " ";
	}
	cout << "\n Done! \n";




	//	getchar();
	fclose(velfile);
	fclose(densfile);
	fclose(gridfile);
	fclose(parfile);
	fclose(solfile);

	return 0;
}





void IC(const momentum_direction& e, direction_density& f, direction_density& ftemp, Solid_list& solid_list) {
	// variables
	// -----------------------------
	int a = 0;
	int ix = 0;		// x 
	int iy = 0;		// y
	int iz = 0;		// z
					//------------------------------
					// function body
					//		cout << "\n IC: ------------------------------------" << "\n";
	for (a = 0; a < 27; a++) {
		//		cout << "\n a = " << a << "\n";
		for (iz = 0; iz < Nz; iz++) {
			//				cout << "\n";
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					if (solid_list(ix, iy, iz) == -1) {
						if (a == 13) {
							f(ix, iy, iz, a) = 1.;
							ftemp(ix, iy, iz, a) = 1.;
						}
						else {
							f(ix, iy, iz, a) = 0.;
							ftemp(ix, iy, iz, a) = 0.;
						}
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

	updatePBC(f);
	updatePBC(ftemp);

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
	int ix = 0;		// only used for is_node_solid, x
	int iy = 0;		// only used for is_node_solid, y
	int iz = 0;		// only used for is_node_solid, z
	int a = 0;		// counter for the 27 directions from one reference point.
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

					if (solid_list(ix, iy, iz) == -1) {	//LOOK AT SOLID_LIST ELEMENTS.
						rho(ix, iy, iz) += f(ix, iy, iz, a);
						u(ix, iy, iz, 0) += e(a, 0) * f(ix, iy, iz, a);
						u(ix, iy, iz, 1) += e(a, 1) * f(ix, iy, iz, a);
						u(ix, iy, iz, 2) += e(a, 2) * f(ix, iy, iz, a);
					}
//						cout << "x: " << u(ix,iy,iz,0) << ", y: " << u(ix,iy,iz,1) << ", z: " << u(ix,iy,iz,2) << "\n";
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
					u(ix, iy, iz, 0) /= rho(ix, iy, iz);
					u(ix, iy, iz, 1) /= rho(ix, iy, iz);
					u(ix, iy, iz, 2) /= rho(ix, iy, iz);
				}
			}
		}
	}
};


void edf(Solid_list& solid_list, velocity& u, density& rho, EDF& feq, momentum_direction& e) {
	// variables
	double f1 = 3.;
	double f2 = 9.*0.5;
	double f3 = 3.*0.5;
	double weights[4];
	int cellist[27] = { 3, 2, 3, 2, 1, 2, 3, 2, 3  ,  2, 1, 2, 1, 0, 1, 2, 1, 2  ,  3, 2, 3, 2, 1, 2, 3, 2, 3 };
	double w = 0.;
	double usq = 0;
	double ueq[3];
	double csq = 0;
	double dotprod = 0;
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;	// Direction counter. Used for feq and w. 
	weights[0] = 8. / 27.;	// Change name to weights
	weights[1] = 2. / 27.;
	weights[2] = 1. / 54.;
	weights[3] = 1. / 216.;
	// function body
//	cout << "\n EDF: -------------------------------------- \n";
	for (a = 0; a < 27; a++) {
//		cout << "\n a = " << a << ": \n \n";
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {

					if (solid_list(ix, iy, iz) == -1) {
						ueq[0] = u(ix, iy, iz, 0);
						ueq[1] = u(ix, iy, iz, 1);
						ueq[2] = u(ix, iy, iz, 2) + (tau*gg) / rho(ix, iy, iz);	//Add subroutine for this
						w = weights[cellist[a]];
						//							cout << "\n " << ux << " " << uy << " " << uz << "\n";
						usq = ueq[0] * ueq[0] + ueq[1] * ueq[1] + ueq[2] * ueq[2];
						csq = c*c;
						dotprod = edfvecdot(e, ueq, a);
						feq(ix, iy, iz, a) = rho(ix, iy, iz) * w * (1. + f1*dotprod / (csq)+f2*dotprod*dotprod / (csq*csq) - f3*usq / csq);
					}
//					cout << " " << feq(ix, iy, iz, a) << " ";
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


void stream(Solid_list& solid_list, direction_density& f, direction_density& ftemp, momentum_direction& e) {
	// variable defs
	int ix = 0;			// x
	int iy = 0;			// y
	int iz = 0;			// z
	int a = 0;			// direction counter
						// function body
//	cout << "\n Stream: ------------------------------------ \n";
	for (a = 0; a < 27; a++) {
//		cout << "\n a = " << a << "\n";
		for (iz = 0; iz < Nz; iz++) {
//			cout << "\n";
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {
					if (solid_list(ix, iy, iz) != 1) {
						// cout << " " << e[26 - a][0] << " " << e[26 - a][1] << " " << e[26 - a][2] << "\n";
						ftemp(ix, iy, iz, a) = f(ix + e(26 - a, 0), iy + e(26 - a, 1), iz + e(26 - a, 2), a);	//the e:s shifts ix,iy,iz to the node streaming into the node we're looking at. 
					}
//						cout << " " << ftemp(ix, iy, iz, a) << " ";
				}
			}
		}
	}
	updatePBC(ftemp);
	/*	for (a = 0; a < 27; a++) {
	cout << "\n a = " << a << "\n";
	for (iz = -1; iz < Nz + 1; iz++) {
	cout << "\n";
	for (iy = -1; iy < Ny + 1; iy++) {
	for (ix = -1; ix < Nx+1; ix++) {
	cout << " " << ftemp(ix, iy, iz, a) << " ";
	}
	}
	}
	}*/
}


void collision(Solid_list& solid_list, direction_density& f, direction_density& ftemp, EDF& feq) {				// Kolla så att collision function fungerar!
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int a = 0;
	int ashift = 0;
//		cout << "\n Collision: ------------------------------------\n \n";
	for (a = 0; a < 27; a++) {
		ashift = 2 * (13 - a);
//			cout << "\n a = " << a << "\n";
		for (iz = 0; iz < Nz; iz++) {
//				cout << "\n";
			for (iy = 0; iy < Ny; iy++) {
				for (ix = 0; ix < Nx; ix++) {					// Include a switch here.
					switch (solid_list(ix, iy, iz)) {
					case 1:	// Interior node
						break;
					case 0:	// Surface node
						f(ix, iy, iz, a + ashift) = ftemp(ix, iy, iz, a);
						break;
					case -1://	Fluid node
						f(ix, iy, iz, a) = ftemp(ix, iy, iz, a) - (ftemp(ix, iy, iz, a) - feq(ix, iy, iz, a)) / tau;
						break;
					}
//					cout << " " << f(ix, iy, iz, a) << " ";


					/*					if (solid_list(ix, iy, iz) != 1) {			// if it is not an interior solid node.
					if (solid_list(ix, iy, iz) == -1) {		// if it is a fluid node.
					// if fluid point, compute collisions.
					f(ix, iy, iz, a) = ftemp(ix, iy, iz, a) - (ftemp(ix, iy, iz, a) - feq(ix, iy, iz, a)) / tau;
					}
					else {
					//	If not fluid point, then surface node and bounceback BC.
					f(ix, iy, iz, a + ashift) = ftemp(ix, iy, iz, a);
					}
					}
					//						cout << " " << f[m] << " ";
					*/
				}
			}
		}
	}
	updatePBC(f);
};


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
			f.element[f.index(Nx, iy, iz + Nz, a)] = f(-1, iy, iz, a);

		// copying to forward upper edge
		for (ix = 0; ix < Nx; ix++)
			f.element[f.index(ix, -1, iz + Nz, a)] = f(ix, Ny - 1, iz, a);

		// copying to back upper edge
		for (ix = 0; ix < Nx; ix++)
			f.element[f.index(ix, Ny, iz + Nz, a)] = f(ix, -1, iz, a);

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
		f.element[f.index(-1, Ny, -1, a)] = f(Nx - 1, 0, Nz - 1, a);		// bottom left-back
																			//----------------------------------------------------------------------
																			// copy left side
		ix = Nx;
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				f.element[f.index(ix, iy, iz, a)] = f(0, iy, iz, a);
			}
		}
		// copy right side
		ix = -1;
		for (iz = 0; iz < Nz; iz++) {
			for (iy = 0; iy < Ny; iy++) {
				f.element[f.index(ix, iy, iz, a)] = f(ix + Nx, iy, iz, a);
			}
		}
		// copy forward side
		iy = Ny;
		for (iz = 0; iz < Nz; iz++) {
			for (ix = 0; ix < Nx; ix++) {
				f.element[f.index(ix, iy, iz, a)] = f(ix, 0, iz, a);
			}
		}
		// copy back side
		iy = -1;
		for (iz = 0; iz < Nz; iz++) {
			for (ix = 0; ix < Nx; ix++) {
				f.element[f.index(ix, iy, iz, a)] = f(ix, iy + Ny, iz, a);
			}
		}
	}
	//		break;
}


void printstuff(FILE * velfile, FILE * densfile, FILE * parfile, int Renum, velocity& u, density& rho) {
	int ix = 0;
	int iy = 0;
	int iz = 0;
	if (Renum == 0) {
		cout << "\n Re = " << Re << "\n";
		cout << "\n tau = " << tau << "\n";
		cout << "\n L = " << L << "\n";
		cout << "\n umax_theo = " << umax_theo << "\n";
		cout << "\n gg = " << gg << "\n";
		fprintf(parfile, "%e %e %e %e %e", Re, tau, L, umax_theo, gg);
	}
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				fprintf(velfile, "%e %e %e ", u(ix, iy, iz, 0), u(ix, iy, iz, 1), u(ix, iy, iz, 2));
				fprintf(densfile, "%e ", rho(ix, iy, iz));
			}
		}
	}
	fprintf(velfile, "\n");
	fprintf(densfile, "\n");
}


