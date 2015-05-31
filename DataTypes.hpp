/*
 * DataTypes.hpp
 *
 *  Created on: 19-Mar-2015
 *      Author: theo_xubuntu
 */

#include <vector>
#include <math.h>
#include <random>
#include <mpi.h>
using namespace std;

#ifndef DATATYPES_HPP_
#define DATATYPES_HPP_

//std::mt19937 mt_gen;
//std::uniform_real_distribution<double> random_theta(-M_PI, M_PI);

struct Particle {
	unsigned long tstep; //1
	unsigned int part_id, proc_id, cell_id, cell_id_axis[2]; //5
	double r[2], v[2], theta, f[2]; //7 double
};

struct MinimalParticle {
	unsigned long tstep;
	double r[2], theta;
};

typedef struct {
	unsigned int cell_id;
	unsigned int cell_id_axis[2];
	vector<Particle> particles;
	vector<Particle> particles_new;
} Cell;

typedef struct {
	unsigned int rank;
	unsigned int coordinate[2];
	unsigned int n_cell[2]; //n_cell_x & n_cell_y
	vector<Cell> cells; // vector containg the cells in the block
	double xlim[2]; //xlim[0]->start , xlim[1]->end
	double ylim[2]; //ylim[0]->start , ylim[1]->end
	unsigned int total_particle; // # particle in the Block
} Block;

typedef struct {
	unsigned int nx, ny, n_blocks_x, n_blocks_y, is_confined; //5
	unsigned long tmax, wait_for_equilibration, config_save_interval; //3
	double rho, sigma, rc, Pe, dt; //5
} Header;


#define IAmRoot myrank==0


#endif /* DATATYPES_HPP_ */
