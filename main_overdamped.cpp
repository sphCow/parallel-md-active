/*
 ============================================================================
 Name        : mpi_md_mkII.cpp
 Author      : parswa
 Version     : -1 !!
 Copyright   : Copyleft
 Description : Parallel Molecular Dynamics
 ============================================================================
 */

/* http://stackoverflow.com/questions/10419990/creating-an-mpi-datatype-for-a-structure-containing-pointers
 *
 *
 *
 */

#include "DataTypes.hpp"
#include "BoundaryDataTransfer.h"

#include <math.h>
#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <ctime>
#include <vector>
#include <algorithm>
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <sys/types.h>

#include <random>
#include <chrono>

using namespace std;

extern string print_block(Block b);
extern string print_particle_vector(vector<Particle>& p);
extern string print_cells_vector(vector<Cell>& cv);
extern string print_particle(Particle &p);

void CopyParticles(Block& b, vector<Particle> in, int dir);

vector<int> get_neghibour_cell(int ix, int iy, int n_cell_x, int n_cell_y);
vector<int> get_neghibour_cell2(Cell& ic, Block& b);

void get_force_LJ(Particle& i, Particle& j, string s);
void get_force_LJ2(Particle& i, Particle& j);

inline bool i_am_not_ghost(Cell ic, Block b) {
	return (ic.cell_id_axis[0] > 0 && ic.cell_id_axis[0] <= b.n_cell[0] && ic.cell_id_axis[1] > 0
			&& ic.cell_id_axis[1] <= b.n_cell[1]);
}

bool sort_comp(Particle i, Particle j) {
	return (i.part_id < j.part_id);
}

long seedgen(void);
map<string, string> read_input();

extern MPI_Datatype make_ParticleType();
extern MPI_Datatype make_HeaderType();
extern MPI_Datatype make_MinimalParticleType();

//------------------- GENERAL VARIABLES ------------------------//
unsigned int nx, ny, n_blocks_x, n_blocks_y;
unsigned long n, max_time_steps;
double rho, dt, sigma, rc, Pe;
bool is_confined;

string input_data_path, output_data_path;
int debug_type = 0;

double ax, ay, lx, ly, lxh, lyh, rc2, D, D_r, tau, v_p, F_p;
double noise_var, nv, nvr, lx_cell, ly_cell;
unsigned int n_cell, n_cell_x, n_cell_y, n_part_per_proc;
unsigned int myrank, maxrank, mypid;

unsigned long config_save_interval, wait_for_equilibration;
int start_mode;
unsigned long t, config_saved = 0;

Header h_in;

//MPI_sizes
int particle_size, header_size, minimal_particle_size;

//debug

ofstream flog; //(log_fname);
//-------------------------------------------------------------//
/* --------------- SETUP PARTICLE DATATYPE FOR MPI -------------------	 */
/*	------------------------------------------------------------------	 */

int main(int argc, char *argv[]) {

	//std::mt19937 mt_gen;
	//std::random_device rd;
	//mt_gen.seed(28569);

	//----------------------------------------------------------------//
	MPI::Init(argc, argv);
	maxrank = MPI::COMM_WORLD.Get_size();
	myrank = MPI::COMM_WORLD.Get_rank();
	mypid = getpid();
	//cout << myrank << " " << mypid << endl;

	//-------------------- COMMIT MPI_DATATYPES ---------------------//
	MPI_Datatype ParticleType = make_ParticleType();
	MPI_Datatype HeaderType = make_HeaderType();
	MPI_Datatype MinimalParticleType = make_MinimalParticleType();
	MPI_Type_commit(&ParticleType);
	MPI_Type_commit(&HeaderType);
	MPI_Type_commit(&MinimalParticleType);

	MPI_Type_size(ParticleType, &particle_size);
	MPI_Type_size(HeaderType, &header_size);
	MPI_Type_size(MinimalParticleType, &minimal_particle_size);

	//--------------------------- file I/O ---------------------------------//
	string log_fname = "log_" + std::to_string(myrank) + ".txt";

	if (IAmRoot) flog.open(log_fname.c_str(), std::ofstream::out);

	string debug_fname = "debug_" + std::to_string(myrank) + ".txt";
	ofstream fdebug;

	//
	//MPI_File mfile;
//	MPI_File mfile;
//	MPI_File_open(MPI_COMM_WORLD, "out", MPI_MODE_CREATE | MPI_MODE_WRONLY,
//	MPI_INFO_NULL, &mfile);

	unsigned long offset = 0;

	//int psize, header_size;
	//MPI_Type_size(ParticleType, &psize);
	//MPI_Type_size(HeaderType, &header_size);

	//random
	std::mt19937 mt_gen;
	long myseed = seedgen();
	mt_gen.seed(myseed);
	std::uniform_real_distribution<double> random_theta(-M_PI, M_PI);
	std::normal_distribution<double> normal_rand(0.0, 1.0); // N(mean, stddeviation)

	vector<int> total_particle_vector;
	total_particle_vector.resize(maxrank);

	int start_mode; // 0= scratch, 1=load

	int in_fname_len_;
	int out_fname_len_;
	vector<char> infname_vchar;
	vector<char> outfname_vchar;
	string outfname_str;
	string infname_str;

	/* ------------------------------------------------------------------------
	 * ----------------READ/SET PARAMETERS FROM ROOT PROCESS-------------------
	 ---------------------------------------------------------------------- */
	if (myrank == 0) {
		map<string, string> in_map;
		in_map = read_input(); //from cin

		n_blocks_x = stoi(in_map["n_blocks_x"]);
		max_time_steps = stoi(in_map["max_time_steps"]);
		config_save_interval = stoi(in_map["config_save_interval"]);
		wait_for_equilibration = stoi(in_map["wait_for_equilibration"]);

		output_data_path = in_map["output_data_path"];
		std::copy(output_data_path.begin(), output_data_path.end(), back_inserter(outfname_vchar));
		out_fname_len_ = outfname_vchar.size();

		if (in_map["debug_type"] == "force") debug_type = 1;
		else if (in_map["debug_type"] == "none") debug_type = 0;

		if (in_map["start_mode"] == "new") {
			start_mode = 0;

			nx = stoi(in_map["nx"]);
			ny = stoi(in_map["ny"]);
//			n_blocks_x = stoi(in_map["n_blocks_x"]);
			rho = stod(in_map["rho"]);
			sigma = stod(in_map["sigma"]);
			rc = stod(in_map["rc"]);
			Pe = stod(in_map["Pe"]);
			is_confined = false; //in_map["is_confined"]==

			dt = stod(in_map["dt"]);

			cout << "Simulating from scratch..." << endl;

		} else if (in_map["start_mode"] == "load") {
			start_mode = 1;
			input_data_path = in_map["input_data_path"];
			cout << "Loading data from " << input_data_path << endl;

			std::copy(input_data_path.begin(), input_data_path.end(), back_inserter(infname_vchar));

			//ca = &v[0]; // pointer to start of vector
			in_fname_len_ = infname_vchar.size();

		} else {
			cerr << "start_mode error..." << endl;
			return -1;
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		//len_ = output_data_path.length();
		//ca = strcpy((char*)malloc(len_+1), output_data_path.c_str());

	}

	//-------------------- SET PARAMETRS ON ALL PROCESSES -----------------------//

	//broadcast mode
	MPI::COMM_WORLD.Bcast(&start_mode, 1, MPI::INT, 0);
	MPI::COMM_WORLD.Bcast(&debug_type, 1, MPI::INT, 0);

	//broadcast output filename
	MPI::COMM_WORLD.Bcast(&out_fname_len_, 1, MPI::INT, 0);
	if (myrank != 0) outfname_vchar.resize(out_fname_len_);

	MPI::COMM_WORLD.Bcast(&outfname_vchar.front(), outfname_vchar.size(), MPI::CHAR, 0);
	outfname_vchar.push_back('\0'); // null terminator

	outfname_str = string(outfname_vchar.begin(), outfname_vchar.end());
	char* char_out_fname;
	char* char_in_fname;
	char_out_fname = &outfname_vchar.front();
	//-				-- x --

	//open output file
	MPI_File mfile;
	MPI_File_open(MPI_COMM_WORLD,
			char_out_fname,
			MPI_MODE_CREATE | MPI_MODE_WRONLY,
			MPI_INFO_NULL,
			&mfile);

	//input file
	MPI_File minfile;

	//if loading from prev file
	if (start_mode == 1) {
		MPI::COMM_WORLD.Bcast(&in_fname_len_, 1, MPI::INT, 0);
		if (myrank != 0) infname_vchar.resize(in_fname_len_);

		MPI::COMM_WORLD.Bcast(&infname_vchar.front(), infname_vchar.size(), MPI::CHAR, 0);
		infname_vchar.push_back('\0'); // null terminator

		infname_str = string(infname_vchar.begin(), infname_vchar.end());
		char_in_fname = &infname_vchar.front();

		//open input mpi_file
		MPI_File_open(MPI_COMM_WORLD,
				char_in_fname,
				MPI_MODE_RDONLY,
				MPI_INFO_NULL, &minfile);

		//read the header
		MPI_File_read_at(minfile, 0, &h_in, 1, HeaderType, MPI_STATUS_IGNORE);

		//cout << h_in.dt << " " << h_in.ny << endl;
		nx = h_in.nx;
		ny = h_in.ny;
		rho = h_in.rho;
		dt = h_in.dt;
		sigma = h_in.sigma;
		rc = h_in.rc;
		Pe = h_in.Pe;

		//int in_wait_for_eq = h_in.wait_for_equilibration
		//cout << h_in.config_save_interval << " " << h_in.wait_for_equilibration << " "
		//		<< h_in.tmax<<endl;

	}

	//if from scratch
	else if (start_mode == 0) {
		MPI::COMM_WORLD.Bcast(&nx, 1, MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(&ny, 1, MPI::DOUBLE, 0);
//		MPI::COMM_WORLD.Bcast(&n_blocks_x, 1, MPI::INT, 0);
		MPI::COMM_WORLD.Bcast(&rho, 1, MPI::DOUBLE, 0);
//		MPI::COMM_WORLD.Bcast(&max_time_steps, 1, MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(&dt, 1, MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(&sigma, 1, MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(&rc, 1, MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(&Pe, 1, MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(&is_confined, 1, MPI::BOOL, 0);
//		MPI::COMM_WORLD.Bcast(&config_save_interval, 1, MPI::INT, 0);
//		MPI::COMM_WORLD.Bcast(&wait_for_equilibration, 1, MPI::INT, 0);
	}

	MPI::COMM_WORLD.Bcast(&n_blocks_x, 1, MPI::INT, 0);
	MPI::COMM_WORLD.Bcast(&max_time_steps, 1, MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast(&config_save_interval, 1, MPI::INT, 0);
	MPI::COMM_WORLD.Bcast(&wait_for_equilibration, 1, MPI::INT, 0);

	//cout << particle_size << " " << header_size << endl;
	//---------------------------- set derived variables ----------------------------//
	n_blocks_y = maxrank / n_blocks_x;
	n = nx * ny;
	ax = sqrt(2.0 / (rho * sqrt(3.0)));	//sqrt(2.0/(rho*sqrt(3.0)));
	ay = ax * sqrt(3.0) / 2.0;
	lx = nx * ax;
	ly = ny * ay;
	lxh = lx / 2.0;
	lyh = ly / 2.0;

	//cut-offs
	rc2 = rc * rc;

	//active
	D = 1.0;
	D_r = 3.0 * D / (sigma * sigma);
	tau = sigma * sigma / D;
	v_p = Pe * sigma / tau;
	F_p = v_p / (D);

	//noise
	noise_var = sqrt(2.0 * dt);	//np.sqrt(2.0 / dt)
	nv = sqrt(2.0 * D * dt);
	nvr = sqrt(2.0 * D_r * dt);

	//cells
	n_cell_x = floor(lx / rc);		//2;
	n_cell_y = floor(ly / rc);		//2;
	n_cell = n_cell_x * n_cell_y;
	lx_cell = lx / double(n_cell_x);
	ly_cell = ly / double(n_cell_y);
	int n_cells_per_block_x = ceil(n_cell_x / float(n_blocks_x));
	int n_cells_per_block_y = ceil(n_cell_y / float(n_blocks_y));

	t = 0;

	//header
	Header header;

	MPI::COMM_WORLD.Barrier();

	/* ---------------------------------------------------------------------------
	 * ------------------------ MAP CELLS TO PROCESSORS --------------------------
	 * --------------------------------------------------------------------------*/

	Block b; //block is a collection of cells;
	b.rank = myrank;

	//rank -> coordinate
	b.coordinate[0] = myrank % n_blocks_x;
	b.coordinate[1] = myrank / n_blocks_x;

	// # of cells in the block
	b.n_cell[0] =
			(b.coordinate[0] == (n_blocks_x - 1)) ?
													(n_cell_x - n_cells_per_block_x * (n_blocks_x - 1)) :
													(n_cells_per_block_x);
	b.n_cell[1] =
			(b.coordinate[1] == (n_blocks_y - 1)) ?
													(n_cell_y - n_cells_per_block_y * (n_blocks_y - 1)) :
													(n_cells_per_block_y);

	//cout << myrank << " {" << b.coordinate[0] << "," << b.coordinate[1] << "} " << b.n_cell[0] << " " << b.n_cell[1] << endl;

	int total_cell_in_block = (b.n_cell[0] + 2) * (b.n_cell[1] + 2);

	b.cells.resize(total_cell_in_block);

	//assign cell id
	for (unsigned i = 0; i < b.n_cell[0] + 2; i++) { // +2 for ghosts cells
		for (unsigned j = 0; j < b.n_cell[1] + 2; j++) {
			int cell_id = i + j * (b.n_cell[0] + 2); //cx + cy * b.n_cell[0];
			b.cells[cell_id].cell_id = cell_id;
			b.cells[cell_id].cell_id_axis[0] = i;
			b.cells[cell_id].cell_id_axis[1] = j;
		}
	}

	b.xlim[0] = -lxh + n_cells_per_block_x * b.coordinate[0] * lx_cell;
	b.ylim[0] = -lyh + n_cells_per_block_y * b.coordinate[1] * ly_cell;
	b.xlim[1] = b.xlim[0] + b.n_cell[0] * lx_cell;
	b.ylim[1] = b.ylim[0] + b.n_cell[1] * ly_cell;

	double xbin = n_cells_per_block_x * lx_cell; //need modification
	double ybin = n_cells_per_block_y * ly_cell;

	MPI::COMM_WORLD.Barrier();

	if (myrank == 0) {
		//print input to cout
		cout << setw(20) << setprecision(5) << "system size -> " << nx << " x " << ny << endl;
		cout << setw(20) << setprecision(5) << "decomposition -> " << n_blocks_x << " x " << n_blocks_y << endl;
		cout << setw(20) << setprecision(5) << "rho -> " << rho << endl;
		cout << setw(20) << setprecision(5) << "sigma -> " << sigma << endl;
		cout << setw(20) << setprecision(5) << "rc -> " << rc << endl;
		cout << setw(20) << setprecision(5) << "Pe -> " << Pe << endl;
		cout << setw(20) << setprecision(5) << "max_time_steps -> " << max_time_steps << endl;
		cout << setw(20) << setprecision(5) << "dt -> " << dt << endl;
		cout << setw(20) << setprecision(5) << "config_save_interval -> " << config_save_interval << endl;
		cout << setw(20) << setprecision(5) << "wait_for_equilibration -> " << wait_for_equilibration << endl;
		cout << endl;
		cout << setw(20) << setprecision(5) << "ax -> " << ax << endl;
		cout << setw(20) << setprecision(5) << "ay -> " << ay << endl;
		cout << setw(20) << setprecision(5) << "lx -> " << lx << endl;
		cout << setw(20) << setprecision(5) << "ly -> " << ly << endl;
		cout << setw(20) << setprecision(5) << "cells -> " << n_cell_x << " x " << n_cell_y << endl;
		cout << setw(20) << setprecision(5) << "lx_cell -> " << lx_cell << endl;
		cout << setw(20) << setprecision(5) << "ly_cell -> " << ly_cell << endl;
		cout << endl << endl;
		cout << setw(10) << setprecision(5) << "out_file -> " << output_data_path << endl;
		cout << setw(10) << setprecision(5) << "debug_type -> " << debug_type << endl;
		cout << endl << endl;

		//set header
		header.nx = nx;
		header.ny = ny;
		header.n_blocks_x = n_blocks_x;
		header.n_blocks_y = n_blocks_y;
		header.is_confined = is_confined;

		header.tmax = max_time_steps;
		header.wait_for_equilibration = wait_for_equilibration;
		header.config_save_interval = config_save_interval;

		header.rho = rho;
		header.sigma = sigma;
		header.rc = rc;
		header.Pe = Pe;
		header.dt = dt;

		//write header
		MPI_File_write_at(mfile, 0, &header, 1, HeaderType, MPI_STATUS_IGNORE);
	}

	//MPI_File_seek(mfile,header_size,MPI_SEEK_SET);
	MPI::COMM_WORLD.Barrier();

	/* ---------------------------------------------------------------------------------
	 * ------------------------ GENERATE/LOAD LATTICE  ------------------------
	 ----------------------------------------------------------------------------------- */

//	//generate the full lattice in all ranks
	vector<vector<Particle>> lattice;
	lattice.resize(maxrank);

	if (start_mode == 0) {
		//generate lattice here
		for (unsigned i = 0; i < nx; i++) {
			for (unsigned j = 0; j < ny; j++) {

				double x = -lxh + i * ax;
				double y = -lyh + j * ay;
				if (j % 2 != 0)
					x += ax / 2.0;

				Particle *p = new Particle;
				p->tstep = 0;
				p->part_id = i + j * nx;

				p->r[0] = x;
				p->r[1] = y;

				p->theta = random_theta(mt_gen); //0.0; //RANDOM

				int bx = int((x + lxh) / xbin);
				int by = int((y + lyh) / ybin);
				p->proc_id = by * n_blocks_x + bx;

				p->f[0] = 0.0;
				p->f[1] = 0.0;

				lattice.at(p->proc_id).push_back(*p);

				//_debug_
				//lat_.at(p->part_id) = (*p);

				delete p; //Annihilate the particle!!!
			}
		}

		//-------> load last config from input MPI_File <-------------
	} else if (start_mode == 1) {
		unsigned long nconfig = (h_in.tmax - h_in.wait_for_equilibration) / (h_in.config_save_interval) - 1;
		//cout << nconfig << endl;

		MPI_Offset last_config_offset = (nconfig - 1) * h_in.nx * h_in.ny * particle_size + header_size;
		vector<Particle> config(h_in.nx * h_in.ny);

		MPI_File_open(MPI_COMM_WORLD,
				char_in_fname,
				MPI_MODE_RDONLY,
				MPI_INFO_NULL, &minfile);

		MPI_File_read_at(minfile,
				last_config_offset,
				&config.front(),
				h_in.nx * h_in.ny,
				ParticleType,
				MPI_STATUS_IGNORE);

		sort(config.begin(), config.end(), sort_comp);

		for (auto &p : config) {
			p.tstep = 0;

			int bx = int((p.r[0] + lxh) / xbin);
			int by = int((p.r[1] + lyh) / ybin);
			p.proc_id = by * n_blocks_x + bx;

			p.f[0] = 0.0;
			p.f[1] = 0.0;

			lattice.at(p.proc_id).push_back(p);

		}

		//print loaded lattice
		if (IAmRoot) {
			ofstream floaded("loaded_lattice.txt", ios::out);
			floaded << "# total particles -> " << config.size() << "serial interaction f" << endl << endl;
			//blocks plot
			//A blank line in the data file -> break in the line connecting data points in gnuplot.
			for (unsigned i = 0; i < n_blocks_x; i++) {
				floaded << -lxh + i * n_cells_per_block_x * lx_cell << " " << -lyh << endl;
				floaded << -lxh + i * n_cells_per_block_x * lx_cell << " " << +lyh << endl << endl;
			}

			floaded << -lxh << " " << lyh << endl;
			floaded << lxh << " " << lyh << endl;
			floaded << lxh << " " << -lyh << endl << endl;

			for (unsigned i = 0; i < n_blocks_y; i++) {
				floaded << -lxh << " " << -lyh + i * n_cells_per_block_y * ly_cell << endl;
				floaded << lxh << " " << -lyh + i * n_cells_per_block_y * ly_cell << endl << endl;
			}

			//cells plot
			floaded << endl << endl;
			//A blank line in the data file -> break in the line connecting data points in gnuplot.
			for (unsigned i = 0; i <= n_cell_x; i++) {
				floaded << -lxh + i * lx_cell << " " << -lyh << endl;
				floaded << -lxh + i * lx_cell << " " << +lyh << endl << endl;
			}

			for (unsigned i = 0; i <= n_cell_y; i++) {
				floaded << -lxh << " " << -lyh + i * ly_cell << endl;
				floaded << lxh << " " << -lyh + i * ly_cell << endl << endl;
			}
			floaded << endl << endl;

			//---------- serial force debug---------------------//
			if (debug_type == 1) {
				cout << "debugging interaction force..." << endl;
				for (auto &i : config) {
					for (auto &j : config) {
						if (i.part_id != j.part_id) {
							get_force_LJ(i, j, "d");
						}
					}
				}
			}
			//---------------------------------------------------//

			//print to file
			for (auto &p : config) {
				floaded << p.part_id << "\t" << p.r[0] << "\t" << p.r[1] << "\t"
						<< p.f[0] << "\t" << p.f[1] << endl;
			}

			floaded.close();
		}

	}

	MPI::COMM_WORLD.Barrier();
//----------------------------------------------------------------------------//

	//put parts in cell
	for (auto i = lattice.at(myrank).begin(); i != lattice.at(myrank).end(); i++) {
		int cx = (i->r[0] - b.xlim[0]) / lx_cell + 1;
		int cy = (i->r[1] - b.ylim[0]) / ly_cell + 1;
		int cell_id = (cx) + (cy) * (b.n_cell[0] + 2); //n_cell_x;
		i->cell_id = cell_id; //change
		i->cell_id_axis[0] = cx;
		i->cell_id_axis[1] = cy;

		b.cells[cell_id].particles.push_back(*i);
		//b.cells[cell_id].cell_id=cell_id;
	}

	//get total # of particles
	b.total_particle = 0;
	for (auto i = b.cells.begin(); i != b.cells.end(); i++) {
		b.total_particle += i->particles.size();
	}

	MPI::COMM_WORLD.Barrier();

	/* ---------------------------------------------------------------------------------
	 * -------------------------- SEND/RECEIVE "GHOST" CELLS !! ------------------------
	 ----------------------------------------------------------------------------------- */
	//Neghibour nb(myrank % n_blocks_x, myrank / n_blocks_x, n_blocks_x, n_blocks_y);
	BoundaryDataTransfer bcdt(b, n_blocks_x, n_blocks_y, ParticleType); //(b, myrank); //get blocks where this block b should send data to

	//BoundaryCellsSend bc;

	vector<MPI::Status> status;
	vector<MPI_Request> req;
	req.resize(8);

	MPI_Request *send_req;

	vector<int> rec_count;

	vector<vector<Particle>> send_buffer1;
	vector<vector<Particle>> send_buffer;

	while (t < max_time_steps) {

		vector<vector<Particle>> recv_buffer;

		//set forces to 0 for all particles
		for (auto &ic : b.cells) {
			for (auto &ip : ic.particles) {
				ip.f[0] = 0;
				ip.f[1] = 0;
				ip.tstep = t;
			}
		}

		//send
		bcdt.FillSendBuffer(b);
		bcdt.SendToNeghibours();
		//recv_buffer = bcdt.SendRecv(bcdt.boundary_send_buffer);

		//calculate interaction forces in own cell
		for (auto & ic : b.cells) { //loop over all cells ic in block b
			if (i_am_not_ghost(ic, b)) { // exclude ghost cells
				for (auto i = ic.particles.begin(); i != ic.particles.end(); i++) {
					for (auto j = i + 1; j != ic.particles.end(); j++) {
						get_force_LJ2(*i, *j);
					}
				}
			}
		}

		send_req = bcdt.getSendRequest();
		MPI_Waitall(8, send_req, MPI_STATUSES_IGNORE);

		//RECEIVE
		bcdt.RecvFromNeghibours();
		recv_buffer = bcdt.getBoundaryRecvBuffer();

		for (int i = 0; i < 8; i++) {
			CopyParticles(b, recv_buffer[i], i);
			recv_buffer[i].clear();
		}

		MPI::COMM_WORLD.Barrier();

		/* ---------------------------------------------------------------------------------
		 * -------------------------- CALCULATE INTERACTION FORCES -------------------------
		 ----------------------------------------------------------------------------------- */

		//calculate interaction forces with neghibour cells
		for (auto & ic : b.cells) { //loop over all cells ic in block b

			vector<int> ic_neghibour = get_neghibour_cell2(ic, b);

			for (auto & i : ic.particles) {
				for (auto & kc : ic_neghibour) {
					for (auto & j : b.cells.at(kc).particles) {
						get_force_LJ2(i, j);
					}
				}
			}
		}

		//MPI::COMM_WORLD.Barrier();

		//--------- clear ghost cells ----------------
		for (auto & ic : b.cells) { //loop over all cells ic in block b
			if (!i_am_not_ghost(ic, b)) {
				ic.particles.clear(); // remove particles in ghost cells
			}
		}

		MPI::COMM_WORLD.Barrier();

		//DEBUG :: STORE INTERACTION FORCES TO TEXT FILE & EXIT
		if (debug_type == 1) {
			vector<int> n_part_all_blocks;
			n_part_all_blocks.resize(maxrank);

			MPI::COMM_WORLD.Allgather(&b.total_particle, 1, MPI::INT, &n_part_all_blocks.front(), 1, MPI::INT);

			int total_part_in_system = 0;
			for (unsigned i = 0; i < maxrank; i++)
				total_part_in_system += n_part_all_blocks[i];

			vector<Particle> all_parts_system(total_part_in_system);
			int displ[maxrank];

			if (IAmRoot) {
				int sum = 0;
				for (unsigned i = 0; i < maxrank; ++i) {
					displ[i] = sum;
					sum += n_part_all_blocks[i];
				}
			}

			//store all particles in a vector in each block
			vector<Particle> all_parts_block;
			for (auto &ic : b.cells) {
				for (auto &ip : ic.particles) {
					all_parts_block.push_back(ip);
				}
			}

			MPI::COMM_WORLD.Gatherv(&all_parts_block.front(),
					n_part_all_blocks[myrank],
					ParticleType,
					&all_parts_system.front(),
					&n_part_all_blocks.front(),
					&displ[0],
					ParticleType, 0);

			if (IAmRoot) {
				ofstream fparallel("force_parallel.txt", ios::out);
				sort(all_parts_system.begin(), all_parts_system.end(), sort_comp);
				for (auto &p : all_parts_system) {
					fparallel << p.part_id << "\t" << p.r[0] << "\t" << p.r[1] << "\t"
							<< p.f[0] << "\t" << p.f[1] << endl;
				}

				fparallel.close();
				cout << "interaction force debug has completed..." << endl;
			}

			MPI::COMM_WORLD.Barrier();

			return 0;
			//
			MPI_Finalize();
			MPI_Abort(MPI_COMM_WORLD, 0);
		}

		//----- END DEBUG -----------//

		/* ---------------------------------------------------------------------------------
		 * ------------------------ INTEGRATE & UPDATE CELL LISTS --------------------------
		 ----------------------------------------------------------------------------------- */

		send_buffer.clear();
		send_buffer.resize(8);

		for (auto & ic : b.cells) {
			//ic.particles_new.clear();
			//if (i_am_not_ghost(ic, b)) {
			for (auto ip = ic.particles.begin(); ip != ic.particles.end(); ip++) {

				unsigned int old_bid = ip->proc_id;

				ip->r[0] += D * (ip->f[0] + F_p * cos(ip->theta)) * dt + nv * normal_rand(mt_gen);
				ip->r[1] += D * (ip->f[1] + F_p * sin(ip->theta)) * dt + nv * normal_rand(mt_gen);
				ip->theta += nvr * normal_rand(mt_gen); //sqrt(2 * D_r) * noise_var * normal_rand(mt_gen);

				//pbc
				ip->r[0] -= lx * rint(ip->r[0] / lx);
				ip->r[1] -= ly * rint(ip->r[1] / ly);

				//check if the particle is in block b
				unsigned int new_bx = int((ip->r[0] + lxh) / xbin);
				unsigned int new_by = int((ip->r[1] + lyh) / ybin);
				unsigned int new_bid = new_by * n_blocks_x + new_bx;

				//particle moved to another block
				if (old_bid != new_bid) {
					ip->proc_id = new_bid;
					int dir;

					//find its neghibour id: e=0; ne=1; n=2; nw=3; w=4; sw=5; s=6; se=7
					vector<int> nb = bcdt.neghibour_ranks;
					//vector<int> nb = get_neghibour_cell(old_bx, old_by, n_blocks_x, n_blocks_y);

					auto pos = find(nb.begin(), nb.end(), new_bid); // - bc.begin();

					if (pos == nb.end()) {
						cerr << t << " i've moved too far: " << old_bid << " -> " << new_bid << endl;
						flog << t << " i've moved too far: " << old_bid << " -> " << new_bid << endl;

						MPI::COMM_WORLD.Abort(0);
						MPI::Finalize();
						return -1;
					} else {
						dir = std::distance(nb.begin(), pos);
						send_buffer.at(dir).push_back(*ip);
					}
				}

				else {
					unsigned int old_cid = ic.cell_id;  //ip.cell_id;
					unsigned int new_cx = (ip->r[0] - b.xlim[0]) / lx_cell + 1;
					unsigned int new_cy = (ip->r[1] - b.ylim[0]) / ly_cell + 1;
					unsigned int new_cid = (new_cx) + (new_cy) * (b.n_cell[0] + 2); //n_cell_x;

					if (old_cid != new_cid) {
						//add p in new cell
						ip->cell_id = new_cid;
						ip->cell_id_axis[0] = new_cx;
						ip->cell_id_axis[1] = new_cy;
						b.cells[new_cid].particles_new.push_back(*ip);
						//flog<<"upd"<< print_particle(*ip) << old_cid<<endl;
					} else {
						ic.particles_new.push_back(*ip);
					}

				}
			}
			//}
		}

		//MPI::COMM_WORLD.Barrier();

		//SEND
		bcdt.setBlock(b);
		bcdt.setBoundarySendBuffer(send_buffer);
		bcdt.SendToNeghibours();

		//do something here

		send_req = bcdt.getSendRequest();
		MPI_Waitall(8, send_req, MPI_STATUSES_IGNORE);

		//RECEIVE
		bcdt.RecvFromNeghibours();
		recv_buffer = bcdt.getBoundaryRecvBuffer();

		//fallback
//		recv_buffer = bcdt.SendRecv(send_buffer);

		MPI::COMM_WORLD.Barrier();

		//put received particles in respective cells + update ip.cell_id
		for (auto &j : recv_buffer) {
			for (auto &p : j) {
				unsigned int new_cx = (p.r[0] - b.xlim[0]) / lx_cell + 1;
				unsigned int new_cy = (p.r[1] - b.ylim[0]) / ly_cell + 1;
				unsigned int new_cid = (new_cx) + (new_cy) * (b.n_cell[0] + 2); //n_cell_x;

				p.cell_id = new_cid;
				p.cell_id_axis[0] = new_cx;
				p.cell_id_axis[1] = new_cy;
				b.cells[new_cid].particles_new.push_back(p);
			}
		}

		//swap
		for (auto &ic : b.cells) {
			ic.particles.swap(ic.particles_new);
			ic.particles_new.clear();
		}

		//update total_number
		b.total_particle = 0;
		for (auto i = b.cells.begin(); i != b.cells.end(); i++) {
			b.total_particle += i->particles.size();
		}

//		int debug_t_part;
//		MPI::COMM_WORLD.Reduce(&b.total_particle,&debug_t_part,1,MPI::INT,MPI::SUM,0);
//		if(IAmRoot) cout << debug_t_part << endl;

		if (myrank == 0) {
			//cout.flush();
			//if (fmod(percent_done, 5.0) == 0) flog << percent_done << endl;
			if (fmod(t, max_time_steps * 0.1 / 100) == 0) {
				float percent_done = t * 100.0 / float(max_time_steps);

				time_t now = time(0);

				tm *ltm = localtime(&now);

				cout << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec
						<< " :->> "
						<< setprecision(4) << setw(4)
						<< percent_done << " :->> "
						<< setw(8) << t << "/"
						<< max_time_steps
						<< "\n";
			}

		}

		//debug
		MPI::COMM_WORLD.Barrier();

		/* ---------------------------------------------------------------------------------
		 * ------------------------------------ FILE I/O ------------------------------------
		 ----------------------------------------------------------------------------------- */

		//MPI_Offset ofs;
		//MPI_File_get_position(mfile,&ofs);
		//cout << myrank << " " << "ofs -> " << ofs << endl;
		if (t > wait_for_equilibration && fmod(t, config_save_interval) == 0) {

			MPI::COMM_WORLD.Allgather(&b.total_particle,
					1, MPI_INT,
					&total_particle_vector.front(),
					1, MPI_INT);

			MPI_Status stat;

			offset = 0;			//0;

			for (unsigned i = 0; i < myrank; i++) {
				offset += total_particle_vector[i] * particle_size;
			}

			MPI::COMM_WORLD.Barrier();

			vector<Particle> write_buffer;
			for (auto ib = b.cells.begin(); ib != b.cells.end(); ib++) {
				if (i_am_not_ghost(*ib, b)) {
					for (auto ip = ib->particles.begin(); ip != ib->particles.end(); ip++) {
						write_buffer.push_back(*ip);
					}
				}
			}

			MPI_Offset rof = n * particle_size * config_saved + offset + header_size;
			MPI_File_write_at_all(mfile,
					rof, //t * n * particle_size / config_save_interval + offset,
					&write_buffer.front(),
					write_buffer.size(),
					ParticleType,
					&stat);

			MPI::COMM_WORLD.Barrier();
			config_saved++;

		}

		//quits the job if there's file named "quit_md" in the bin dir
		//if (IAmRoot) {
		if (t % 5000 == 0) {
			bool b = std::ifstream("quit_md").good();
			if (b) {
				cout << "qutting...\n";
				MPI_File_close(&mfile);
				MPI::Finalize();
				return -1;
			}
		}
		//}

		t++;
		//flog << endl << endl;

	}

	MPI_File_close(&mfile);

	if (myrank == 0) {
		cout << "Total " << config_saved << " configs were saved... " << endl;
		flog << "Total " << config_saved << " configs were saved... " << endl;
		cout << "Your wait is over... job is done." << endl;
		flog << "Your wait is over... job is done." << endl;
	}

	MPI::Finalize();
	return 0;
}

/* --------------------------------------------------------------------------------------
 * -------------------------------- FUNCTION DEFINITIONS --------------------------------
 ---------------------------------------------------------------------------------------- */

void CopyParticles(Block& b, vector<Particle> in, int dir) {
	/* DIRECTION KEYS
	 *
	 * SENDER	RECEIVER	dir
	 * east		west		0
	 * west		east		1
	 *
	 *
	 */

	for (auto i = in.begin(); i != in.end(); i++) {
		int old_cx = i->cell_id_axis[0];
		int old_cy = i->cell_id_axis[1];
		int new_cid, new_cx, new_cy; // = get_new_cid(cx,cy,dir); //(new_cx) + (new_cy) * (b.n_cell[0] + 2);

		//e->w
		if (dir == 0) {
			new_cx = 0; //
			new_cy = old_cy; //unchanged
		}

		// w -> e
		if (dir == 4) {
			new_cx = b.n_cell[0] + 1; //
			new_cy = old_cy; //unchanged
		}

		// n -> s
		if (dir == 2) {
			new_cx = old_cx; //
			new_cy = 0; //unchanged
		}

		// s -> n
		if (dir == 6) {
			new_cx = old_cx; //
			new_cy = b.n_cell[1] + 1; //unchanged
		}

		// ne -> sw
		if (dir == 1) {
			new_cx = 0; //
			new_cy = 0; //old_cy - b.n_cell[1]; //unchanged
		}

		//nw->se
		if (dir == 3) {
			new_cx = b.n_cell[0] + 1; //
			new_cy = 0; //old_cy - b.n_cell[1]; //unchanged
		}

		//sw->ne
		if (dir == 5) {
			new_cx = b.n_cell[0] + 1; //
			new_cy = b.n_cell[1] + 1; //old_cy - b.n_cell[1]; //unchanged
		}

		//se->nw
		if (dir == 7) {
			new_cx = 0; //
			new_cy = b.n_cell[1] + 1; //old_cy - b.n_cell[1]; //unchanged
		}

		new_cid = (new_cx) + (new_cy) * (b.n_cell[0] + 2);

		i->cell_id = new_cid;
		//i->proc_id = 255;
		//cout << i->part_id << " " << old_cid << " " << new_cid << endl;
		b.cells[new_cid].particles.push_back(*i);
		//cout << i_am_not_ghost(b.cells[new_cid],b) << endl;
	}
}

vector<int> get_neghibour_cell(int ix, int iy, int n_cell_x, int n_cell_y) {

	vector<int> neghibour_ranks;		//(4);
	neghibour_ranks.resize(8);

//e
	neghibour_ranks[0] = ((ix + 1) % n_cell_x + ((iy + n_cell_y) % n_cell_y) * n_cell_x);		//e
//ne
	neghibour_ranks[1] = ((ix + 1) % n_cell_x + ((iy + n_cell_y + 1) % n_cell_y) * n_cell_x);			//C
//n
	neghibour_ranks[2] = ((ix + n_cell_x) % n_cell_x + ((iy + n_cell_y + 1) % n_cell_y) * n_cell_x);			//F
//nw
	neghibour_ranks[3] = ((ix + n_cell_x - 1) % n_cell_x + ((iy + n_cell_y + 1) % n_cell_y) * n_cell_x);
//w
	neghibour_ranks[4] = ((ix + n_cell_x - 1) % n_cell_x + ((iy + n_cell_y) % n_cell_y) * n_cell_x);			//H
//sw
	neghibour_ranks[5] = ((ix + n_cell_x - 1) % n_cell_x + ((iy + n_cell_y - 1) % n_cell_y) * n_cell_x);			//G
//s
	neghibour_ranks[6] = ((ix + n_cell_x) % n_cell_x + ((iy + n_cell_y - 1) % n_cell_y) * n_cell_x);			//D
//se
	neghibour_ranks[7] = ((ix + 1) % n_cell_x + ((iy + n_cell_y - 1) % n_cell_y) * n_cell_x);			//A
//own cell
//	neghibour_ranks[8] = (ix + iy * n_cell_x); //base_cell

	return neghibour_ranks;
}

vector<int> get_neghibour_cell2(Cell& ic, Block& b) {

	vector<int> nranks;
	int ix = ic.cell_id_axis[0];
	int iy = ic.cell_id_axis[1];
	int ncx = b.n_cell[0] + 2;
	int ncy = b.n_cell[1] + 2;

//bool id_in_lim = [](int id_) {return (id_< lim);};

//e
	int e_ = (ix + 1) + iy * ncx;
//ne
	int ne_ = (ix + 1) + (iy + 1) * ncx;
//n
	int n_ = ix + (iy + 1) * ncx;
//nw
	int nw_ = (ix - 1) + (iy + 1) * ncx;

	if ((ix + 1) < ncx)
		nranks.push_back(e_);

	if ((ix + 1) < ncx && (iy + 1) < ncy)
		nranks.push_back(ne_);

	if ((iy + 1) < ncy)
		nranks.push_back(n_);

	if ((ix - 1) >= 0 && (iy + 1) < ncy)
		nranks.push_back(nw_);

	return nranks;

}

void get_force_LJ(Particle& i, Particle& j, string s) {
	double xr, yr, rr, r2i, r6i, ff;

	xr = i.r[0] - j.r[0];
	xr -= lx * rint(xr / lx);

	yr = i.r[1] - j.r[1];
	yr -= ly * rint(yr / ly);

	rr = xr * xr + yr * yr;

	if (rr < rc2) {
		r2i = 1.0 / rr;
		r6i = r2i * r2i * r2i;
		//ff = 12 * r6i * r6i * r6i;
		//ff = 12.0*r2i*r6i*r6i; //soft
		ff = 48.0 * r2i * r6i * (r6i - 0.5);			//#LJ

		if (ff > 10000) {
			cout << "# warn from " << myrank << " " << s << endl;
			cout << "# " << i.part_id << " " << i.proc_id << " " << i.cell_id << " :-> " << j.part_id << " "
					<< j.proc_id << " " << j.cell_id << endl;
		}

		i.f[0] += ff * xr;
		i.f[1] += ff * yr;
	}

//CONFINEMENT
//	if (is_confined) {
//		double dy, y5, y10;
//
//		for (int i = 0; i < n; i++) {
//			dy = iy - ylo;
//			if (fabs(dy) < yc) {
//				y5 = dy * dy * dy * dy * dy;
//				y10 = y5 * y5;
//				fy += 10.0 / y10 / dy - 4.0 / y5;
//			}
//
//			dy = iy - yup;
//			if (fabs(dy) < yc) {
//				y5 = dy * dy * dy * dy * dy;
//				y10 = y5 * y5;
//				fy += 10.0 / y10 / dy - 4.0 / y5;
//			}
//		}
//	}
}

void get_force_LJ2(Particle& i, Particle& j) {
	double xr, yr, rr, r2i, r6i, ff;

	xr = i.r[0] - j.r[0];
	xr -= lx * rint(xr / lx);

	yr = i.r[1] - j.r[1];
	yr -= ly * rint(yr / ly);

	rr = xr * xr + yr * yr;

	if (rr < rc2) {
		r2i = 1.0 / rr;
		r6i = r2i * r2i * r2i;
		//ff = 12 * r6i * r6i * r6i;
		//ff = 12.0*r2i*r6i*r6i; //soft
		ff = 48.0 * r2i * r6i * (r6i - 0.5);			//#LJ

		if (ff > 10000) {
			cout << "# warn from " << myrank << endl;
			cout << "# " << i.part_id << " " << i.proc_id << " " << i.cell_id << " :-> " << j.part_id << " "
					<< j.proc_id << " " << j.cell_id << endl;
		}

		i.f[0] += ff * xr;
		i.f[1] += ff * yr;
		j.f[0] -= ff * xr;
		j.f[1] -= ff * yr;
	}
}

//Ref http://arxiv.org/pdf/1005.4117v1.pdf
long seedgen(void) {
	long s, seed, pid, seconds;
	pid = getpid(); /* get process ID */

	s = time(&seconds); /* get CPU seconds since 01/01/1970 */
	seed = fabs(((s * 181) * ((pid - 83) * 359)) % 104729);

	return seed;

}

map<string, string> read_input()
{
	map<string, string> in_map;
//string input_data_path, output_data_path;
//double b;

	in_map =
	{	{	"n_blocks_x",""},
		{	"nx",""},
		{	"ny",""},
		{	"rho",""},
		{	"sigma",""},
		{	"rc",""},
		{	"Pe",""},
		{	"is_confined",""},
		{	"max_time_steps",""},
		{	"dt",""},
		{	"config_save_interval",""},
		{	"wait_for_equilibration",""},
		{	"start_mode",""},
		{	"input_data_path",""},
		{	"output_data_path",""},
		{	"debug_type",""}};

	std::string line;
	string input_lines;
	string field;

	while (std::getline(std::cin, line)) {
		if (line == "end")
			break;

		for (auto &i : in_map) {
			int begin_pos = line.find(i.first);

			if (begin_pos != -1) {
				begin_pos += i.first.size() + 3;
				field = line.substr(begin_pos);

				if (i.first == "input_data_path") {
					int begin_quote = field.find("\"");
					int end_quote = field.find("\"", begin_quote + 1);
					field = field.substr(begin_quote + 1, end_quote - begin_quote - 1);
				}

				if (i.first == "output_data_path") {
					int begin_quote = field.find("\"");
					int end_quote = field.find("\"", begin_quote + 1);
					field = field.substr(begin_quote + 1, end_quote - begin_quote - 1);

				}

				if (i.first == "start_mode" ||
						i.first == "debug_type") {
					int begin_quote = field.find("\"");
					int end_quote = field.find("\"", begin_quote + 1);
					field = field.substr(begin_quote + 1, end_quote - begin_quote - 1);
				}

				i.second = (field);
			}
		}

	}

	return in_map;
}
