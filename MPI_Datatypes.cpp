#include "DataTypes.hpp"

MPI_Aint ex_ulong, ex_int, ex_double;

MPI_Datatype make_ParticleType() {

//	unsigned long tstep; //1
//	unsigned int part_id, proc_id, cell_id, cell_id_axis[2]; //5
//	double r[2], v[2], theta, f[2]; //7 double

	struct Particle p[5];

	MPI_Datatype ParticleType;

	MPI_Datatype type[9] = {
	MPI_UNSIGNED_LONG,
	MPI_INT, MPI_INT, MPI_INT, MPI_INT,
	MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };

	int blocklen[9] = {
			1,
			1, 1, 1, 2,
			2, 2, 1, 2 };

	MPI_Aint disp[9];

	MPI_Type_extent(MPI_UNSIGNED_LONG, &ex_ulong);
	MPI_Type_extent(MPI_INT, &ex_int);
	MPI_Type_extent(MPI_DOUBLE, &ex_double);

	disp[0] = 0;
	disp[1] = ex_ulong;
	disp[2] = disp[1] + ex_int;
	disp[3] = disp[2] + ex_int;
	disp[4] = disp[3] + ex_int;
	disp[5] = disp[4] + 2 * ex_int;
	disp[6] = disp[5] + 2 * ex_double;
	disp[7] = disp[6] + 2 * ex_double;
	disp[8] = disp[7] + ex_double;

	MPI_Type_create_struct(9, blocklen, disp, type, &ParticleType);
	return ParticleType;
}

MPI_Datatype make_HeaderType() {
//	unsigned int nx, ny, n_blocks_x, n_blocks_y, is_confined; //5
//	unsigned long tmax, wait_for_equilibration, config_save_interval; //3
//	double rho, sigma, rc, Pe, dt; //5

	MPI_Datatype HeaderType;

	MPI_Datatype type[2] = {
	MPI_UNSIGNED_LONG, MPI_DOUBLE };

	int blocklen[2] = { 8, 5 };

	MPI_Aint disp[2];
	disp[0] = 0;

	MPI_Type_extent(MPI_UNSIGNED_LONG, &ex_ulong);
	MPI_Type_extent(MPI_DOUBLE, &ex_double);

//	for(int i=1; i<13; i++) {
//		MPI_Type_extent(type[i-1], &ext);
//		disp[i] = disp[i-1] + ext;
//	}

	disp[0] = 0;
	disp[1] = 8 * ex_ulong;

	MPI_Type_create_struct(2, blocklen, disp, type, &HeaderType);
	return HeaderType;
}

MPI_Datatype make_MinimalParticleType() {
	MPI_Datatype MinimalParticleType;
	MPI_Datatype type[3] = { MPI_UNSIGNED_LONG, MPI_DOUBLE, MPI_DOUBLE };
	int blocklen[3] = { 1, 2, 1 };
	MPI_Aint disp[3];

	MPI_Type_extent(MPI_UNSIGNED_LONG, &ex_ulong);
	MPI_Type_extent(MPI_DOUBLE, &ex_double);

	disp[0] = 0;
	disp[1] = ex_ulong;
	disp[2] = ex_ulong + 2 * ex_double;

	MPI_Type_create_struct(3,blocklen,disp,type,&MinimalParticleType);

	return MinimalParticleType;
}
