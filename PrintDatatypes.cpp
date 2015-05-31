#include "DataTypes.hpp"

#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

string print_particle(Particle & p) {
	std::stringstream ss;
	ss << fixed;
//	ss << setw(8) << "#pid" << " " << setw(12) << "cid" << " " << setw(8) << "r[0]" << " " << setw(8) << "r[1]" << " "
//			<< setw(8) << "theta" << " " << setw(8) << "f[0]" << " " << setw(8) << "f[1]" << " " << '\n';

//for (int i = 0; i < p.size(); i++) {
	ss << setprecision(3) << setw(4) << p.part_id << " " << setprecision(3) << setw(8) << p.proc_id << " "
			<< setprecision(3) << setw(8) << p.cell_id << " {" << setprecision(3)
			//<< setw(3) << p.cell_id_axis[0] << " " << setprecision(3) << setw(3) << p.cell_id_axis[1] << " }" << setprecision(3) << setw(8) << p.r[0] << " " << setprecision(3)
			<< setw(8) << p.r[1] << " " << setprecision(3) << setw(8) << p.theta << " " << setprecision(15)
			<< setw(20)
			<< p.f[0] << " " << setprecision(15) << setw(20) << p.f[1] << " " << '\n';
//}
//ss << setw(30) << "#-----x-----" << '\n';

	return ss.str();
}

string print_particle_vector(vector<Particle> &p) {
	std::stringstream ss;
	ss << fixed;
//ss << setw(8) << "#pid" << " " << setw(8) << "cid" << " " << setw(8) << "r[0]" << " " << setw(8) << "r[1]" << " " << setw(8) << "theta" << " " << setw(8) << "f[0]" << " "
//		<< setw(8) << "f[1]" << " " << '\n';

	for (uint i = 0; i < p.size(); i++) {
//		ss << setprecision(3) << setw(8) << p[i].part_id << " " << setprecision(3) << setw(8) << p[i].cell_id << " {" << p[i].cell_id_axis[0] << "," << p[i].cell_id_axis[1] << "} "
//				<< setprecision(3) << setw(8) << p[i].r[0] << " " << setprecision(3) << setw(8) << p[i].r[1] << " " << setprecision(3) << setw(8) << p[i].theta << " "
//				<< setprecision(3) << setw(8) << p[i].f[0] << " " << setprecision(3) << setw(8) << p[i].f[1] << " " << '\n';
		ss << print_particle(p[i]);
	}
	ss << setw(30) << "#-----x-----" << '\n';

	return ss.str();
}

string print_cells_vector(vector<Cell> &cv) {
	stringstream ss;
	for (unsigned i = 0; i < cv.size(); i++) {
		//	ss << "# part in cell -> " << cv[i].particles.size() << " cell_id -> " << cv[i].cell_id << " {" << cv[i].cell_id_axis[0] << ", " << cv[i].cell_id_axis[1] << "} " << endl;
		ss << print_particle_vector(cv[i].particles);
	}

	return ss.str();

}

string print_block(Block b) {
	stringstream ss;

	ss << setw(30) << " --- block_begin ---" << endl;
	//ss << "rank -> " << myrank << " " << b.rank << endl;
	ss << "cart_coord -> " << b.coordinate[0] << ", " << b.coordinate[1] << endl;
	ss << "n_cells -> " << b.n_cell[0] << ", " << b.n_cell[1] << endl;
	ss << "x_lim -> " << b.xlim[0] << " to " << b.xlim[1] << endl;
	ss << "y_lim -> " << b.ylim[0] << " to " << b.ylim[1] << endl;
//cout << "offsetx -> " << b.offset[0] << " offset -> " << b.offset[1] << endl;
	ss << "total particles -> " << b.total_particle << endl;
	ss << "particles -> " << endl;
	ss << print_cells_vector(b.cells);

//for(auto b.cells.)
//print_particle_vector(c.particles);
	ss << setw(30) << " --- block_end ---" << endl;

	return ss.str();

}
