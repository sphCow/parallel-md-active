/*
 * BoundaryDataTransfer.cpp
 *
 *  Created on: 19-Mar-2015
 *      Author: theo_xubuntu
 */

#include "BoundaryDataTransfer.h"

//#include ""

vector<int> BoundaryDataTransfer::GetNeghibourRanks(int ix, int iy,
		int n_cell_x, int n_cell_y) {
	//Neghibour(int ix, int iy, int n_cell_x, int n_cell_y) {
	vector<int> neghibour_ranks;
	neghibour_ranks.resize(8);
	neghibour_ranks[0] = ((ix + 1) % n_cell_x
			+ ((iy + n_cell_y) % n_cell_y) * n_cell_x); //B
	neghibour_ranks[1] = ((ix + 1) % n_cell_x
			+ ((iy + n_cell_y + 1) % n_cell_y) * n_cell_x); //C
	neghibour_ranks[2] = ((ix + n_cell_x) % n_cell_x
			+ ((iy + n_cell_y + 1) % n_cell_y) * n_cell_x); //F
	neghibour_ranks[3] = ((ix + n_cell_x - 1) % n_cell_x
			+ ((iy + n_cell_y + 1) % n_cell_y) * n_cell_x);
	neghibour_ranks[4] = ((ix + n_cell_x - 1) % n_cell_x
			+ ((iy + n_cell_y) % n_cell_y) * n_cell_x); //H
	neghibour_ranks[5] = ((ix + n_cell_x - 1) % n_cell_x
			+ ((iy + n_cell_y - 1) % n_cell_y) * n_cell_x); //G
	neghibour_ranks[6] = ((ix + n_cell_x) % n_cell_x
			+ ((iy + n_cell_y - 1) % n_cell_y) * n_cell_x); //D
	neghibour_ranks[7] = ((ix + 1) % n_cell_x
			+ ((iy + n_cell_y - 1) % n_cell_y) * n_cell_x); //A

	return neghibour_ranks;
}

//constructor
BoundaryDataTransfer::BoundaryDataTransfer(Block& b_, const int nbx,
		const int nby, const MPI_Datatype& ParticleType_) {
	b = b_;
	ParticleType = ParticleType_;
	neghibour_ranks = GetNeghibourRanks(b_.rank % nbx, b_.rank / nbx, nbx, nby);
}

void BoundaryDataTransfer::populate_particles(Cell &i, vector<Particle>& v) {
	//v.clear();
	for (auto j = i.particles.begin(); j != i.particles.end(); j++) {
		v.push_back(*j);
	}
}

void BoundaryDataTransfer::FillSendBuffer(Block& b) {
	boundary_send_buffer.clear();
	boundary_send_buffer.resize(8);

	for (auto i = b.cells.begin(); i != b.cells.end(); i++) {
		uint cx = i->cell_id_axis[0];
		uint cy = i->cell_id_axis[1];

		//send north row
		if (cx > 0 && cx <= b.n_cell[0] && cy == b.n_cell[1]) {
			populate_particles(*i, boundary_send_buffer[2]);
		}

		//send south row
		if (cx > 0 && cx <= b.n_cell[0] && cy == 1) {
			populate_particles(*i, boundary_send_buffer[6]);
		}

		//send east column
		if (cy > 0 && cy <= b.n_cell[1] && cx == b.n_cell[0]) {
			populate_particles(*i, boundary_send_buffer[0]);
		}

		//send west column
		if (cy > 0 && cy <= b.n_cell[1] && cx == 1) {
			populate_particles(*i, boundary_send_buffer[4]); //0

		}

		//send north_east cell
		if (cx == b.n_cell[0] && cy == b.n_cell[1]) {
			//cout << myrank << " ne " << i->cell_id << endl;
			populate_particles(*i, boundary_send_buffer[1]);
		}

		//send north_west cell
		if (cx == 1 && cy == b.n_cell[1]) {
			populate_particles(*i, boundary_send_buffer[3]);
		}

		//send south_east cell
		if (cx == b.n_cell[0] && cy == 1) {
			populate_particles(*i, boundary_send_buffer[7]);
		}

		//send south_west cell
		if (cx == 1 && cy == 1) {
			populate_particles(*i, boundary_send_buffer[5]);
		}
	}
}

void BoundaryDataTransfer::SendToNeghibours() {
	for (auto myn = neghibour_ranks.begin(); myn != neghibour_ranks.end();
			myn++) {
		int pos = myn - neghibour_ranks.begin();

		//MPI::COMM_WORLD.Send(&boundary_send_buffer[pos].front(),
		//		boundary_send_buffer[pos].size(), ParticleType, *myn, pos);
		send_request[pos] = MPI::COMM_WORLD.Isend(&boundary_send_buffer[pos].front(), boundary_send_buffer[pos].size(), ParticleType, *myn, pos);
	}
}

void BoundaryDataTransfer::RecvFromNeghibours() {
	boundary_recv_buffer.clear();
	boundary_recv_buffer.resize(8);
	vector<MPI::Status> status(8);
	vector<int> rec_count(8);

	for (auto myn = neghibour_ranks.begin(); myn != neghibour_ranks.end();
			myn++) {
		int pos = myn - neghibour_ranks.begin();
		int sender_rank;

		if (pos < 4)
			sender_rank = *(myn + 4); //0-3
		else
			sender_rank = *(myn - 4); //4-7

		//MPI::COMM_WORLD.Send(&send_buffer[pos].front(), send_buffer[pos].size(), ParticleType, *myn, pos);

		MPI::COMM_WORLD.Probe(sender_rank, pos, status[pos]);
		rec_count[pos] = status[pos].Get_count(ParticleType);
		boundary_recv_buffer[pos].resize(rec_count[pos]);
		MPI::COMM_WORLD.Recv(&boundary_recv_buffer[pos].front(), rec_count[pos],
				ParticleType, sender_rank, pos);
	}

	//iprobe
	//		for (auto myn = bc.neghibour_ranks.begin(); myn < bc.neghibour_ranks.begin() + 4; myn++) {
	//			int pos = myn - bc.neghibour_ranks.begin();
	//			bool flag = false;
	//
	//			while (!flag) {
	//				flag = MPI::COMM_WORLD.Iprobe(*(myn + 4), pos, status[pos]);
	//				rec_count[pos] = status[pos].Get_count(ParticleType);
	//				recv_buffer[pos].resize(rec_count[pos]);
	//			}
	//
	//			MPI::COMM_WORLD.Recv(&recv_buffer[pos].front(), rec_count[pos], ParticleType, *(myn + 4), pos);
	//		}

}

void BoundaryDataTransfer::CopyParticles(Block& b, vector<Particle> in,
		int dir) {
	for (auto i = in.begin(); i != in.end(); i++) {
		int old_cx = i->cell_id_axis[0];
		int old_cy = i->cell_id_axis[1];
		int new_cid, new_cx, new_cy;

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
		b.cells[new_cid].particles.push_back(*i);
	}
}

void BoundaryDataTransfer::FillGhostCells(Block& b) {
	for (int i = 0; i < 8; i++) {
		CopyParticles(b, boundary_recv_buffer[i], i);
		boundary_recv_buffer[i].clear();
	}

}

const vector<vector<Particle>> BoundaryDataTransfer::SendRecv(const vector<vector<Particle>>& send_buffer) {

	vector<vector<Particle>> recv_buffer(8);
	vector<MPI::Status> status(8);
	vector<int> rec_count(8);

	for (auto myn = neghibour_ranks.begin();
			myn < neghibour_ranks.begin() + 4; myn++) {
		int pos = myn - neghibour_ranks.begin();

		MPI::COMM_WORLD.Send(&send_buffer[pos].front(),
				send_buffer[pos].size(), ParticleType, *myn, pos);

		MPI::COMM_WORLD.Probe(*(myn + 4), pos, status[pos]);
		rec_count[pos] = status[pos].Get_count(ParticleType);
		recv_buffer[pos].resize(rec_count[pos]);
		MPI::COMM_WORLD.Recv(&recv_buffer[pos].front(), rec_count[pos],
				ParticleType, *(myn + 4), pos);

		//CopyParticles(b, recv_buffer[pos], pos);
	}

	//}

	for (auto myn = neghibour_ranks.begin() + 4;
			myn != neghibour_ranks.end(); myn++) {
		int pos = myn - neghibour_ranks.begin();

		MPI::COMM_WORLD.Send(&send_buffer[pos].front(),
				send_buffer[pos].size(), ParticleType, *myn, pos);

		MPI::COMM_WORLD.Probe(*(myn - 4), pos, status[pos]);
		rec_count[pos] = status[pos].Get_count(ParticleType);
		recv_buffer[pos].resize(rec_count[pos]);
		MPI::COMM_WORLD.Recv(&recv_buffer[pos].front(), rec_count[pos],
				ParticleType, *(myn - 4), pos);

		//CopyParticles(b, recv_buffer[pos], pos);
	}

	return recv_buffer;

}
