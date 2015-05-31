/*
 * BoundaryDataTransfer.h
 *
 *  Created on: 19-Mar-2015
 *      Author: theo_xubuntu
 */


#include <vector>
#include <mpi.h>
#include <math.h>
#include <random>

#include "DataTypes.hpp"

using namespace std;

#ifndef BOUNDARYDATATRANSFER_H_
#define BOUNDARYDATATRANSFER_H_

class BoundaryDataTransfer {
public:
	BoundaryDataTransfer(Block& b_, const int nbx, const int nby,
			const MPI_Datatype& ParticleType_);

	void FillSendBuffer(Block& b);

	void SendToNeghibours();
	void RecvFromNeghibours();

	const vector<vector<Particle>> SendRecv(const vector<vector<Particle>>& send_buffer);


	vector<int> neghibour_ranks;

	void FillGhostCells(Block &b);

	//getter & setter for send buffer
	const vector<vector<Particle> >& getBoundarySendBuffer() const {
		return boundary_send_buffer;
	}

	void setBoundarySendBuffer(const vector<vector<Particle> >& boundarySendBuffer) {
		boundary_send_buffer = boundarySendBuffer;
	}

	//setter for Block
	void setBlock(Block &b_) {
		b = b_;
	}

	const Block& getBlock() const {
		return b;
	}

	//getter for send_request
	MPI_Request* getSendRequest() {
		return send_request;
	}

	const vector<vector<Particle>>& getBoundaryRecvBuffer() const {
		return boundary_recv_buffer;
	}

	vector<vector<Particle>> boundary_send_buffer;

private:
	void populate_particles(Cell &i, vector<Particle>& v);
	vector<int> GetNeghibourRanks(int ix, int iy, int n_cell_x, int n_cell_y);
	void CopyParticles(Block& b, vector<Particle> in, int dir);

	Block b;
	int myrank;
	MPI_Datatype ParticleType;
	MPI_Request send_request[8];

	vector<vector<Particle>> boundary_recv_buffer;


};

#endif /* BOUNDARYDATATRANSFER_H_ */
