#ifndef MPI_WRAPPER_2D
#define MPI_WRAPPER_2D

#define _SCL_SECURE_NO_WARNINGS

#include "mpi.h"
#include <algorithm>
#include <iostream>

namespace MPI_Wrapper{

	// MPI Initialization wrapper
	int Parallel_Init(int* argc = nullptr, char*** argv = nullptr);

	// MPI Finalize wrapper
	int Parallel_Finalize();

	/*
	get the corresponding MPI_DataType.
	*/
	template<class Type>
	MPI_Datatype get_MPI_Datatype(){
		std::cerr << "get_MPI_Datatype: no specilization is used." << std::endl;
		return MPI_BYTE;
	}

	/*
	Declarations of template specification.
	*/
	template<>
	MPI_Datatype get_MPI_Datatype<short>();
	template<>
	MPI_Datatype get_MPI_Datatype<unsigned short>();
	template<>
	MPI_Datatype get_MPI_Datatype<int>();
	template<>
	MPI_Datatype get_MPI_Datatype<unsigned int>();
	template<>
	MPI_Datatype get_MPI_Datatype<long>();
	template<>
	MPI_Datatype get_MPI_Datatype<long long>();
	template<>
	MPI_Datatype get_MPI_Datatype<float>();
	template<>
	MPI_Datatype get_MPI_Datatype<double>();
	template<>
	MPI_Datatype get_MPI_Datatype<char>();

	/*
	class Parallel2D:
	a class wrapper of MPI, arrange the processes in a 2D grid.
	The indexing convention is [x,y], with y changes first.
	*/
	class Parallel2D {
	public:
		typedef unsigned int ArraySize;
		//CONSTRUCTOR AND DESTRUCTOR ========================
		/*
		The constructor will take care of all the initialization that is not related to 2D grid structure.
		*/
		Parallel2D(MPI_Comm world_comm = MPI_COMM_WORLD);
		Parallel2D(const int n_rows,
					const int n_cols,
					MPI_Comm world_comm = MPI_COMM_WORLD);
		~Parallel2D();

		//INITIALIZATION ====================================
		/*
		Initialize all the 2D grid structure parameters.
		\ row_size : size of the first dimension (dim0) of the MPI process grid.
		\ col_size : size of the second dimension (dim1) of the MPI process grid.
		*/
		void InitializeGrid(const int row_size, const int col_size);

		//ABORT =============================================
		int ForceAbort();
		//void AbortRequest();

		//BARRIER ===========================================
		int Barrier();

		//GLOBAL AND DIRECTIONAL PROCESSES COMMUNICATIONS ===

		//BROADCAST =========================================
		/*
		broadcast a message from process with world_rank_ == from to all the other processes.
		*/
		template<class Type> int Broadcast(Type& message, const int from) const;

		/*
		broadcast an array from process with world_rank_ == from to all the other processes.
		*/
		template<class Type> int Broadcast(Type* array, const int len, const int from) const;

		/*
		broadcast message/array in the same row.
		Processes with grid_rank_[0] == from will broadcast the variable to processes 
		with same grid_rank_[1] (i.e. the same row).
		*/
		template<class Type> int Broadcast_row(Type& message, const int from) const;
		template<class Type> int Broadcast_row(Type* array, const int len, const int from) const;

		/*
		broadcast message/array in same column.
		Processes with grid_rank_[1] == from will broadcast the variable to processes
		with same grid_rank_[0] (i.e. the same column).
		*/
		template<class Type> int Broadcast_col(Type& message, const int from) const;
		template<class Type> int Broadcast_col(Type* array, const int len, const int from) const;

		//GATHER ============================================
		/*
		Gather data to world_rank == gather_to.
		len: the size of array for each sender.
		*/
		template<class Type> 
		int Gather(const Type* send_array, const int len,
					Type* recv_array, const int gather_to) const;

		//SCATTER ===========================================
		/*
		Scatter data from world_rank == scatter_from.
		len: the size of array for each receiver.
		*/
		template<class Type>
		int Scatter(const Type* send_array, const int len,
					Type* recv_array, const int scatter_from) const;
		
		//REDUCE ============================================
		/*
		Perform reduce operation in all the processes.
		The result is substitued in-place.
		*/
		template<class Type> int Reduce(Type& message, const MPI_Op& op) const;
		template<class Type> int Reduce(Type* array, const int len, const MPI_Op& op) const;
		/*
		perform reduce operation in the same row.
		The reduce operation will be performed within the same grid_rank_[1] (i.e. same row).
		*/
		template<class Type> int Reduce_row(Type& message, const MPI_Op& op) const;
		template<class Type> int Reduce_row(Type* array, const int len, const MPI_Op& op) const;
		/*
		perform reduce operation in the same column.
		The reduce operation will be performed within the same grid_rank_[0] (i.e. same column).
		*/
		template<class Type> int Reduce_col(Type& message, const MPI_Op& op) const;
		template<class Type> int Reduce_col(Type* array, const int len, const MPI_Op& op) const;

		/*
		The SUM, MAX, MIN are derived from REDUCE.
		*/
		//SUM ===============================================
		template<class Type> int Sum(Type& message) const;
		template<class Type> int Sum(Type* array, const int len) const;
		template<class Type> int Sum_row(Type& message) const;
		template<class Type> int Sum_row(Type* array, const int len) const;
		template<class Type> int Sum_col(Type& message) const;
		template<class Type> int Sum_col(Type* array, const int len) const;

		//MAX ===============================================
		template<class Type> int Max(Type& message) const;
		template<class Type> int Max(Type* array, const int len) const;
		template<class Type> int Max_row(Type& message) const;
		template<class Type> int Max_row(Type* array, const int len) const;
		template<class Type> int Max_col(Type& message) const;
		template<class Type> int Max_col(Type* array, const int len) const;

		//MIN ===============================================
		template<class Type> int Min(Type& message) const;
		template<class Type> int Min(Type* array, const int len) const;
		template<class Type> int Min_row(Type& message) const;
		template<class Type> int Min_row(Type* array, const int len) const;
		template<class Type> int Min_col(Type& message) const;
		template<class Type> int Min_col(Type* array, const int len) const;

		//SEND ==============================================
		/*
		send a message/an array to the process (receiver) with world_rank_ == to.
		*/
		template<class Type> int Send(Type& message, const int to) const;
		template<class Type> int Send(Type* array, const int len, const int to) const;
		/*
		send an array to the process with grid_rank_ = to_gridrank.
		*/
		template<class Type> int Send(Type* array, const int len, const int* to_gridrank) const;
		/*
		send a message/an array to the processes (receiver) with grid_rank_[0] = to (in the same row).
		*/
		template<class Type> int Send_row(Type& message, const int to) const;
		template<class Type> int Send_row(Type* array, const int len, const int to) const;
		/*
		send a message/an array to the processes (receiver) with grid_rank_[1] = to (in the same column).
		*/
		template<class Type> int Send_col(Type& message, const int to) const;
		template<class Type> int Send_col(Type* array, const int len, const int to) const;
		//RECEIVE ===========================================
		/*
		receive a message/an array from the process (sender) with world_rank_ == from.
		*/
		template<class Type> int Receive(Type& message, const int from) const;
		template<class Type> int Receive(Type* array, const int len, const int from) const;
		/*
		receive an array from the process with grid_rank_ = from_gridrank
		*/
		template<class Type> int Receive(Type* array, const int len, const int* from_gridrank) const;
		/*
		receive a message/an array from the processes (sender) with grid_rank_[0] = from 
		(in the same row).
		*/
		template<class Type> int Receive_row(Type& message, const int from) const;
		template<class Type> int Receive_row(Type* array, const int len, const int from) const;
		/*
		receive a message/an array from the processes (sender) with grid_rank_[1] = from 
		(in the same column).
		*/
		template<class Type> int Receive_col(Type& message, const int from) const;
		template<class Type> int Receive_col(Type* array, const int len, const int from) const;
		//SEND_UP ===========================================
		/*
		processes with grid_rank_[0/1] = N will send a message to grid_rank_[0/1] = N+1, with torus topology.
		SendUp_row sends a message in the same row.
		SendUp_col sends a message in the same column.
		/ sendBuffer: message/array to send.
		/ recvBuffer: variable to receive message/array.
		*/
		template<class Type> int SendUp_row(const Type* send_buffer, Type* recv_buffer, const int len) const;
		template<class Type> int SendUp_col(const Type* send_buffer, Type* recv_buffer, const int len) const;
		//SEND_DOWN =========================================
		/*
		processes with grid_rank_[0/1] = N will send a message to grid_rank_[0/1] = N-1, with torus topology.
		SendDown_row sends a message in the same row.
		SendDown_col sends a message in the same column.
		/ sendBuffer: message/array to send.
		/ recvBuffer: variable to receive message/array.
		*/
		template<class Type> int SendDown_row(const Type* send_buffer, Type* recv_buffer, const int len) const;
		template<class Type> int SendDown_col(const Type* send_buffer, Type* recv_buffer, const int len) const;
		//SEND_UPDOWN =======================================
		/*
		processes with grid_rank_[0/1] = N will send a message to grid_rank_[0/1] = N-1 and grid_rank_[0/1] = N+1,
		with torus topology.
		*/
		template<class Type> int SendUpDown_row(const Type* buffer_sendup,
												Type& buffer_recvup, const int len_up, 
												const Type* buffer_senddown,
												Type& buffer_recvdown, const int len_down ) const;
		template<class Type> int SendUpDown_col(const Type* buffer_sendup,
												Type& buffer_recvup, const int len_up, 
												const Type* buffer_senddown,
												Type& buffer_recvdown, const int len_down ) const;

		//TRANSFORMS ========================================
		/*
		transform between world rank and 2D grid rank.
		*/
		int* rank_world2grid(const int world_rank, int* grid_rank) const {
			grid_rank[0] = world_rank % this->grid_size_[1];
			grid_rank[1] = world_rank / this->grid_size_[1];
			return grid_rank;
		}
		int rank_grid2world(const int* grid_rank) const {
			return grid_rank[0] * this->stride_[1] + grid_rank[1] * this->stride_[0];
		}
		//MISCELLANEOUS =====================================
		int get_world_size() const { return this->world_size_;}
		int get_world_rank() const { return this->world_rank_;}
		const int* get_grid_size() const { return this->grid_size_;}
		const int* get_grid_rank() const { return this->grid_rank_;}
		const int* get_stride() const { return this->stride_;}
		int get_root() const { return this->root_;}
		bool is_root() const { return this->is_root_;}
		const MPI_Comm& get_world_comm() const { return this->world_comm_;}
		const MPI_Comm& get_row_comm() const { return this->row_comm_;}
		const MPI_Comm& get_col_comm() const { return this->col_comm_;}
		const MPI_Group& get_world_group() const { return this->world_group_;}
		const MPI_Group& get_row_group() const { return this->row_group_;}
		const MPI_Group& get_col_group() const { return this->col_group_;}

	private:
		//Number of processes
		int world_size_;
		//Process ID
		int world_rank_;
		//Number of processes for row(dim0) and col(dim1)
		int grid_size_[2];
		int& row_size_ = grid_size_[0];
		int& col_size_ = grid_size_[1];
		//Process ID in the 2D grid
		/*
		/grid_rank_[0] give the column number;
		/grid_rank_[1] gives the row number.
		*/
		int grid_rank_[2];
		int& row_number_ = grid_rank_[0];
		int& col_number_ = grid_rank_[1];

		//stride of the 2D process array
		//last index changes first.
		int stride_[2];

		int root_;
		bool is_root_;
		// bool last_proc_[2];

		// world_comm_: communicator that controls all the processes, including IO processes.
		// row_comm_: communicator array with size grid_size_[0].
		// col_comm_: communicator array with size grid_size_[1].
		MPI_Comm world_comm_, row_comm_, col_comm_;
		MPI_Group world_group_, row_group_, col_group_, corner_group_;

		// internal functions ==================
		/*
		MPI_Allreduce wrapper
		/TODO: should test if the second one is faster.
		*/
		template<class Type>
		int _allreduce_(Type* array, const int len, const MPI_Op& op, const MPI_Comm& comm) const;
		template<class Type, int LEN>
		int _allreduce_(Type array[LEN], const MPI_Op& op, const MPI_Comm& comm) const;
		
		/*
		row and columns advance in a 2D grid.
		/steps: advance steps, can be negative. Assume the grid has torus topology.
		/ref_rank: the reference process rank of the advance. can be world_rank, 
			or 2D grid rank.
		*/
		int _row_advance_(const int steps, const int ref_rank) const;
		int _row_advance_(const int steps, const int* ref_rank) const;
		int _col_advance_(const int steps, const int ref_rank) const;
		int _col_advance_(const int steps, const int* ref_rank) const;
	};

	template<class Type>
	int Parallel2D::Broadcast(Type& message, const int from) const {
		return MPI_Bcast(&message, sizeof(Type), MPI_BYTE, from, this->world_comm_);
	}

	template<class Type>
	int Parallel2D::Broadcast(Type* array, const int len, const int from) const {
		return MPI_Bcast(array, len*sizeof(Type), MPI_BYTE, from, this->world_comm_);
	}

	template<class Type>
	int Parallel2D::Broadcast_row(Type& message, const int from) const{
		return MPI_Bcast(&message, sizeof(Type), MPI_BYTE, from, this->row_comm_);
	}

	template<class Type>
	int Parallel2D::Broadcast_row(Type* array, const int len, const int from) const{
		return MPI_Bcast(array, len*sizeof(Type), MPI_BYTE, from, this->row_comm_);
	}

	template<class Type>
	int Parallel2D::Broadcast_col(Type& message, const int from) const{
		return MPI_Bcast(&message, sizeof(Type), MPI_BYTE, from, this->col_comm_);
	}

	template<class Type>
	int Parallel2D::Broadcast_col(Type* array, const int len, const int from) const{
		return MPI_Bcast(array, len*sizeof(Type), MPI_BYTE, from, this->col_comm_);
	}

	template<class Type>
	int Parallel2D::Gather(const Type* send_array, const int len,
					Type* recv_array, const int gather_to) const {
		return MPI_Gather(send_array, len, 
						get_MPI_Datatype<Type>(), 
						recv_array, len,
						get_MPI_Datatype<Type>(), gather_to, this->world_comm_);
	}

	template<class Type>
	int Parallel2D::Scatter(const Type* send_array, const int len,
					Type* recv_array, const int scatter_from) const {
		return MPI_Scatter(send_array, len, 
						get_MPI_Datatype<Type>(),
						recv_array, len, 
						get_MPI_Datatype<Type>(), 
						scatter_from, this->world_comm_);
	}

	template<class Type>
	int Parallel2D::_allreduce_(Type* array, 
								const int len, 
								const MPI_Op& op, 
								const MPI_Comm& comm) const {
		auto* rec_buffer = new Type[len];
		auto status = MPI_Allreduce(array, rec_buffer, 
			len, get_MPI_Datatype<Type>(), op, comm); 
		auto rec_end = rec_buffer + len;
		std::copy(rec_buffer, rec_end, array);
		delete[] rec_buffer;
		return status;
	}

	template<class Type, int LEN>
	int Parallel2D::_allreduce_(Type array[LEN],
								const MPI_Op& op,
								const MPI_Comm& comm) const {
		Type rec_buffer[LEN];
		auto status = MPI_Allreduce(array, rec_buffer, LEN, 
			get_MPI_Datatype<Type>(), op, comm);
		auto rec_end = rec_buffer + LEN;
		std::copy(rec_buffer, rec_end, array);
		return status;
	}

	template<class Type>
	int Parallel2D::Reduce(Type& message, const MPI_Op& op) const {
		return this->Reduce(&message, 1, op);
	}

	template<class Type>
	int Parallel2D::Reduce(Type* array, const int len, const MPI_Op& op) const {
		return this->_allreduce_(array, len, op, this->world_comm_);
	}

	template<class Type>
	int Parallel2D::Reduce_row(Type& message, const MPI_Op& op) const {
		return this->Reduce_row(&message, 1, op);
	}

	template<class Type>
	int Parallel2D::Reduce_row(Type* array, const int len, const MPI_Op& op) const {
		return this->_allreduce_(array, len, op, this->row_comm_);
	}

	template<class Type>
	int Parallel2D::Reduce_col(Type& message, const MPI_Op& op) const {
		return this->Reduce_col(&message, 1, op);
	} 

	template<class Type>
	int Parallel2D::Reduce_col(Type* array, const int len, const MPI_Op& op) const {
		return this->_allreduce_(array, len, op, this->col_comm_);
	}

	template<class Type>
	int Parallel2D::Sum(Type& message) const {
		return this->Reduce(message, MPI_SUM);
	}

	template<class Type>
	int Parallel2D::Sum(Type* array, const int len) const {
		return this->Reduce(array, len, MPI_SUM);
	}

	template<class Type>
	int Parallel2D::Sum_row(Type& message) const {
		return this->Reduce_row(message, MPI_SUM);
	}

	template<class Type>
	int Parallel2D::Sum_row(Type* array, const int len) const {
		return this->Reduce_row(array, len, MPI_SUM);
	}

	template<class Type>
	int Parallel2D::Sum_col(Type& message) const {
		return this->Reduce_col(message, MPI_SUM);
	}

	template<class Type>
	int Parallel2D::Sum_col(Type* array, const int len) const {
		return this->Reduce_col(array, len, MPI_SUM);
	}

	template<class Type>
	int Parallel2D::Max(Type& message) const {
		return this->Reduce(message, MPI_MAX);
	}

	template<class Type>
	int Parallel2D::Max(Type* array, const int len) const {
		return this->Reduce(array, len, MPI_MAX);
	}

	template<class Type>
	int Parallel2D::Max_row(Type& message) const {
		return this->Reduce_row(message, MPI_MAX);
	}

	template<class Type>
	int Parallel2D::Max_row(Type* array, const int len) const {
		return this->Reduce_row(array, len, MPI_MAX);
	}

	template<class Type>
	int Parallel2D::Max_col(Type& message) const {
		return this->Reduce_col(message, MPI_MAX);
	}

	template<class Type>
	int Parallel2D::Max_col(Type* array, const int len) const { 
		return this->Reduce_col(array, len, MPI_MAX);
	}

	template<class Type>
	int Parallel2D::Min(Type& message) const {
		return this->Reduce(message, MPI_MIN);
	}

	template<class Type>
	int Parallel2D::Min(Type* array, const int len) const {
		return this->Reduce(array, len, MPI_MIN);
	}

	template<class Type>
	int Parallel2D::Min_row(Type& message) const {
		return this->Reduce_row(message, MPI_MIN);
	}

	template<class Type>
	int Parallel2D::Min_row(Type* array, const int len) const {
		return this->Reduce_row(array, len, MPI_MIN);
	}

	template<class Type>
	int Parallel2D::Min_col(Type& message) const {
		return this->Reduce_col(message, MPI_MIN);
	}

	template<class Type>
	int Parallel2D::Min_col(Type* array, const int len) const {
		return this->Reduce_col(array, len, MPI_MIN);
	}

	template<class Type>
	int Parallel2D::Send(Type& message, const int to) const {
		return MPI_Send(&message, sizeof(Type), MPI_BYTE, to, 0, this->world_comm_);
	}

	template<class Type>
	int Parallel2D::Send(Type* array, const int len, const int to) const {
		return MPI_Send(array, len*sizeof(Type), MPI_BYTE, to, 0, this->world_comm_);
	}

	template<class Type> 
	int Parallel2D::Send(Type* array, const int len, const int* to_gridrank) const{
		return MPI_Send(array, len*sizeof(Type), MPI_BYTE, 
						rank_grid2world(to_gridrank), 0, this->world_comm_);
	}

	template<class Type>
	int Parallel2D::Send_row(Type& message, const int to) const {
		return MPI_Send(&message, sizeof(Type), MPI_BYTE, to, 0, this->row_comm_);
	}

	template<class Type>
	int Parallel2D::Send_row(Type* array, const int len, const int to) const {
		return MPI_Send(array, len*sizeof(Type), MPI_BYTE, to, 0, this->row_comm_);
	}

	template<class Type>
	int Parallel2D::Send_col(Type& message, const int to) const {
		return MPI_Send(&message, sizeof(Type), MPI_BYTE, to, 0, this->col_comm_);
	}

	template<class Type>
	int Parallel2D::Send_col(Type* array, const int len, const int to) const {
		return MPI_Send(array, len*sizeof(Type), MPI_BYTE, to, 0, this->col_comm_);
	}

	template<class Type>
	int Parallel2D::Receive(Type& message, const int from) const {
		return MPI_Recv(&message, sizeof(Type), MPI_BYTE, from, 0, this->world_comm_, MPI_STATUS_IGNORE);
	}

	template<class Type>
	int Parallel2D::Receive(Type* array, const int len, const int from) const {
		return MPI_Recv(array, len*sizeof(Type), MPI_BYTE, from, 0, this->world_comm_, MPI_STATUS_IGNORE);
	}

	template<class Type>
	int Parallel2D::Receive(Type* array, const int len, const int* from_gridrank) const{
		return MPI_Recv(array, len*sizeof(Type), MPI_BYTE, 
						rank_grid2world(from_gridrank), 0, this->world_comm_, MPI_STATUS_IGNORE);
	}

	template<class Type>
	int Parallel2D::Receive_row(Type& message, const int from) const {
		return MPI_Recv(&message, sizeof(Type), MPI_BYTE, from, 0, this->row_comm_, MPI_STATUS_IGNORE);
	}

	template<class Type>
	int Parallel2D::Receive_row(Type* array, const int len, const int from) const {
		return MPI_Recv(array, len*sizeof(Type), MPI_BYTE, from, 0, this->row_comm_, MPI_STATUS_IGNORE);
	}

	template<class Type>
	int Parallel2D::Receive_col(Type& message, const int from) const {
		return MPI_Recv(&message, sizeof(Type), MPI_BYTE, from, 0, this->col_comm_, MPI_STATUS_IGNORE);
	}

	template<class Type>
	int Parallel2D::Receive_col(Type* array, const int len, const int from) const { 
		return MPI_Recv(array, len*sizeof(Type), MPI_BYTE, from, 0, this->col_comm_, MPI_STATUS_IGNORE);
	}

	template<class Type>
	int Parallel2D::SendUp_row(const Type* send_buffer, Type* recv_buffer, const int len) const {
		auto& my_row_rank = this->grid_rank_[0]; //aka. my_col_number
		auto send_dest_row_rank = (my_row_rank + 1) % this->grid_size_[1];
		auto recv_from_row_rank = ( (my_row_rank - 1) % this->grid_size_[1] 
								+ this->grid_size_[1] ) % this->grid_size_[1];
		//auto send_status = MPI_Send(send_buffer, len*sizeof(Type), MPI_BYTE, 
		//							send_dest_row_rank, 0, this->row_comm_);
		//auto recv_status = MPI_Recv(recv_buffer, len*sizeof(Type), MPI_BYTE, 
		//							recv_from_row_rank, 0, this->row_comm_, 
		//							MPI_STATUS_IGNORE);
		auto status = MPI_Sendrecv(send_buffer, len*sizeof(Type), 
									MPI_BYTE, send_dest_row_rank, 0, 
									recv_buffer, len*sizeof(Type), 
									MPI_BYTE, recv_from_row_rank, 0, this->row_comm_, 
									MPI_STATUS_IGNORE);
		return status;
	}

	template<class Type>
	int Parallel2D::SendUp_col(const Type* send_buffer, Type* recv_buffer, const int len) const {
		auto& my_col_rank = this->grid_rank_[1]; //aka. my_row_number
		auto send_dest_col_rank = (my_col_rank + 1) % this->grid_size_[0];
		auto recv_from_col_rank = ( (my_col_rank - 1) % this->grid_size_[0] + this->grid_size_[0] )
		 						% this->grid_size_[0];
		//auto send_status = MPI_Send(send_buffer, len*sizeof(Type), MPI_BYTE, 
		//							send_dest_col_rank, 0, this->col_comm_);
		//auto recv_status = MPI_Recv(recv_buffer, len*sizeof(Type), MPI_BYTE, 
		//							recv_from_col_rank, 0, this->col_comm_,
		//							MPI_STATUS_IGNORE);
		auto status = MPI_Sendrecv(send_buffer, len*sizeof(Type), 
									MPI_BYTE, send_dest_col_rank, 0, 
									recv_buffer, len*sizeof(Type), 
									MPI_BYTE, recv_from_col_rank, 0, this->col_comm_, 
									MPI_STATUS_IGNORE);
		return status;
	}

	template<class Type>
	int Parallel2D::SendDown_row(const Type* send_buffer, Type* recv_buffer, const int len) const {
		auto& my_row_rank = this->grid_rank_[0]; //aka. my_col_number
		auto send_dest_row_rank = ( (my_row_rank - 1) % this->grid_size_[1] 
								+ this->grid_size_[1] ) % this->grid_size_[1];
		auto recv_from_row_rank = (my_row_rank + 1) % this->grid_size_[1];
		//auto send_status = MPI_Send(send_buffer, len*sizeof(Type), MPI_BYTE, 
		//							send_dest_row_rank, 0, this->row_comm_);
		//auto recv_status = MPI_Recv(recv_buffer, len*sizeof(Type), MPI_BYTE, 
		//							recv_from_row_rank, 0, this->row_comm_,
		//							MPI_STATUS_IGNORE);
		auto status = MPI_Sendrecv(send_buffer, len*sizeof(Type), 
									MPI_BYTE, send_dest_row_rank, 0, 
									recv_buffer, len*sizeof(Type), 
									MPI_BYTE, recv_from_row_rank, 0, this->row_comm_, 
									MPI_STATUS_IGNORE);
		return status;
	}

	template<class Type>
	int Parallel2D::SendDown_col(const Type* send_buffer, Type* recv_buffer, const int len) const {
		auto& my_col_rank = this->grid_rank_[1]; //aka. my_row_number
		auto send_dest_col_rank = ( (my_col_rank - 1) % this->grid_size_[0] + this->grid_size_[0] )
		 						% this->grid_size_[0];
		auto recv_from_col_rank = (my_col_rank + 1) % this->grid_size_[0];
		//auto send_status = MPI_Send(send_buffer, len*sizeof(Type), MPI_BYTE, 
		//							send_dest_col_rank, 0, this->col_comm_);
		//auto recv_status = MPI_Recv(recv_buffer, len*sizeof(Type), MPI_BYTE, 
		//							recv_from_col_rank, 0, this->col_comm_,
		//							MPI_STATUS_IGNORE);
		auto status = MPI_Sendrecv(send_buffer, len*sizeof(Type), 
									MPI_BYTE, send_dest_col_rank, 0, 
									recv_buffer, len*sizeof(Type), 
									MPI_BYTE, recv_from_col_rank, 0, this->col_comm_, 
									MPI_STATUS_IGNORE);
		return status;
	}
}

#endif
