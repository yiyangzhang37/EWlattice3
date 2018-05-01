#ifndef MPI_WRAPPER_2D
#define MPI_WRAPPER_2D

#include "mpi.h"


class Parallel2D {
public:
	//CONSTRUCTOR AND DESTRUCTOR ========================
	Parallel2D();
	~Parallel2D();

	//INITIALIZATION ====================================
	/*
	\ proc_size0 : size of the first dimension of the MPI process grid.
	\ proc_size1 : size of the second dimension of the MPI process grid.
	*/
	void Initialize(const int proc_size0, const int proc_size1);

	//ABORT =============================================
	void ForceAbort();
	void AbortRequest();

	//BARRIER ===========================================
	void Barrier();

	//GLOBAL AND DIRECTIONAL PROCESSES COMMUNICATIONS ===

	//BROADCAST =========================================
	/*
	broadcast a message from process with world_rank_ == from to all the other processes.
	*/
	template<class Type> int Broadcast(Type& message, const int from);

	/*
	broadcast an array from process with world_rank_ == from to all the other processes.
	*/
	template<class Type> int Broadcast(Type* array, const int len, const int from);

	/*
	broadcast message/array in the dim1 direction. 
	Processes with grid_rank_[0] == from will broadcast the variable to processes 
	with same grid_rank_[1].
	*/
	template<class Type> int Broadcast_dim0(Type& message, const int from);
	template<class Type> int Broadcast_dim0(Type* array, const int len, const int from);

	/*
	broadcast message/array in the dim0 direction.
	Processes with grid_rank_[1] == from will broadcast the variable to processes
	with same grid_rank_[0].
	*/
	template<class Type> int Broadcast_dim1(Type& message, const int from);
	template<class Type> int Broadcast_dim1(Type* array, const int len, const int from);
	
	//REDUCE ============================================
	template<class Type> int Reduce(Type& message, const MPI_Op& op);
	template<class Type> int Reduce(Type& array, const int len, const MPI_Op& op);
	template<class Type> int Reduce_dim0(Type& message, const MPI_Op& op);
	template<class Type> int Reduce_dim0(Type& array, const int len, const MPI_Op& op);
	template<class Type> int Reduce_dim1(Type& message, const MPI_Op& op);
	template<class Type> int Reduce_dim1(Type& array, const int len, const MPI_Op& op);

	//SUM ===============================================
	template<class Type> int Sum(Type& message);
	template<class Type> int Sum(Type& array, const int len);
	template<class Type> int Sum_dim0(Type& message);
	template<class Type> int Sum_dim0(Type& array, const int len);
	template<class Type> int Sum_dim1(Type& message);
	template<class Type> int Sum_dim1(Type& array, const int len);

	//MAX ===============================================
	template<class Type> int Max(Type& message);
	template<class Type> int Max(Type& array, const int len);
	template<class Type> int Max_dim0(Type& message);
	template<class Type> int Max_dim0(Type& array, const int len);
	template<class Type> int Max_dim1(Type& message);
	template<class Type> int Max_dim1(Type& array, const int len);

	//MIN ===============================================
	template<class Type> int Min(Type& message);
	template<class Type> int Min(Type& array, const int len);
	template<class Type> int Min_dim0(Type& message);
	template<class Type> int Min_dim0(Type& array, const int len);
	template<class Type> int Min_dim1(Type& message);
	template<class Type> int Min_dim1(Type& array, const int len);

	//SEND ==============================================
	/*
	send a message/an array to the process (receiver) with world_rank_ == to.
	*/
	template<class Type> int Send(Type& message, const int to);
	template<class Type> int Send(Type* array, const int len, const int to);
	/*
	*/
	template<class Type> int Send_dim0(Type& message, const int to);
	//RECEIVE ===========================================
	//SEND_UP ===========================================
	//SEND_DOWN =========================================
	//SEND_UPDOWN =======================================

	//MISCELLANEOUS =====================================
private:
	//Number of processes
	int world_size_;
	//Process ID
	int world_rank_;
	//Number of processes for dim0 and dim1
	int grid_size_[2];
	//Process ID in the 2D grid
	int grid_rank_[2];

	int root_;
	bool isRoot_;
	bool last_proc_[2];

	// world_comm_: communicator that controls all the processes, including IO processes.
	// dim0_comm_: communicator array with size grid_size_[0].
	// dim1_comm_: communicator array with size grid_size_[1].
	MPI_Comm world_comm_, *dim0_comm_, *dim1_comm_;
	MPI_Group world_group_, *dim0_group_, *dim1_group_;
};

#endif
