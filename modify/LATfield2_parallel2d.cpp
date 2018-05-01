#include "LATfield2_parallel2d.h"

Parallel2d::Parallel2d()
{


	int argc = 1;
	char** argv = new char*[argc];
	for (int i = 0; i<argc; i++) { argv[i] = new char[20]; }
#ifndef EXTERNAL_IO	
	lat_world_comm_ = MPI_COMM_WORLD;
	world_comm_ = MPI_COMM_WORLD;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(lat_world_comm_, &lat_world_rank_);
	MPI_Comm_size(lat_world_comm_, &lat_world_size_);
	MPI_Comm_rank(world_comm_, &world_rank_);
	MPI_Comm_size(world_comm_, &world_size_);
#else
	world_comm_ = MPI_COMM_WORLD;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(world_comm_, &world_rank_);
	MPI_Comm_size(world_comm_, &world_size_);

#endif

}
#ifdef EXTERNAL_IO
void Parallel2d::initialize(int proc_size0, int proc_size1, int IO_total_size, int IO_node_size)
#else
void Parallel2d::initialize(int proc_size0, int proc_size1)
#endif
{

	grid_size_[0] = proc_size0;
	grid_size_[1] = proc_size1;

	dim0_comm_ = (MPI_Comm *)malloc(grid_size_[1] * sizeof(MPI_Comm));
	dim1_comm_ = (MPI_Comm *)malloc(grid_size_[0] * sizeof(MPI_Comm));

	dim0_group_ = (MPI_Group *)malloc(grid_size_[1] * sizeof(MPI_Group));
	dim1_group_ = (MPI_Group *)malloc(grid_size_[0] * sizeof(MPI_Group));

	int rang[3], i, j, comm_rank;

#ifdef EXTERNAL_IO

	if (world_rank_ == 0)
	{
		if (proc_size0*proc_size1 + IO_total_size != world_size_)
		{
			cerr << "Latfield::Parallel2d::initialization - wrong number of process" << endl;
			cerr << "Latfield::Parallel2d::initialization - Number of total process must be equal to proc_size0*proc_size1+IO_total_size" << endl;
			cerr << "Latfield::Parallel2d::initialization - Within the call : Parallel2d(int proc_size0, int proc_size1, int IO_total_size)" << endl;
			this->abortForce();
		}



	}

	MPI_Comm_group(world_comm_, &world_group_);

	rang[0] = 0;
	rang[1] = proc_size0 * proc_size1 - 1;
	rang[2] = 1;

	MPI_Group_range_incl(world_group_, 1, &rang, &lat_world_group_);
	MPI_Comm_create(world_comm_, lat_world_group_, &lat_world_comm_);

	//lat_world_comm_ = MPI_COMM_WORLD;
	//MPI_Comm_group(lat_world_comm_,&lat_world_group_);


	MPI_Group_rank(lat_world_group_, &comm_rank);
	if (comm_rank != MPI_UNDEFINED)
	{
		lat_world_rank_ = comm_rank;
		MPI_Comm_size(lat_world_comm_, &lat_world_size_);

		rang[2] = 1;
		for (j = 0; j<grid_size_[1]; j++)
		{
			rang[0] = j * grid_size_[0];
			rang[1] = grid_size_[0] - 1 + j * grid_size_[0];
			MPI_Group_range_incl(lat_world_group_, 1, &rang, &dim0_group_[j]);
			MPI_Comm_create(lat_world_comm_, dim0_group_[j], &dim0_comm_[j]);
		}


		rang[2] = grid_size_[0];
		for (i = 0; i<grid_size_[0]; i++)
		{
			rang[0] = i;
			rang[1] = i + (grid_size_[1] - 1)*grid_size_[0];
			MPI_Group_range_incl(lat_world_group_, 1, &rang, &dim1_group_[i]);
			MPI_Comm_create(lat_world_comm_, dim1_group_[i], &dim1_comm_[i]);
		}


		for (i = 0; i<grid_size_[0]; i++)
		{
			MPI_Group_rank(dim1_group_[i], &comm_rank);
			if (comm_rank != MPI_UNDEFINED)grid_rank_[1] = comm_rank;
		}

		for (j = 0; j<grid_size_[1]; j++)
		{
			MPI_Group_rank(dim0_group_[j], &comm_rank);
			if (comm_rank != MPI_UNDEFINED)grid_rank_[0] = comm_rank;
		}


		root_ = 0;
		isIO_ = false;
	}
	else
	{
		lat_world_rank_ = -1;
		root_ = 0;
		grid_rank_[1] = -1;
		grid_rank_[0] = -1;
		isIO_ = true;
	}

	if (grid_rank_[0] == grid_size_[0] - 1)last_proc_[0] = true;
	else last_proc_[0] = false;
	if (grid_rank_[1] == grid_size_[1] - 1)last_proc_[1] = true;
	else last_proc_[1] = false;




	IO_Server.initialize(proc_size0, proc_size1, IO_total_size, IO_node_size);

#else

	if (lat_world_rank_ == 0)
	{
		if (proc_size0*proc_size1 != lat_world_size_)
		{
			std::cerr << "Latfield::Parallel2d::initialization - wrong number of process" << std::endl;
			std::cerr << "Latfield::Parallel2d::initialization - Number of total process must be equal to proc_size0*proc_size1" << std::endl;
			std::cerr << "Latfield::Parallel2d::initialization - Within the call : Parallel2d(int proc_size0, int proc_size1)" << std::endl;
			this->abortForce();
		}
	}

	MPI_Comm_group(lat_world_comm_, &lat_world_group_);

	rang[2] = 1;
	for (j = 0; j<grid_size_[1]; j++)
	{
		rang[0] = j * grid_size_[0];
		rang[1] = grid_size_[0] - 1 + j * grid_size_[0];
		MPI_Group_range_incl(lat_world_group_, 1, &rang, &dim0_group_[j]);
		MPI_Comm_create(lat_world_comm_, dim0_group_[j], &dim0_comm_[j]);
	}

	rang[2] = grid_size_[0];
	for (i = 0; i<grid_size_[0]; i++)
	{
		rang[0] = i;
		rang[1] = i + (grid_size_[1] - 1)*grid_size_[0];
		MPI_Group_range_incl(lat_world_group_, 1, &rang, &dim1_group_[i]);
		MPI_Comm_create(lat_world_comm_, dim1_group_[i], &dim1_comm_[i]);
	}



	for (i = 0; i<grid_size_[0]; i++)
	{
		MPI_Group_rank(dim1_group_[i], &comm_rank);
		if (comm_rank != MPI_UNDEFINED)grid_rank_[1] = comm_rank;
	}

	for (j = 0; j<grid_size_[1]; j++)
	{
		MPI_Group_rank(dim0_group_[j], &comm_rank);
		if (comm_rank != MPI_UNDEFINED)grid_rank_[0] = comm_rank;
	}


	root_ = 0;

	if (grid_rank_[0] == grid_size_[0] - 1)last_proc_[0] = true;
	else last_proc_[0] = false;
	if (grid_rank_[1] == grid_size_[1] - 1)last_proc_[1] = true;
	else last_proc_[1] = false;

#endif

	if (root_ == lat_world_rank_)isRoot_ = true;
	else isRoot_ = false;

}

Parallel2d::~Parallel2d()
{


	/*free(dim0_comm_);
	free(dim1_comm_);
	free(dim0_group_);
	free(dim1_group_);*/
	int finalized;
	MPI_Finalized(&finalized);
	if (!finalized) { MPI_Finalize(); }
}

//ABORT AND BARRIER===============================

void Parallel2d::abortForce()
{
	MPI_Abort(world_comm_, EXIT_FAILURE);
}

void Parallel2d::abortRequest()
{
	char failure;
	if (isRoot())
	{
		failure = char(1);
		broadcast(failure, root_);
		exit(EXIT_FAILURE);
	}
	else
	{
		std::cout << "Parallel::abortRequest() called from non-Root process." << std::endl;
		std::cout << "Parallel::abortRequest() calling abortForce()..." << std::endl;
		abortForce();
	}
}

void Parallel2d::barrier()
{

	MPI_Barrier(lat_world_comm_);
}