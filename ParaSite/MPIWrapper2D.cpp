#include"MPIWrapper2D.h"
#include<iostream>
#include<cassert>


namespace MPI_Wrapper{

    int Parallel_Init(int* argc, char*** argv){
        return MPI_Init(argc, argv);
    }

    int Parallel_Finalize(){
        return MPI_Finalize();
    }

	template<>
	MPI_Datatype get_MPI_Datatype<short>(){
		return MPI_SHORT;
	}
	template<>
	MPI_Datatype get_MPI_Datatype<unsigned short>(){
		return MPI_UNSIGNED_SHORT;
	}
	template<>
	MPI_Datatype get_MPI_Datatype<int>(){
		return MPI_INT;
	}
	template<>
	MPI_Datatype get_MPI_Datatype<unsigned int>(){
		return MPI_UNSIGNED;
	}
	template<>
	MPI_Datatype get_MPI_Datatype<long>(){
		return MPI_LONG;
	}
	template<>
	MPI_Datatype get_MPI_Datatype<long long>(){
		return MPI_LONG_LONG;
	}
	template<>
	MPI_Datatype get_MPI_Datatype<float>(){
		return MPI_FLOAT;
	}
	template<>
	MPI_Datatype get_MPI_Datatype<double>(){
		return MPI_DOUBLE;
	}
	template<>
	MPI_Datatype get_MPI_Datatype<char>(){
		return MPI_CHAR;
	}


    Parallel2D::Parallel2D(MPI_Comm world_comm)
                            :
                            grid_size_{0,0},
                            grid_rank_{0,0},
                            stride_{0,0},
                            row_comm_(MPI_COMM_NULL),
                            col_comm_(MPI_COMM_NULL),
                            row_group_(MPI_GROUP_NULL),
                            col_group_(MPI_GROUP_NULL),
                            corner_group_(MPI_GROUP_NULL){
		auto status = MPI_Comm_dup( world_comm, &(this->world_comm_) );
        assert(status == MPI_SUCCESS);
        MPI_Comm_rank(this->world_comm_, &(this->world_rank_));
        MPI_Comm_size(this->world_comm_, &(this->world_size_));
        MPI_Comm_group(this->world_comm_, &(this->world_group_));
        this->root_ = 0;
        if( this->world_rank_ == this->root_ ) this->is_root_ = true;
        else this->is_root_ = false; 
    }

	Parallel2D::Parallel2D(const int n_rows,
							const int n_cols,
							MPI_Comm world_comm)
		:
		Parallel2D(world_comm){
		this->InitializeGrid(n_rows, n_cols);
	}


    Parallel2D::~Parallel2D(){
        MPI_Group_free(&(this->world_group_));
        if(this->row_group_ != MPI_GROUP_NULL) MPI_Group_free(&(this->row_group_));
        if(this->col_group_ != MPI_GROUP_NULL) MPI_Group_free(&(this->col_group_));
        MPI_Comm_free(&(this->world_comm_));
        if(this->row_comm_ != MPI_COMM_NULL) MPI_Comm_free(&(this->row_comm_));
        if(this->col_comm_ != MPI_COMM_NULL) MPI_Comm_free(&(this->col_comm_));
    }

    void Parallel2D::InitializeGrid(const int row_size, const int col_size){
        this->grid_size_[0] = row_size;
        this->grid_size_[1] = col_size;
        this->stride_[0] = this->grid_size_[1];
        this->stride_[1] = 1;
        //assign row_comm_
        auto color = this->world_rank_ / this->grid_size_[1]; //row number
        MPI_Comm_split(this->world_comm_, color, this->world_rank_, &(this->row_comm_));
        MPI_Comm_rank(this->row_comm_, &(this->grid_rank_[0]));
        MPI_Comm_group(this->row_comm_, &(this->row_group_));
        //assign col_comm_
        color = this->world_rank_ % this->grid_size_[1]; //col number
        MPI_Comm_split(this->world_comm_, color, this->world_rank_, &(this->col_comm_));
        MPI_Comm_rank(this->col_comm_, &(this->grid_rank_[1]));
        MPI_Comm_group(this->col_comm_, &(this->col_group_));
        //assgin corner_group_
        
        return;
    }

    int Parallel2D::ForceAbort(){
        return MPI_Abort(this->world_comm_, EXIT_FAILURE);
    }

    int Parallel2D::Barrier() const {
        return MPI_Barrier(this->world_comm_);
    }

	// internal functions ===============================
	int Parallel2D::_row_advance_(const int steps, const int ref_rank) const {
		int grid_rank[2] = { 0, 0 };
		return _row_advance_(steps, rank_world2grid(ref_rank, grid_rank));
	}

	int Parallel2D::_row_advance_(const int steps, const int* ref_rank) const {
		// column number (ref_rank[0]) is unchanged.
		int advanced_rank[] = { ref_rank[0], 
			((ref_rank[1] + steps) % row_size_ + row_size_) % row_size_ };
		return rank_grid2world(advanced_rank);
	}

	int Parallel2D::_col_advance_(const int steps, const int ref_rank) const {
		int grid_rank[2] = { 0,0 };
		return _col_advance_(steps, rank_world2grid(ref_rank, grid_rank));
	}

	int Parallel2D::_col_advance_(const int steps, const int* ref_rank) const {
		// row number (ref_rank[1]) is unchanged. 
		int advanced_rank[] = { ((ref_rank[0] + steps) % col_size_ + col_size_) % col_size_, 
			ref_rank[1] };
		return rank_grid2world(advanced_rank);
	}

//No longer needed once we defined get_MPI_datatype().
/*
    // _allreduce_ specializations
    template<>
	int Parallel2D::_allreduce_<unsigned int>(unsigned int* array, 
								const int len, 
								const MPI_Op& op, 
								const MPI_Comm& comm) const {
		auto* rec_buffer = new unsigned int[len];
		auto status = MPI_Allreduce(array, rec_buffer, 
			len, MPI_UNSIGNED, op, comm); 
		auto rec_end = rec_buffer + len;
		std::copy(rec_buffer, rec_end, array);
		delete[] rec_buffer;
		return status;
	}

    template<>
	int Parallel2D::_allreduce_<int>(int* array, 
								const int len, 
								const MPI_Op& op, 
								const MPI_Comm& comm) const {
		auto* rec_buffer = new int[len];
		auto status = MPI_Allreduce(array, rec_buffer, 
			len, MPI_INT, op, comm); 
		auto rec_end = rec_buffer + len;
		std::copy(rec_buffer, rec_end, array);
		delete[] rec_buffer;
		return status;
	}

    template<>
	int Parallel2D::_allreduce_<float>(float* array, 
								const int len, 
								const MPI_Op& op, 
								const MPI_Comm& comm) const {
		auto* rec_buffer = new float[len];
		auto status = MPI_Allreduce(array, rec_buffer, 
			len, MPI_FLOAT, op, comm); 
		auto rec_end = rec_buffer + len;
		std::copy(rec_buffer, rec_end, array);
		delete[] rec_buffer;
		return status;
	}

	template<>
	int Parallel2D::_allreduce_<double>(double* array, 
								const int len, 
								const MPI_Op& op, 
								const MPI_Comm& comm) const {
		auto* rec_buffer = new double[len];
		auto status = MPI_Allreduce(array, rec_buffer, 
			len, MPI_DOUBLE, op, comm); 
		auto rec_end = rec_buffer + len;
		std::copy(rec_buffer, rec_end, array);
		delete[] rec_buffer;
		return status;
	}
*/
}