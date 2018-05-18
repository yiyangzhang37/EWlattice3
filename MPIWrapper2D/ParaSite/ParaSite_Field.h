#ifndef PARASITE_FIELD_H
#define PARASITE_FIELD_H

#include "ParaSite_Site.h"
#include "MPIWrapper2D.h"
#include "HDF5Wrapper.h"
#include <string>
#include <fstream>
#include <iostream>
#include <memory>
#include <tuple>
#include <type_traits>


namespace ParaSite{

    using ParallelObject = MPI_Wrapper::Parallel2D;
    using FileProcessor = HDF5_Wrapper::HDF5Wrapper;

    template<class FieldType, int DIM>
    class Field{
    private:
        //allocator
        void allocator();
		void deallocator();
        
        //make update_halo_table_
		//This is independent of mat_rows_ and mat_cols_.
        void make_update_halo_table();

        //compute the configurations needed for HDF5 write and save.
        //This function will take care of the convention transition between this 
        //Field class and the HDF5 convention.
        //All the arrrys in this function has length DIM+1.
        void make_fileio_settings(
            hsize_t* dataset_size,
            hsize_t* file_offset,
            hsize_t* file_block_size,
            hsize_t* mem_size,
            hsize_t* mem_offset,
            hsize_t* mem_block_size) const;
    
        //std::unique_ptr<FieldType> data_ptr_ = nullptr;
        FieldType* data_ptr_ = nullptr;

        const Lattice<DIM>* lattice_ = nullptr;
        /*
        stores the shape of the matrix field.
        if it is a vector field, then the components will be stored in mat_rows_,
        while mat_cols will be set to 1.
        /// It will follow row-first convention.
        */
        int components_;
        int mat_rows_;
        int mat_cols_;
        int stride_[2];

        /*
        a pointer pointed to the Parallel2D object.
        */
        const ParallelObject* parallel_ptr;

        /*
        this table saves the information about how to update halos in a parallel configuration.
        The size of the std::vector<> is the number of halo sites in each process.
        The tuple has the information 
        (halo_index, mapped_index, 
        mapped_grid_loc[0] relative to current grid_loc, 
        mapped_grid_loc[1] relative to current grid_loc).
        The index are in terms of local_mem_index.
        */
        std::vector<std::tuple<IndexType, IndexType, int, int>> update_halo_table_;

    public:
        /*
        define a vector field, with dimension = components.
        */
        Field(const Lattice<DIM>& lattice, const int mat_rows);
        Field(const Lattice<DIM>& lattice, const int mat_rows, const int mat_cols);
		Field(const Lattice<DIM>& lattice, const int mat_rows, 
			const ParallelObject& parallel_object);
		Field(const Lattice<DIM>& lattice, const int mat_rows, const int mat_cols,
			const ParallelObject& parallel_object);
        ~Field();

        void assign_parallel_object(const ParallelObject& parallel_object);
		void reinit(const int mat_rows, const int mat_cols);

        //data retrieval and modification
        constexpr int get_component(const int i_row, const int j_col) const{
            return i_row * stride_[0] + j_col * stride_[1];
        }

        FieldType& operator()(const Site<DIM>& site) const;
        FieldType& operator()(const Site<DIM>& site, const int i_row) const;
        FieldType& operator()(const Site<DIM>& site, const int i_row, const int j_col) const;

        FieldType& operator()(const IndexType index) const;
        FieldType& operator()(const IndexType index, const int i_row) const;
        FieldType& operator()(const IndexType index, const int i_row, const int j_col) const;
        
        //file operations
        void read(const std::string& file_name, 
                const std::string& dataset_name = "") const;
        void write(const std::string& file_name,
                const std::string& dataset_name = "") const;

        void update_halo() const;
        
        const Lattice<DIM>& get_lattice() const {return *(this->lattice_);}
        //Lattice<DIM>& get_lattice() const {return const_cast<Lattice<DIM>>(*(this->lattice_));}
        constexpr int get_rows() const {return this->mat_rows_;}
        constexpr int get_cols() const {return this->mat_cols_;}
        constexpr int get_components() const {return this->components_;}

        auto get_update_halo_table() const -> const decltype(this->update_halo_table_)&
            {return this->update_halo_table_;} 

    };

    template<class FieldType, int DIM>
    Field<FieldType, DIM>::Field(const Lattice<DIM>& lattice, const int mat_rows)
        :
        Field(lattice, mat_rows, 1)
    {
    }

    template<class FieldType, int DIM>
    Field<FieldType, DIM>::Field(
        const Lattice<DIM>& lattice, 
        const int mat_rows, 
        const int mat_cols)
        :
        lattice_(&lattice),
        mat_rows_(mat_rows),
        mat_cols_(mat_cols)
    {
        assert(mat_rows_ > 0);
        assert(mat_cols_ > 0);
        components_ = mat_rows_ * mat_cols_;
        stride_[0] = 1;
        stride_[1] = mat_rows;
        allocator();
        make_update_halo_table();
    }

	template<class FieldType, int DIM>
	Field<FieldType, DIM>::Field(
		const Lattice<DIM>& lattice,
		const int mat_rows,
		const ParallelObject& parallel_object)
		:
		Field(lattice, mat_rows, 1, parallel_object)
	{}


	template<class FieldType, int DIM>
	Field<FieldType, DIM>::Field(
		const Lattice<DIM>& lattice,
		const int mat_rows,
		const int mat_cols,
		const ParallelObject& parallel_object)
		:
		Field(lattice, mat_rows, mat_cols)
	{
		this->assign_parallel_object(parallel_object);
	}

    template<class FieldType, int DIM>
    Field<FieldType, DIM>::~Field(){
		this->deallocator();
    }

    template<class FieldType, int DIM>
    void Field<FieldType, DIM>::allocator(){
        //data_ptr_ = std::unique_ptr<FieldType>(
        //    new FieldType[this->lattice_->get_local_mem_sites() * components_] );
        this->data_ptr_ = new FieldType[this->lattice_->get_local_mem_sites() * components_];
		return;
    }

	template<class FieldType, int DIM>
	void Field<FieldType, DIM>::deallocator() {
		if (this->data_ptr_ != nullptr) {
			delete[] this->data_ptr_;
		}
		return;
	}

    template<class FieldType, int DIM>
    void Field<FieldType, DIM>::assign_parallel_object(const ParallelObject& parallel_object){
        this->parallel_ptr = &parallel_object;
    }

	template<class FieldType, int DIM>
	void Field<FieldType, DIM>::reinit(const int mat_rows, const int mat_cols) {
		assert(mat_rows > 0);
		assert(mat_cols > 0);
        this->deallocator();
		this->mat_rows_ = mat_rows;
		this->mat_cols_ = mat_cols;
		components_ = this->mat_rows_ * this->mat_cols_;
		stride_[0] = 1;
		stride_[1] = this->mat_rows_;
		allocator();
		return;
	}

    template<class FieldType, int DIM>
    FieldType& Field<FieldType, DIM>::operator()(const Site<DIM>& site) const {
        return *(data_ptr_ + site.get_index()*this->components_);
    }

    template<class FieldType, int DIM>
    FieldType& Field<FieldType, DIM>::operator()(const Site<DIM>& site, const int i_row) const {
        return *(data_ptr_ + site.get_index()*this->components_ + i_row);
    }

    template<class FieldType, int DIM>
    FieldType& 
    Field<FieldType, DIM>::operator()(const Site<DIM>& site, const int i_row, const int j_col) const {
        return *(data_ptr_ + site.get_index()*this->components_ + get_component(i_row, j_col));
    }

    template<class FieldType, int DIM>
    void Field<FieldType, DIM>::make_update_halo_table() {
        update_halo_table_ = std::vector<std::tuple<IndexType, IndexType, int, int>>(
                            this->lattice_->get_local_mem_sites() - this->lattice_->get_visible_local_sites());
        auto first_halo = this->lattice_->get_local_halo_first();
        auto next_to_last_halo = this->lattice_->get_local_halo_next_to_last();
        auto h = this->lattice_->get_local_halo_first();
        IndexType halo_coord[DIM], mapped_coord[DIM];
        IndexType mapped_index;
        int mapped_rel_grid_loc[2];
        IndexType count = 0;
        while(h < next_to_last_halo){
            this->lattice_->local_mem_index2coord(h, halo_coord);
            this->lattice_->get_mapped_coord(halo_coord, mapped_coord);
            mapped_index = this->lattice_->local_mem_coord2index(mapped_coord);
            this->lattice_->get_mapped_relative_grid_loc(halo_coord, mapped_rel_grid_loc);
            update_halo_table_[count] = 
                    std::make_tuple(h, mapped_index, mapped_rel_grid_loc[0], mapped_rel_grid_loc[1]);
            h = this->lattice_->get_local_halo_next(h);
            count++;
        }
    }


    template<class FieldType, int DIM>
    void Field<FieldType, DIM>::update_halo() const {

        for(const auto& t : this->update_halo_table_){
            auto recv_start = this->data_ptr_ + std::get<0>(t)*components_;
            auto send_start = this->data_ptr_ + std::get<1>(t)*components_;

            if(std::get<2>(t) == 0 && std::get<3>(t) == 0){
                //self sync
                std::copy_n(send_start, components_, recv_start);
            } else if(std::get<2>(t) == 1 && std::get<3>(t) == 0){
                //y axis send-down
                this->parallel_ptr->SendDown_col(send_start, recv_start, components_);
            } else if(std::get<2>(t) == -1 && std::get<3>(t) == 0){
                //y axis send-up
                this->parallel_ptr->SendUp_col(send_start, recv_start, components_);
            } else if(std::get<2>(t) == 0 && std::get<3>(t) == 1){
                //z axis send-down
                this->parallel_ptr->SendDown_row(send_start, recv_start, components_);
            } else if(std::get<2>(t) == 0 && std::get<3>(t) == -1){
                //z axis send-up
                this->parallel_ptr->SendUp_row(send_start, recv_start, components_);
            } else{
                GridIndexType recv_from_grid_loc[2], recv_from_grid_rank[2];
                GridIndexType send_to_grid_loc[2], send_to_grid_rank[2];
                GridIndexType rel_grid_loc[] = {std::get<2>(t), std::get<3>(t)};
                GridIndexType rel_opp_grid_loc[] = {-std::get<2>(t), -std::get<3>(t)};
                this->lattice_->get_mapped_grid_loc_from_rel_grid_loc(rel_grid_loc, recv_from_grid_loc);
                this->lattice_->get_mapped_grid_loc_from_rel_grid_loc(rel_opp_grid_loc, send_to_grid_loc);
                transform_gridloc_to_gridrank(recv_from_grid_loc, recv_from_grid_rank);
                transform_gridloc_to_gridrank(send_to_grid_loc, send_to_grid_rank);
                this->parallel_ptr->Send(send_start, components_, send_to_grid_rank);
                this->parallel_ptr->Receive(recv_start, components_, recv_from_grid_rank);
            }
        }
    }

    template<class FieldType, int DIM>
    void Field<FieldType, DIM>::make_fileio_settings(
            hsize_t* dataset_size,
            hsize_t* file_offset,
            hsize_t* file_block_size,
            hsize_t* mem_size,
            hsize_t* mem_offset,
            hsize_t* mem_block_size) const {
        // dataset size
        auto ptr_glb_size = this->lattice_->get_global_size();
        std::vector<IndexType> lat_size(ptr_glb_size, ptr_glb_size + DIM);
        // The save order is reversed: (z, y, x, components) 
        std::reverse_copy(lat_size.begin(), lat_size.end(), dataset_size);
        dataset_size[DIM] = this->components_;

        // mem_size
        auto ptr_mem_size = this->lattice_->get_local_mem_size();
        std::vector<IndexType> vec_mem_size(ptr_mem_size, ptr_mem_size + DIM);
        std::reverse_copy(vec_mem_size.begin(), vec_mem_size.end(), mem_size);
        mem_size[DIM] = this->components_;

        //file_offset
        IndexType tmp_lc[DIM], tmp_gc[DIM];
        std::fill_n(tmp_lc, DIM, 0);
        this->lattice_->local_vis_coord_to_global_coord(tmp_lc, tmp_gc);
        std::reverse_copy(tmp_gc, tmp_gc+DIM, file_offset);
        file_offset[DIM] = 0;

        //mem_offset
        std::fill_n(mem_offset, DIM, this->lattice_->get_halo());
        mem_offset[DIM] = 0;

        //file_block_size (also mem_block_size)
        auto ptr_vis_size = this->lattice_->get_local_size();
        std::vector<IndexType> vis_size(ptr_vis_size, ptr_vis_size + DIM);
        std::reverse_copy(vis_size.begin(), vis_size.end(), file_block_size);
        file_block_size[DIM] = this->components_;

        std::copy_n(file_block_size, DIM+1, mem_block_size);

        return;
    }

    template<class FieldType, int DIM>
    void Field<FieldType, DIM>::write(
        const std::string& file_name,
        const std::string& dataset_name) const{
        FileProcessor fp;

        std::string dst_name;
        if(dataset_name == "") dst_name = file_name;
        else dst_name = dataset_name;

        hsize_t dataset_size[DIM+1];
        hsize_t file_offset[DIM+1];
        hsize_t file_block_size[DIM+1];
        hsize_t mem_size[DIM+1];
        hsize_t mem_offset[DIM+1];
        hsize_t mem_block_size[DIM+1];

        this->make_fileio_settings(
            dataset_size,
            file_offset,
            file_block_size,
            mem_size,
            mem_offset,
            mem_block_size);
        
        int msg_send = 1;
        int msg_recv = 0;
        const auto& my_rank = this->parallel_ptr->get_world_rank();
        const auto& world_size = this->parallel_ptr->get_world_size();
        if(my_rank == 0){
            fp.SaveSingleDatasetFile<FieldType, DIM+1>(
                                file_name,
                                dst_name,
                                dataset_size,
                                file_offset,
                                file_block_size,
                                mem_size,
                                mem_offset,
                                mem_block_size, 
                                this->data_ptr_,
                                true);
            if(world_size > 1)
                this->parallel_ptr->Send(msg_send, 1);
        } else{
            this->parallel_ptr->Receive(msg_recv, my_rank - 1);
            fp.SaveSingleDatasetFile<FieldType, DIM+1>(
                                file_name,
                                dst_name,
                                dataset_size,
                                file_offset,
                                file_block_size,
                                mem_size,
                                mem_offset,
                                mem_block_size, 
                                this->data_ptr_,
                                false);
            if(my_rank < world_size - 1)
                this->parallel_ptr->Send(msg_send, my_rank + 1);
        }
        return;
    }

    template<class FieldType, int DIM>
    void Field<FieldType, DIM>::read(
        const std::string& file_name,
        const std::string& dataset_name) const{
        FileProcessor fp;

        std::string dst_name;
        if(dataset_name == "") dst_name = file_name;
        else dst_name = dataset_name;

        hsize_t dataset_size[DIM+1];
        hsize_t file_offset[DIM+1];
        hsize_t file_block_size[DIM+1];
        hsize_t mem_size[DIM+1];
        hsize_t mem_offset[DIM+1];
        hsize_t mem_block_size[DIM+1];

        this->make_fileio_settings(
            dataset_size,
            file_offset,
            file_block_size,
            mem_size,
            mem_offset,
            mem_block_size);

        int msg_send = 1;
        int msg_recv = 0;
        const auto& my_rank = this->parallel_ptr->get_world_rank();
        const auto& world_size = this->parallel_ptr->get_world_size();
        if(my_rank == 0){
            fp.ReadSingleDatasetFile<FieldType, DIM+1>(
                                file_name,
                                dst_name,
                                dataset_size,
                                file_offset,
                                file_block_size,
                                mem_size,
                                mem_offset,
                                mem_block_size, 
                                this->data_ptr_,
                                true);
            if(world_size > 1)
                this->parallel_ptr->Send(msg_send, 1);
        } else{
            this->parallel_ptr->Receive(msg_recv, my_rank - 1);
            fp.ReadSingleDatasetFile<FieldType, DIM+1>(
                                file_name,
                                dst_name,
                                dataset_size,
                                file_offset,
                                file_block_size,
                                mem_size,
                                mem_offset,
                                mem_block_size, 
                                this->data_ptr_,
                                false);
            if(my_rank < world_size - 1)
                this->parallel_ptr->Send(msg_send, my_rank + 1);
        }
        this->update_halo();
        return;
    }

}


#endif