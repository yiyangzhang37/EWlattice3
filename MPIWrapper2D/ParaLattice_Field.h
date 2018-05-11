#ifndef PARA_LATTICE_FIELD_H
#define PARA_LATTICE_FIELD_H

#include "ParaLattice_Site.h"
#include "MPIWrapper2D.h"
#include <string>
#include <fstream>
#include <iostream>
#include <memory>
#include <tuple>
#include <type_traits>


namespace ParaLattice{

    using ParallelObject = MPI_Wrapper::Parallel2D;

    template<class FieldType, int DIM>
    class Field{
    private:
        //allocator
        void allocator();
        
        //make update_halo_table_
        void make_update_halo_table();
    
        //std::unique_ptr<FieldType> data_ptr_ = nullptr;
        FieldType* data_ptr_;

        const Lattice<DIM>* lattice_;
        /*
        stores the shape of the matrix field.
        if it is a vector field, then the components will be stored in mat_rows_,
        while mat_cols will be set to 1.
        / It will follow row-first convention.
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
        ~Field();

        void assign_parallel_object(const ParallelObject& parallel_object);

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
        void read(const std::string& file_name) const;
        void write(const std::string& file_name) const;

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
    Field<FieldType, DIM>::~Field(){
        delete[] this->data_ptr_;
    }

    template<class FieldType, int DIM>
    void Field<FieldType, DIM>::allocator(){
        //data_ptr_ = std::unique_ptr<FieldType>(
        //    new FieldType[this->lattice_->get_local_mem_sites() * components_] );
        data_ptr_ = new FieldType[this->lattice_->get_local_mem_sites() * components_];
    }

    template<class FieldType, int DIM>
    void Field<FieldType, DIM>::assign_parallel_object(const ParallelObject& parallel_object){
        this->parallel_ptr = &parallel_object;
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
            IndexType mapped_index = this->lattice_->local_mem_coord2index(mapped_coord);
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
                this->lattice_->get_mapped_grid_loc(rel_grid_loc, recv_from_grid_loc);
                this->lattice_->get_mapped_grid_loc(rel_opp_grid_loc, send_to_grid_loc);
                transform_gridloc_to_gridrank(recv_from_grid_loc, recv_from_grid_rank);
                transform_gridloc_to_gridrank(send_to_grid_loc, send_to_grid_rank);
                this->parallel_ptr->Send(send_start, components_, send_to_grid_rank);
                this->parallel_ptr->Receive(recv_start, components_, recv_from_grid_rank);
            }
        }
    }


}


#endif