#ifndef LATFIELD3_FIELD_H
#define LATFIELD3_FIELD_H

#include "LATfield3_Site.h"
#include <string>
#include <fstream>
#include <iostream>
#include <memory>

namespace LATfield3{

    template<class FieldType, int DIM>
    class Field{
    public:
        Field() = delete;
        /*
        define a vector field, with dimension = components.
        */
        Field(const Lattice<DIM>& lattice, const int mat_rows = 1);
        Field(const Lattice<DIM>& lattice, const int mat_rows = 1, const int mat_cols = 1);
        ~Field() = default;

        //data retrieval and modification
        constexpr int get_component(const int i_row, const int j_col) const{
            return i_row * stride_[0] + j_col * stride_[1];
        }

        FieldType& operator()(const Site& site) const;
        FieldType& operator()(const Site& site, const int i_row) const;
        FieldType& operator()(const Site& site, const int i_row, const int j_col) const;
        const FieldType& operator()(const Site& site) const;
        const FieldType& operator()(const Site& site, const int i_row) const;
        const FieldType& operator()(const Site& site, const int i_row, const int j_col) const;

        FieldType& operator()(const IndexType index) const;
        FieldType& operator()(const IndexType index, const int i_row) const;
        FieldType& operator()(const IndexType index, const int i_row, const int j_col) const;
        const FieldType& operator()(const IndexType index) const;
        const FieldType& operator()(const IndexType index, const int i_row) const;
        const FieldType& operator()(const IndexType index, const int i_row, const int j_col) const;
        
        //file operations
        void read(const std::string& file_name) const;
        void write(const std::string& file_name) const;

        void update_halo() const;
        
        const Lattice& get_lattice() const {return *(this->lattice_);}
        Lattice& get_lattice() const {return const_cast<Lattice>(*(this->lattice_));}
        constexpr int get_rows() const {return this->mat_rows_;}
        constexpr int get_cols() const {return this->mat_cols_;}
    private:
        //allocator
        allocator();

        std::unique_ptr<FieldType> data_ptr_ = nullptr;

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
    };

    template<class FieldType, int DIM>
    Field<FieldType, DIM>::Field(const Lattice<DIM>& lattice, const int mat_rows){
        Field(lattice, mat_rows, 1);
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
        assert(mat_cols_ > 0)
        components_ = mat_rows_ * mat_cols_;
        stride_[0] = 1;
        stride_[1] = mat_rows;
        allocator();
    }

    template<class FieldType, int DIM>
    void Field<FieldType, DIM>::allocator(){
        data_ptr_ = std::unique_ptr<FieldType>(
            new FieldType[this->lattice_->get_local_mem_sites() * mat_rows * mat_cols] );
    }

    template<class FieldType, int DIM>
    FieldType& Field<FieldType, DIM>::operator()(const Site& site) const {
        return data_ptr_[site.get_index()*this->components_];
    }

    template<class FieldType, int DIM>
    const FieldType& Field<FieldType, DIM>::operator()(const Site& site) const {
        return data_ptr_[site.get_index()*this->components_];
    }

    template<class FieldType, int DIM>
    FieldType& Field<FieldType, DIM>::operator()(const Site& site, const int i_row) const {
        return data_ptr_[site.get_index()*this->components_ + i_row];
    }

    template<class FieldType, int DIM>
    const FieldType& Field<FieldType, DIM>::operator()(const Site& site, const int i_row) const {
        return data_ptr_[site.get_index()*this->components_ + i_row];
    }

    template<class FieldType, int DIM>
    FieldType& 
    Field<FieldType, DIM>::operator()(const Site& site, const int i_row, const int j_col) const {
        return data_ptr_[site.get_index()*this->components_ + get_component(i_row, j_col)];
    }

    template<class FieldType, int DIM>
    const FieldType& 
    Field<FieldType, DIM>::operator()(const Site& site, const int i_row, const int j_col) const {
        return data_ptr_[site.get_index()*this->components_ + get_component(i_row, j_col)];
    }

    template<class FieldType, int DIM>
    void Field<FieldType, DIM>::update_halo() const {
        //TODO
    }





}


#endif