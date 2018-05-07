#ifndef LATFIELD3_FIELD_H
#define LATFIELD3_FIELD_H

#include "LATfield3_Site.h"
#include <string>
#include <fstream>
#include <iostream>

namespace LATfield3{

    template<class FieldType>
    class Field{
    public:
        Field() = delete;
        /*
        define a vector field, with dimension = components.
        */
        Field(const Lattice& lattice, const int components /*mat_rows*/ = 1);
        Field(const Lattice& lattice, const int mat_rows = 1, const int mat_cols = 1);
        ~Field();

        //allocator
        allocator();
        deallocator();

        //data retrieval and modification
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
        void read(const std::string& file_name);
        void write(const std::string& file_name) const;

        void update_halo() const;
        template<int DIM>
        const Lattice& get_lattice() const {return *(this->lattice_);}
        
        Lattice& get_lattice() const {return const_cast<Lattice>(*(this->lattice_));}
        int get_rows() const {return this->mat_rows_;}
        int get_cols() const {return this->mat_cols_;}
    private:
        FieldType* data_ = nullptr;

        const Lattice* lattice_;
        /*
        stores the shape of the matrix field.
        if it is a vector field, then the components will be stored in mat_rows_,
        while mat_cols will be set to 1.
        */
        int mat_rows_;
        int mat_cols_;
    };
}


#endif