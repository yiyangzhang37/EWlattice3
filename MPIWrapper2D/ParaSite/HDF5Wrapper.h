#ifndef HDF5_WRAPPER_H
#define HDF5_WRAPPER_H

#include "MPIWrapper2D.h"
#include <hdf5_hl.h>
#include <iostream>
#include <algorithm>

namespace HDF5_Wrapper{

    typedef MPI_Wrapper::Parallel2D ParallelObject;

    template<class Type>
    hid_t get_H5_datatype(){
        return -1;
    }

    class HDF5Wrapper{
    private:
        //const ParallelObject* parallel_ = nullptr;
        std::string file_name_;
        hid_t file_id_ = -1;
        hid_t dataset_id_ = -1;

        hid_t current_loc_ = -1; /*file or group identifier*/

    public:
        HDF5Wrapper() = default;

        /*
        when the constructor is called, the file "file_name" is opened.
        If the file does not exist (or the file is not h5 format), it will create one.
        one object will only open one HDF5 file.
        */
        HDF5Wrapper(const std::string& file_name);
        //HDF5Wrapper(const ParallelObject& parallel);
        
        /*
        when the destructor is called, the h5 file is closed.
        */
        ~HDF5Wrapper();

        /*************** High-level functions ***************/
        /*
        Save a single dataset file.
        The file will be created in TRUNC mode. All previous info will be erased.
        */
        template<class Type, int DIM>
        void SaveSingleDatasetFile(
            const std::string& file_name,
            const std::string& dataset_name,
            const hsize_t* dataset_size,
            const Type* data);
        /*
        Save a single dataset file.
        Each call of this function will add data to partial of the dataset.
        The first call (specified by is_first == true) will create a new file in TRUNC mode,
        the following calls will simply open that file.
        */
        
        template<class Type, int DIM>
        void SaveSingleDatasetFile(
            const std::string& file_name,
            const std::string& dataset_name,
            const hsize_t* dataset_size,
            const hsize_t* file_offset,
            const hsize_t* file_block_size,
            const hsize_t* mem_size,
            const hsize_t* mem_offset,
            const hsize_t* mem_block_size,
            const Type* data, 
            const bool is_first);
        
        /*
        Partially read data from one file by each call.
        */
        template<class Type, int DIM>
        void ReadSingleDatasetFile(
            const std::string& file_name,
            const std::string& dataset_name,
            const hsize_t* dataset_size,
            Type* data);

        template<class Type, int DIM>
        void ReadSingleDatasetFile(
            const std::string& file_name,
            const std::string& dataset_name,
            const hsize_t* dataset_size,
            const hsize_t* file_offset,
            const hsize_t* file_block_size,
            const hsize_t* mem_size,
            const hsize_t* mem_offset,
            const hsize_t* mem_block_size,
            Type* data);

        /*************** File objects ***************/
        /*
        create a new h5 file, and move the current_loc_ to this->file_id_.
        accepted flags: 
        H5F_ACC_TRUNC: truncate file.
        H5F_ACC_EXCL: fail if it already exists.
        If creation fails, return a negative value.
        */
        hid_t create_file(
            const std::string& file_name, 
            const unsigned flag = H5F_ACC_EXCL);

        /*
        open a h5 file, and move the current_loc_ to this->file_id_.
        accepted flags:
        H5F_ACC_RDONLY: read_only.
        H5F_ACC_RDWR: read and write.
        If open fails, return a negative value.
        */
        hid_t open_file(
            const std::string& file_name, 
            const unsigned flag = H5F_ACC_RDWR);
        /*
        check whether a h5 file exists.
        If the file exists and is HDF5 format, open the file.
        If it does not exist, create the file with truncation mode.
        */
        hid_t create_or_open_file(
            const std::string& file_name,
            const unsigned create_flag = H5F_ACC_TRUNC,
            const unsigned open_flag = H5F_ACC_RDWR);
        
        herr_t close_file();

        /*************** Group objects ***************/

        /*************** Dataset objects ***************/
        /*
        create dataset in this->current_loc_.
        */
        template<class Type, int DIM>
        hid_t create_dataset(
           const std::string& dataset_name,
           const hsize_t* size);
        
        hid_t open_dataset(
            const std::string& dataset_name);

        /*
        In this->current_loc_,
        create a dataset if it does not exist, otherwise open the existing dataset.
        */
        template<class Type, int DIM>
        hid_t create_or_open_dataset(
            const std::string& dataset_name,
            const hsize_t* size);
        
        template<int DIM>
        herr_t reset_dataset_size(
            const hid_t dataset_id,
            const hsize_t* size) const;

        herr_t close_dataset();
        
        /*
        at least by default, the array indexing convention is Last-index-run-first.
        */
        /* write/read the whole data buffer to/from the whole dataset */
        template<class Type>
        herr_t write_dataset(const Type* data) const; 
        template<class Type>
        herr_t read_dataset(Type* data) const;

        /* write/read the whole data buffer to selected region in the dataset*/
        template<class Type, int DIM>
        herr_t write_dataset(
            const Type* data, 
            const hsize_t* file_offset,
            const hsize_t* file_block_size,
            const hsize_t* mem_size,
            const hsize_t* mem_offset,
            const hsize_t* mem_block_size) const;
        template<class Type, int DIM>
        herr_t read_dataset(
            Type* data, 
            const hsize_t* file_offset,
            const hsize_t* file_block_size,
            const hsize_t* mem_size,
            const hsize_t* mem_offset,
            const hsize_t* mem_block_size) const;


        /*************** Dataspace objects ***************/
        /*
        create a simple dataspace.
        */
        template<int DIM>
        hid_t create_dataspace(const hsize_t* size) const;

        /*
        select one subregion given the dataspace_id.
        The stride and block count is 1 in this function.
        */
        template<int DIM>
        herr_t select_subregion(
            const hid_t dataspace_id,
            const hsize_t* offset, 
            const hsize_t* block_size) const;

        /*
        get the dataspace object associated with this->dataset_id_,
        The this->dataspace_id_ will be refreshed by this function.
        H5Sclose() needs to be called once this object is no longer needed.
        */
        hid_t get_dataspace() const;

        herr_t close_dataspace(const hid_t dataspace_id) const;

        /*************** Attributes objects ***************/
        template<class Type>
        herr_t add_attribute(
            const hid_t loc_id,
            const std::string& name,
            const Type val) const;
        
        template<class Type>
        herr_t add_attribute(
            const hid_t loc_id,
            const std::string& name,
            const Type* vals,
            const hsize_t len) const;
        
        template<class Type>
        herr_t attach_attribute_to_dataset(
            const std::string& name,
            const Type val) const;
        
        template<class Type>
        herr_t attach_attribute_to_dataset(
            const std::string& name,
            const Type* vals,
            const hsize_t len) const;

        herr_t delete_attribute(
            const hid_t loc_id,
            const std::string& name) const;

        /*Miscellaneous*/
        std::string get_file_name() const; 
        hsize_t get_file_size() const;
        hid_t get_file_id() const {return this->file_id_;}
        hid_t get_dataset_id() const {return this->dataset_id_;}
        hid_t get_current_loc() const {return this->current_loc_;}
    

    };

    template<class Type, int DIM>
    hid_t HDF5Wrapper::create_dataset(
        const std::string& dataset_name,
        const hsize_t* size){
        auto dataspace_id = this->create_dataspace<DIM>(size);
        auto status = H5Dcreate(
            this->current_loc_, dataset_name.c_str(),
            get_H5_datatype<Type>(),
            dataspace_id,
            H5P_DEFAULT /*Link creation property list*/,
            H5P_DEFAULT /*Dataset creation property list*/,
            H5P_DEFAULT /*Dataset access property list*/);
        if(status < 0) std::cerr << "Dataset creation error." << std::endl;
        this->dataset_id_ = status;
        this->close_dataspace(dataspace_id);
        return this->dataset_id_;
    }

    template<class Type, int DIM>
    hid_t HDF5Wrapper::create_or_open_dataset(
        const std::string& dataset_name,
        const hsize_t* size) {
        //check if a dataset exists //todo
        hid_t status = -1;
        if(1){
            status = this->create_dataset<Type, DIM>(dataset_name, size);
        }
        else{
            status = this->open_dataset(dataset_name);
            status = this->reset_dataset_size<DIM>(this->dataset_id_, size);
        }
        return status;
    }

    template<int DIM>
    herr_t HDF5Wrapper::reset_dataset_size(
        const hid_t dataset_id,
        const hsize_t* size) const {
        return H5Dset_extent(dataset_id, size);
    }

    template<int DIM>
    hid_t HDF5Wrapper::create_dataspace(const hsize_t* size) const {
        return H5Screate_simple(DIM, size, nullptr);
    }


    template<class Type>
    herr_t HDF5Wrapper::write_dataset(const Type* data) const {
        return H5Dwrite(this->dataset_id_, 
                        get_H5_datatype<Type>(), 
                        H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    }

    template<class Type>
    herr_t HDF5Wrapper::read_dataset(Type* data) const {
        return H5Dread(this->dataset_id_, 
                        get_H5_datatype<Type>(), 
                        H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    }

    template<class Type, int DIM>
    herr_t HDF5Wrapper::write_dataset(
        const Type* data, 
        const hsize_t* file_offset,
        const hsize_t* file_block_size,
        const hsize_t* mem_size,
        const hsize_t* mem_offset,
        const hsize_t* mem_block_size) const {
        auto file_dsid = this->get_dataspace();
        auto status = this->select_subregion<DIM>(file_dsid, file_offset, file_block_size);
        auto mem_dsid = H5Screate_simple(DIM, mem_size, nullptr);
        status = this->select_subregion<DIM>(mem_dsid, mem_offset, mem_block_size);
        if(status < 0) std::cerr << "Select subregion error." << std::endl;
        status = H5Dwrite(this->dataset_id_,
                        get_H5_datatype<Type>(),
                        mem_dsid,
                        file_dsid,
                        H5P_DEFAULT,
                        data);
        this->close_dataspace(file_dsid);
        this->close_dataspace(mem_dsid);
        return status;
    }

    template<class Type, int DIM>
    herr_t HDF5Wrapper::read_dataset(
                    Type* data, 
                    const hsize_t* file_offset,
                    const hsize_t* file_block_size,
                    const hsize_t* mem_size,
                    const hsize_t* mem_offset,
                    const hsize_t* mem_block_size) const {
        auto file_dsid = this->get_dataspace();
        auto status = this->select_subregion<DIM>(file_dsid, file_offset, file_block_size);
        auto mem_dsid = H5Screate_simple(DIM, mem_size, nullptr);
        status = this->select_subregion<DIM>(mem_dsid, mem_offset, mem_block_size);
        if(status < 0) std::cerr << "Select subregion error." << std::endl;
        status = H5Dread(this->dataset_id_,
                        get_H5_datatype<Type>(),
                        mem_dsid,
                        file_dsid,
                        H5P_DEFAULT,
                        data);
        this->close_dataspace(file_dsid);
        this->close_dataspace(mem_dsid);
        return status;
    }


    template<int DIM>
    herr_t HDF5Wrapper::select_subregion(
        const hid_t dataspace_id,
        const hsize_t* offset, 
        const hsize_t* block_size) const{
        hsize_t stride[DIM];
        hsize_t count[DIM];
        std::fill_n(stride, DIM, 1);
        std::fill_n(count, DIM, 1);
        return H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET,
                                    offset, 
                                    stride,
                                    count,
                                    block_size);
    }

    template<class Type>
    herr_t HDF5Wrapper::add_attribute(
            const hid_t loc_id,
            const std::string& name,
            const Type val) const {
        return this->add_attribute(loc_id, name, &val, 1);
    }

    template<class Type>
    herr_t HDF5Wrapper::add_attribute(
            const hid_t loc_id,
            const std::string& name,
            const Type* vals,
            const hsize_t len) const {
        auto is_exist = H5Aexists(loc_id, name.c_str());
        if(is_exist > 0 ){
            std::cerr << "Attribute already exists." << std::endl;
            return -1;
        }
        hsize_t dims[] = {len};
        auto attr_dataspace_id = H5Screate_simple(1, dims, nullptr);
        auto attr_id =  H5Acreate(loc_id, 
                            name.c_str(), 
                            get_H5_datatype<Type>(), 
                            attr_dataspace_id,
                            H5P_DEFAULT,
                            H5P_DEFAULT);
        H5Awrite(attr_id, get_H5_datatype<Type>(), vals);
        H5Sclose(attr_dataspace_id);
        return H5Aclose(attr_id);
    }

    template<class Type>
    herr_t HDF5Wrapper::attach_attribute_to_dataset(
        const std::string& name,
        const Type val) const {
        return this->add_attribute(this->dataset_id_, name, val);
    }

    template<class Type>
    herr_t HDF5Wrapper::attach_attribute_to_dataset(
        const std::string& name,
        const Type* vals,
        const hsize_t len) const {
        return this->add_attribute(this->dataset_id_, name, vals, len);
    }

    template<class Type, int DIM>
    void HDF5Wrapper::SaveSingleDatasetFile(
        const std::string& file_name,
        const std::string& dataset_name,
        const hsize_t* size,
        const Type* data) {
        //If there is any file opened, close it first.
        if(this->file_id_ >= 0) this->close_file();

        this->create_file(file_name, H5F_ACC_TRUNC);
        this->create_dataset<Type, DIM>(dataset_name, size);
        this->write_dataset<Type>(data);
        this->close_dataset();
        this->close_file();
        return;
    }

    template<class Type, int DIM>
    void HDF5Wrapper::SaveSingleDatasetFile(
            const std::string& file_name,
            const std::string& dataset_name,
            const hsize_t* dataset_size,
            const hsize_t* file_offset,
            const hsize_t* file_block_size,
            const hsize_t* mem_size,
            const hsize_t* mem_offset,
            const hsize_t* mem_block_size,
            const Type* data, 
            const bool is_first) {
        //If there is any file opened, close it first.
        if(this->file_id_ >= 0) this->close_file();

        if(is_first){
            this->create_file(file_name, H5F_ACC_TRUNC);
            this->create_dataset<Type,DIM>(dataset_name, dataset_size);
        } else{
            this->open_file(file_name, H5F_ACC_RDWR);
            this->open_dataset(dataset_name);
        }

        this->write_dataset<Type, DIM>(data, 
                                    file_offset, 
                                    file_block_size,
                                    mem_size, 
                                    mem_offset, 
                                    mem_block_size);

        this->close_dataset();
        this->close_file();
        return;
    }

    template<class Type, int DIM>
    void HDF5Wrapper::ReadSingleDatasetFile(
        const std::string& file_name,
        const std::string& dataset_name,
        const hsize_t* dataset_size,
        const hsize_t* file_offset,
        const hsize_t* file_block_size,
        const hsize_t* mem_size,
        const hsize_t* mem_offset,
        const hsize_t* mem_block_size,
        Type* data){
        //If there is any file opened, close it first.
        if(this->file_id_ >= 0) this->close_file();

        this->open_file(file_name, H5F_ACC_RDONLY);
        this->open_dataset(dataset_name);
        this->read_dataset<Type, DIM>(data, 
                                    file_offset, 
                                    file_block_size,
                                    mem_size,
                                    mem_offset,
                                    mem_block_size);

        this->close_dataset();
        this->close_file();
        return;
    }

    template<class Type, int DIM>
    void HDF5Wrapper::ReadSingleDatasetFile(
        const std::string& file_name,
        const std::string& dataset_name,
        const hsize_t* size,
        Type* data){
        hsize_t offset[DIM];
        hsize_t block_size[DIM];
        std::fill_n(offset, DIM, 0);
        std::copy_n(size, DIM, block_size);
        this->ReadSingleDatasetFile<Type, DIM>(file_name, 
                                                dataset_name, 
                                                size, 
                                                offset,
                                                block_size,
                                                block_size,
                                                offset,
                                                block_size,
                                                data);
        return;
    }

    

    








}


#endif