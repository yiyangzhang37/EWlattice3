#include "./HDF5Wrapper.h"

namespace HDF5_Wrapper{

    template<>
    hid_t get_H5_datatype<int>(){
        return H5T_NATIVE_INT;
    }

    template<>
    hid_t get_H5_datatype<unsigned int>(){
        return H5T_NATIVE_UINT;
    }

    template<>
    hid_t get_H5_datatype<float>(){
        return H5T_NATIVE_FLOAT;
    }
    
    template<>
    hid_t get_H5_datatype<double>(){
        return H5T_NATIVE_DOUBLE;
    }

    template<>
    hid_t get_H5_datatype<char*>(){
        return H5T_C_S1;
    }

    HDF5Wrapper::HDF5Wrapper(const std::string& file_name){
        this->create_or_open_file(file_name);
    }

    HDF5Wrapper::~HDF5Wrapper(){
        //close_file() cannot be put here.
        //The destructor may free the file_id_ first.
    }
    
    hid_t HDF5Wrapper::create_file(
        const std::string& file_name, 
        const unsigned flag) {
        this->file_name_ = file_name;
        this->file_id_ = H5Fcreate(file_name.c_str(), flag, H5P_DEFAULT, H5P_DEFAULT);
        this->current_loc_ = this->file_id_;
        return this->file_id_;
    }

    hid_t HDF5Wrapper::open_file(
        const std::string& file_name, 
        const unsigned flag){
        this->file_name_ = file_name;
        this->file_id_ = H5Fopen(file_name.c_str(), flag, H5P_DEFAULT);
        this->current_loc_ = this->file_id_;
        return this->file_id_;
    }

    hid_t HDF5Wrapper::create_or_open_file(
        const std::string& file_name,
        const unsigned create_flag,
        const unsigned open_flag){
        //test whether the file exists.
        auto is_hdf5 = H5Fis_hdf5(file_name.c_str());
        auto file_id = this->file_id_;
        if(is_hdf5 < 0){
            file_id = this->create_file(file_name, create_flag);
        }
        else{
            file_id = this->open_file(file_name, open_flag);
        }
        if(file_id < 0) std::cerr << "create or open file error." << std::endl;
        else this->file_id_ = file_id;
        return this->file_id_;
    }

    herr_t HDF5Wrapper::close_file() {
        if(this->file_id_ >= 0){
            auto status = H5Fclose(this->file_id_);
            if(status >= 0){
                this->file_id_ = -1;    
            } else{
                std::cerr << "File close error." << std::endl;
            }
            return status;
        }
        else {
            std::cerr << "No file opened" << std::endl;
            return -1;
        }
    }

    herr_t HDF5Wrapper::close_dataspace(const hid_t dataspace_id) const{
        return H5Sclose(dataspace_id);     
    }
    

    hid_t HDF5Wrapper::open_dataset(
        const std::string& dataset_name){
        auto dataset_id = H5Dopen(this->file_id_, dataset_name.c_str(), H5P_DEFAULT);
        if(dataset_id < 0) std::cerr << "Dataset open failed." << std::endl;
        this->dataset_id_ = dataset_id;
        return this->dataset_id_;
    }

    herr_t HDF5Wrapper::close_dataset(){
        if(this->dataset_id_ >= 0){
            auto status = H5Dclose(this->dataset_id_);
            if(status >= 0)
                this->dataset_id_ = -1;
            else
                std::cerr << "Dataset close error." << std::endl;
            return status;
        }
        else{
            std::cerr << "No datasets opened." << std::endl;
            return -1;
        }
    }

    std::string HDF5Wrapper::get_file_name() const{
        char name[256];
        auto len = H5Fget_name(this->file_id_, name, 256);
        if(len < 0) std::cerr << "get file name error." << std::endl;
        return std::string(name, 0, len);

    }

    hsize_t HDF5Wrapper::get_file_size() const{
        hsize_t size;
        H5Fget_filesize(this->file_id_, &size);
        return size;
    }

    hid_t HDF5Wrapper::get_dataspace() const{
        return H5Dget_space(this->dataset_id_);
    }

    herr_t HDF5Wrapper::delete_attribute(
        const hid_t loc_id,
        const std::string& name) const {
        auto is_exist = H5Aexists(loc_id, name.c_str());
        if(is_exist <= 0){
            std::cerr << "Attribute not exist." << std::endl;
            return -1;
        } else{
            return H5Adelete(loc_id, name.c_str());
        }
    }

    template<>
    herr_t HDF5Wrapper::add_attribute<std::string> (
            const hid_t loc_id,
            const std::string& name,
            const std::string val) const{
        auto is_exist = H5Aexists(loc_id, name.c_str());
        if(is_exist > 0 ){
            std::cerr << "Attribute already exists." << std::endl;
            return -1;
        }
        hsize_t dims[] = {1};
        auto str_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(str_type, H5T_VARIABLE);
        auto attr_dsid = H5Screate_simple(1, dims, nullptr);
        auto attr_id =  H5Acreate(loc_id, 
                            name.c_str(), 
                            str_type, 
                            attr_dsid,
                            H5P_DEFAULT,
                            H5P_DEFAULT);
        char* cstr = new char [val.length()+1];
        std::strcpy(cstr, val.c_str());
        H5Awrite(attr_id, str_type, &cstr);
        H5Sclose(attr_dsid);
        H5Tclose(str_type);
        return H5Aclose(attr_id);
    }






}

