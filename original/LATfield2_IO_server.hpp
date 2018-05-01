#ifndef LATFIELD2_IOSERVER_HPP
#define LATFIELD2_IOSERVER_HPP

/*! \file LATfield2_IO_server.hpp
 \brief LATfield2_IO_server.hpp contains the class IOserver definition.
 \author David Daverio
 */ 

#include "mpi.h"
#include "int2string.hpp"

#define SERVER_STATE_TAG 1

#define SERVER_CONTROL_TAG 2

#define CONTROL_OPEN_OSTREAM 1
#define CONTROL_STOP 2

#define IO_FILE_CONTROL_TAG 3

#define CONTROL_CREATE_FILE 1
#define CONTROL_CLOSE_OSTREAM 2

#define IO_FILE_CONTROL_FILEID_TAG 4
#define IO_FILE_CONTROL_FILENAME_TAG 5


#define OSTREAM_SUCCESS 1
#define OSTREAM_FAIL 0
#define FILE_FAIL -32000
#define MAX_FILE_NUMBER 5
#define IO_FILE_CONTROL_CLOSE_TAG 15

#define FILETYPE_UNSTRUCTURED 1

//byte given;
#define IO_BUFFERS_TOTAL_SIZE 873741824 

typedef int ioserver_file;

//! A structure to describe a file for the I/O server (dedicated MPI processes for writing to disks)
struct file_struct{
    //! path to the file
    string filename;
    //! data array of the file
    char * data;
    //! size of the local part of the file
    long long size; 
    //! type of the file (currently only FILETYPE_UNSTRUCTURED)
    int type;
};



/*! \class IOserver  
 
 \brief A class to handle the I/O using MPI process reserved for IO purpose on which the files are defined
 
 This server is in beta stage, but as such a functionality is very useful, it has been added to the stable part of LATfield2. An example of the usage of this class is given in the IOserver example. User should never instanciate an IOserver object. The IOserver objet (IO_Server) is instanciate within the library header
 
 
 */
class IOserver {
    
public:
      
    ~IOserver();
    
    
  
    
    /*! \brief Server method (only called by server nodes)
     
     Method which is called to start the server.
     */
    void start(); 

    
    

    
    /*! \brief Client method (only called by compute nodes)
     
     Method which is called to stop the server.
     
     */
    void stop();  //client
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to open an Ostream. Meaning a stream from the compute to the server processes.
     
     \return OSTREAM_SUCCESS if the stream is open.
     \return OSTREAM_FAIL  if the stream cannot be open.
     */
    int openOstream();
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to close the current Ostream. After the stream is closed, the server will start to write the files it have in memory.
     
    */
    void closeOstream();
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to create a new file, it return the fileID.
     
     \param  filename: name of the file (including the path...)
     \return  fileID.
     */
    ioserver_file createFile(string filename);
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to close a new file: fileID.
     
     \param ioserver_file fileID: file to close.
     */
    void closeFile(ioserver_file fileID);
    
    /*! \brief Client method (only called by compute nodes)
     
     Method to write to a file.
     !!! Beta, this method work only if fileID have been created and not closed!!!
     
     \param  fileID: file where to write data.
     \param  buffer: pointer to the buffer to add to the file fileID.
     \param  size: size of "buffer", in byte. 
     */
    void writeBuffer(ioserver_file fileID,char * buffer, int size);
  
    
    //2 steps initialisation :
    
    /*!
     Initialize the I/O server, this method is called by parallel.initialize(...). Should never be used!!!
     
     */
    void initialize(int proc_size0,int proc_size1, int IOserver_size, int IO_node_size); //called by avery cores.... initialize global variable
    
    
private:
    
    bool serverOn_flag;
    bool serverReady_flag;
    bool ostreamFile_flag;
    
    MPI_Group world_group_;
    
    MPI_Group IO_Group_;
    MPI_Group computeGroup_;
    MPI_Comm IO_Comm_;
    MPI_Comm computeComm_;
    
    MPI_Group syncLineGroup_; // root IO and compute 
    MPI_Comm  syncLineComm_;
    
    
    MPI_Group masterClientGroup_;
    MPI_Comm  masterClientComm_;
    
    MPI_Group IO_NodeGroup_;
    MPI_Comm IO_NodeComm_; 
    
    
    int IO_Rank_;
    int computeRank_;
    int syncLineRank_;
    int IO_NodeRank_;

    file_struct * files;
    
    
    int IO_ClientSize_;
    int IO_NodeSize_;
    int IO_Node_;
    
    MPI_Request sendRequest;
    
protected:
    
    char * dataBuffer;
    

};


#endif
