#ifndef LATFIELD2_SETTINGSFILE_HPP
#define LATFIELD2_SETTINGSFILE_HPP


/*! \file LATfield2_SettingsFile.hpp
 \brief LATfield2_SettingsFile.hpp contain the class SettingsFile definition.
 \author N. Bevis
 */ 

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

//CLASS PROTOTYPE======================


/*! \class SettingsFile  
 \brief A utility class designed to make reading in runtime parameter values easier.
 

 If the command-line arguments are input via optional inputs on
 either the constructor or open member function, then these
 take preceident: they are effectively last in the file. 
 
 Note, when used with std::string objects, only one word
 is allowed per setting, ie. spaces are not allowed. This
 is because of the way that the >> operator works for
 this class. This fits nicely with the command-line override, however.
 
 Note that the string specified followed by = is searched
 for in the file and then the input read. If one setting
 name is also the end of another that precedes it in the file
 then the wrong one will be read.
 
 Only the primary MPI process is able to create or add to the setting
 file. Further processes will be sent the file contents via
 MPI. To use this class in serial code the preprocessor defintion SERIAL must be set. This flag have not been remove to allow users to use it outside LATfield2. LATfield2 have no serial version, therefor setting prepocessor flag -DSERIAL should never be used with LATfield2.
 
 */
class SettingsFile
{
private:
  std::string filename_;
  std::fstream file_;
  std::stringstream stream_;
  int mode_;
  bool isRoot_;    //is process the root one (false if non-parallel)                     
  //Search=================================
  bool search(const std::string search_string);
	
public:
  static int noCreate;
  static int autoCreate;


	//Constructors=======================
	//! Constructor
    SettingsFile();
    
    /*!
     Constructor + open a file
     \param filename : path to the file.
     \param mode     : noCreate (the read method will exit if the parameter does not exist) or autoCreate (read will add the missing parameter).
     \param argc     : additionnal argument number.
     \param argv     : pointer to the additionnal arguments.
     */
	SettingsFile(const std::string filename, const int mode, const int argc = 0, char** argv = NULL);
	
	//Destructor=========================
	//! desctructor
    ~SettingsFile();
	
	//File open / close / create ========
	/*!
     Open an existinge settings file
     \param filename : path to the file
     \param mode     : noCreate (the read method will exit if the parameter does not exist) or autoCreate (read will add the missing parameter).
     \param argc     : additional argument number.
     \param argv     : pointer to the additional arguments.
     */
    void open(const std::string filename, const int mode, const int argc = 0, char** argv = NULL);
    /*!
     Close the current settings file
     */
	void close();
    
    /*!
     Create a new settings file and open it.
     \param filename: path to the file.
     */
	void create(const std::string filename);
	
	//Settings read / write==================
    
    /*!
     Method to read a parameter.
     
     \param parameter_name : string containing the name of the parameter. If the parameter does not exite and the mode autocreate is set, this method will add the parameter to the settings file with the current value of "parameter". In the case the mode is set to nocreate, then read will exit for security, for this reason in it is always advise to set the mode to nocreate for production runs.
     \param parameter     : pointer to the variable where the parameter will be assigned. 
     */
	template<class TemplateClass>
	void read(const std::string parameter_name, TemplateClass& parameter);
    
    /*!
     Method to add a parameter to the settings file. The new parameter will be just added to the end of the file, even if it already exists.
     
     \param parameter_name: string containing the name of the parameter.
     \param parameter: pointer to the value of the parameter.
     */
	template<class TemplateClass>
	void add(const std::string parameter_name, const TemplateClass& parameter);
	
    /*!
     Method to write a parameter in the settings file. If the parameter_name exist, it will overwrite the parameter. And if it does not exist in the file, it will be added at the end of the file.
     
     \param parameter_name: string containing the name of the parameter
     \param parameter: pointer to the value of the parameter.
     */
    template<class TemplateClass>
	void write(const std::string parameter_name, const TemplateClass& parameter);

};

#endif
