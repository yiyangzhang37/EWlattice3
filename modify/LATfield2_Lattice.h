#ifndef LATFIELD2_LATTICE_HPP
#define LATFIELD2_LATTICE_HPP


/*! \file LATfield2_Lattice.hpp
\brief LATfield2_Lattice.hpp contains the class Lattice definition.
\author David Daverio, Neil Bevis
*/



/*! \class Lattice

\brief The Lattice class describe a cartesian mesh (with 2 or more dimensions). The updateHalo method of the Field class generate the periodicity.


It store the global and local geometry of the mesh. The last 2 dimension of the lattice are scattered into the MPI processes grid.


*/
namespace LATfield2
{
	class Lattice
	{
	public:
		//! Constructor.
		Lattice();

		/*!
		Constructor with initialization
		\sa initialize(int dim, const int* size, int halo);
		\param dim : number of dimension
		\param size : array containing the size of each dimension.
		\param halo : size of the halo (ghost cells, same for each dimension)
		*/
		Lattice(int dim, const int* size, int halo);

		/*!
		Constructor with initialization
		\sa initialize(int dim, const int size, int halo);
		\param dim : number of dimension
		\param size : size of each dimension (same for each dimension)
		\param halo : size of the halo (same for each dimension)
		*/
		Lattice(int dim, const int size, int halo);

		//! Destructor.
		~Lattice();

		/*!
		Initialization of a dim-dimensional lattice, the size of each dimension is set by the second parameter: int *size. The ghost cell number (halo) is the same for each dimension.

		\param dim : number of dimension
		\param size : array containing the size of each dimension.
		\param halo : size of the halo (same for each dimension)
		*/
		void initialize(int dim, const int* size, int halo);

		/*!
		Initialization of a dim-dimensional lattice, each dimension have the same size. The ghost cell number (halo) is the same for each dimension.


		\param dim : number of dimension
		\param size : size of each dimension (same for each dimension)
		\param halo : size of the halo (same for each dimension)
		*/
		void initialize(int dim, const int size, int halo);

		/*!
		\return int. Number of dimensions of the lattice.
		*/
		const int  dim() const;
		/*!
		\return int. Size of the halo (ghost cells).
		*/
		const int  halo();

		/*!
		\return int*. Pointer to the array of the size of each dimension of the lattice.
		*/
		int* size();

		/*!
		Function which returns the size of a given dimension of the lattice.
		\param direction : asked dimension.
		\return int. Global size of the lattice in the given dimension.
		*/
		const int  size(const int direction) const;  //Size in a particular dimension

		/*!
		\return int*. Pointer to the array of the size of each dimension of the sublattice stored in this MPI process.
		*/
		int* sizeLocal() const;               //Local version

		/*!
		Function which returns the size of a given dimension of the sublattice stored in this MPI process.
		\param direction : asked dimension.
		\return int. Global size of the sublattice (of this MPI process) in the given dimension.
		*/
		const int  sizeLocal(const int direction) const;  //Local version

		/*!
		\return long. Number of sites on the lattice (excluding halo sites).
		*/
		const long  sites() const;              //Number of (global) sites
		/*!
		\return long. Number of sites on the lattice (including halo sites).
		*/
		const long  sitesGross() const;         //Number of (global) sites including halo

		/*!
		\return long. Number of sites (excluding halo sites) of the sublattice stored in this MPI process.
		*/
		const long  sitesLocal() const;              //Local version

		/*!
		\return long. Number of sites (including halo sites) of the sublattice stored in this MPI process.
		*/
		const long  sitesLocalGross() const;         //Local version

		/*!
		\return long. Array index of the first site which is not within the halo.
		*/
		const long siteFirst() const;

		/*!
		\return long. Array index of the last site which is not within the halo.
		*/
		const long siteLast() const;



		/*!
		Function which return the number of data_ array elements to jump to move to the next site in the given direction. (does not take into account the number of component of the fields, therefor should be multiplied by Field.components().) Should not be used by user.
		\param direction : asked direction.
		\return long. Number of array elements to jump.
		*/
		const long  jump(const int direction) const;       //Number of sites jumped to move in direction

		/*!
		\return long. Number of sites before first local site in lattice. Should not be used by users.
		*/
		const long  sitesSkip() const;               //Number of sites before first local site in lattice (say in a file)

		/*!
		\return long. Number of sites before first local site in lattice. Should not be used by users.
		*/
		const long  sitesSkip2d() const;

		/*!
		\return *long. Pointer to an array which store the last 2 dimensions coordinate of the first local(in this MPI process) sites. Index 0 is for dim-1, index 1 is for dim-2/
		*/
		long*  coordSkip();              //Number to add to coord[dim_-1] to get global value

		/*!
		Function which save in serial and in ASCII the global and local description of the Lattice. Usefull to read a file writen by fast_save or fast_write methods of the Field class.
		\param filename : filename of the architectur file.
		*/
		void save_arch(const string filename);

		/*!
		\return return true if the description of the lattice has been written on disk.

		\sa save_arch(const string filename)
		*/
		bool is_arch_saved();

	private:
		int        status_;
		static int initialized;

		//Global variables==============
		int  dim_;              //ok//Number of dimensions
		int* size_;             //ok//Number of lattice sites in each direction
		long  sites_;            //ok//Number of sites in lattice
		long  sitesGross_;       //ok//Number of sites in lattice plus halos
		int  halo_;             //ok//Number of sites extra in each direction

								//Local variables===============
		int* sizeLocal_;       //ok//Number of local lattice sites in each direction
		long  sitesLocal_;      //ok//Number of local sites in lattice
		long  sitesLocalGross_; //ok//Number of local sites in lattice plus halo
		long* jump_;            //ok//Jumps needed to move up one in each direction


		long  siteFirst_;       //ok//Index of first local site in lattice

		long  siteLast_;        //ok//Index of last local site in lattice
		long  sitesSkip_;      //Number of global lattice sites before first local site (say in a file)
		long  sitesSkip2d_;
		long  coordSkip_[2];       //Number to add to coord[dim_-1] and coord[dim_-2] to get global value

		//save variable for fast save
		int arch_saved_;

	};

}
#endif