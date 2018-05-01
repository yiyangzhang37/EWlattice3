#include"LATfield2.h"


namespace LATfield2 {

	//CONSTANTS=====================

	int Lattice::initialized = 1;        //Status flag for initialized

	//CONSTRUCTORS===================

	Lattice::Lattice()
	{
		status_ = 0;
		arch_saved_ = false;
	}

	Lattice::Lattice(int dim, const int* size, int halo)
	{
		status_ = 0;
		arch_saved_ = false;
		this->initialize(dim, size, halo);
	}


	Lattice::Lattice(int dim, const int size, int halo)
	{
		status_ = 0;
		arch_saved_ = false;
		int* sizeArray = new int[dim];
		for (int i = 0; i < dim; i++) { sizeArray[i] = size; }
		this->initialize(dim, sizeArray, halo);
		delete[] sizeArray;
	}



	//DESTRUCTOR=========================

	Lattice::~Lattice()
	{
		if ((status_ & initialized) > 0)
		{
			delete[] size_;
			delete[] sizeLocal_;
			delete[] jump_;
		}
	}
	//INITIALIZE=========================

	void Lattice::initialize(int dim, const int size, int halo)
	{
		int* sizeArray = new int[dim];
		for (int i = 0; i < dim; i++) { sizeArray[i] = size; }
		this->initialize(dim, sizeArray, halo);
	}

	void Lattice::initialize(int dim, const int* size, int halo)
	{
		int i;

		if ((status_ & initialized) > 0)
		{
			delete[] size_;
			delete[] sizeLocal_;
			delete[] jump_;
		}
		//Store input lattice properties
		dim_ = dim;
		size_ = new int[dim_];
		for (i = 0; i < dim_; i++) size_[i] = size[i];
		halo_ = halo;


		//Calculate local size
		sizeLocal_ = new int[dim_];
		sizeLocal_[dim_ - 1] = int(ceil((parallel.grid_size()[0] - parallel.grid_rank()[0])*size_[dim_ - 1] / float(parallel.grid_size()[0])));
		sizeLocal_[dim_ - 1] -= int(ceil((parallel.grid_size()[0] - parallel.grid_rank()[0] - 1)*size_[dim_ - 1] / float(parallel.grid_size()[0])));
		sizeLocal_[dim_ - 2] = int(ceil((parallel.grid_size()[1] - parallel.grid_rank()[1])*size_[dim_ - 2] / float(parallel.grid_size()[1])));
		sizeLocal_[dim_ - 2] -= int(ceil((parallel.grid_size()[1] - parallel.grid_rank()[1] - 1)*size_[dim_ - 2] / float(parallel.grid_size()[1])));

		for (i = 0; i < dim_ - 2; i++) sizeLocal_[i] = size_[i];

		//Calculate index jumps
		jump_ = new long[dim_];
		jump_[0] = 1;
		for (i = 1; i < dim_; i++) jump_[i] = jump_[i - 1] * (sizeLocal_[i - 1] + 2 * halo_);

		//Calculate number of sites in lattice
		sitesLocal_ = 1;
		sitesLocalGross_ = 1;
		for (i = 0; i < dim_; i++)
		{
			sitesLocal_ *= sizeLocal_[i];
			sitesLocalGross_ *= sizeLocal_[i] + (2 * halo_);
		}
		sites_ = sitesLocal_;
		sitesGross_ = sitesLocalGross_;
		parallel.sum(sites_);
		parallel.sum(sitesGross_);

		//Calculate index of first and last local sites on lattice
		siteFirst_ = 0;
		siteLast_ = sitesLocalGross_ - 1;
		for (i = 0; i < dim_; i++)
		{
			siteFirst_ += jump_[i] * halo_;
			siteLast_ -= jump_[i] * halo_;
		}

		////calculate coordSkip



		//Get each processor to tell the others in his dim0_group its local sizeLocal_[dim-1])
		int* sizes_dim0 = new int[parallel.grid_size()[0]];
		for (i = 0; i < parallel.grid_size()[0]; i++)
		{
			if (i == parallel.grid_rank()[0]) { sizes_dim0[i] = sizeLocal_[dim_ - 1]; }
			parallel.broadcast_dim0(sizes_dim0[i], i);
		}
		//Sum up sizes for the processors of less than or equal rank
		coordSkip_[0] = 0;
		for (i = 0; i < parallel.grid_rank()[0]; i++)coordSkip_[0] += sizes_dim0[i];



		//Get each processor to tell the others in his dim1_group its local sizeLocal_[dim-2])
		int* sizes_dim1 = new int[parallel.grid_size()[1]];
		for (i = 0; i < parallel.grid_size()[1]; i++)
		{
			if (i == parallel.grid_rank()[1]) { sizes_dim1[i] = sizeLocal_[dim_ - 2]; }
			parallel.broadcast_dim1(sizes_dim1[i], i);
		}
		//Sum up sizes for the processors of less than or equal rank
		coordSkip_[1] = 0;
		for (i = 0; i < parallel.grid_rank()[1]; i++)coordSkip_[1] += sizes_dim1[i];




		////calculate sitesSkip : sitesskip used for fastread , fastload (function witch need to be coded :-) )

		sitesSkip_ = coordSkip_[0];
		for (i = 0; i < dim_ - 1; i++)sitesSkip_ *= size_[i];

		long siteSkiptemp = coordSkip_[1] * sizeLocal_[dim_ - 1];
		for (i = 0; i < dim_ - 2; i++)siteSkiptemp *= size_[i];

		sitesSkip_ += siteSkiptemp;

		//calculate sitesSkip2d : 

		int* sizes1 = new int[parallel.grid_size()[1]];
		int* sizes0 = new int[parallel.grid_size()[0]];
		long* offset1 = new long[parallel.grid_size()[1]];
		long* offset0 = new long[parallel.grid_size()[0]];

		int b = 1;
		int n;
		for (i = 0; i < dim_ - 2; i++)b *= sizeLocal_[i];

		//calulate offset in dim-2
		for (i = 0; i < parallel.grid_size()[1]; i++)sizes1[i] = sizes_dim1[i] * b;
		for (n = 0; n < parallel.grid_size()[1]; n++)offset1[n] = 0;

		for (n = 1; n < parallel.grid_size()[1]; n++)
		{
			for (i = 0; i < n; i++)offset1[n] += sizes1[i];
		}

		//calulate offset in dim-1
		for (i = 0; i < parallel.grid_size()[0]; i++)sizes0[i] = size_[dim_ - 2] * sizes_dim0[i] * b;
		for (n = 0; n < parallel.grid_size()[0]; n++)offset0[n] = 0;

		for (n = 1; n < parallel.grid_size()[0]; n++)
		{
			for (i = 0; i < n; i++)offset0[n] += sizes0[i];
		}

		sitesSkip2d_ = offset0[parallel.grid_rank()[0]] + offset1[parallel.grid_rank()[1]];

		//Set status
		status_ = status_ | initialized;

		//Free memory
		delete[] sizes_dim0;
		delete[] sizes_dim1;


	}


	void Lattice::save_arch(const string filename)
	{
		fstream file;
		int p, i;

		for (p = 0; p < parallel.size(); p++)
		{
			if (parallel.rank() == p)
			{
				if (parallel.rank() == 0)
				{
					//Truncate file if first process
					file.open(filename.c_str(), fstream::out | fstream::trunc);
				}
				else
				{
					//Append to file if another process
					file.open(filename.c_str(), fstream::out | fstream::app);
				}
				if (!file.is_open())
				{
					cerr << "Latfield::Lattice::save_arch - Could not open file for writing" << endl;
					cerr << "Latfield::Lattice::save_arch - File: " << filename << endl;
					parallel.abortRequest();
				}

				if (parallel.rank() == 0)
				{
					file << "# Architerctur of the lattice" << endl;
					file << "# Number of dimension :" << endl;
					file << this->dim() << endl;
					file << "############################" << endl;
				}
				file << "############################" << endl;
				file << "#  Ranks of processor, world, dim-1, dim-2 :" << endl;
				file << parallel.rank() << " " << parallel.grid_rank()[0] << " " << parallel.grid_rank()[1] << endl;
				file << "#  Local size :" << endl;
				for (i = 0; i < this->dim(); i++)
				{
					file << this->sizeLocal()[i] << endl;
				}
				file << "############################" << endl;


				file.close();
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		arch_saved_ = true;
		if (parallel.rank() == 0)cout << "Architecture saved in : " << filename << endl;

	}

	//MISCELLANEOUS======================

	bool Lattice::is_arch_saved() { return arch_saved_; };
	const int  Lattice::dim() const { return dim_; };
	int* Lattice::size() { return size_; };
	const int  Lattice::size(const int i) const { return size_[i]; }
	const long  Lattice::sites() const { return sites_; }
	const long  Lattice::sitesGross() const { return sitesGross_; }
	const int  Lattice::halo() { return halo_; }

	int* Lattice::sizeLocal() const { return sizeLocal_; };
	const int  Lattice::sizeLocal(const int i) const { return sizeLocal_[i]; }
	const long  Lattice::sitesLocal() const { return sitesLocal_; }
	const long  Lattice::sitesLocalGross() const { return sitesLocalGross_; }

	const long  Lattice::jump(const int i) const { return jump_[i]; }
	const long  Lattice::sitesSkip() const { return sitesSkip_; }
	const long  Lattice::sitesSkip2d() const { return sitesSkip2d_; }
	long*  Lattice::coordSkip() { return coordSkip_; }
	const long Lattice::siteFirst() const { return siteFirst_; }
	const long Lattice::siteLast() const { return siteLast_; }

}