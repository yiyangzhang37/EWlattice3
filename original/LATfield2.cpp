#include "LATfield2.hpp"

namespace LATfield2{

	//int2string.hpp
	string int2string(int number, int max, bool zeropad)
	{
		string output;
		char c;
		int i;

		//Get number of digits needed for max
		int digits = 1;
		for (i = 10; (i - 1)<max; i *= 10) digits++;

		//Create string of length digits (padded with zeros)
		for (i = 0; i<digits; i++)
		{
			c = '0' + number % 10;
			if (c != '0' || zeropad)
			{
				output = c + output;
			}
			number /= 10;
		}

		return output;
	}

	//LATfield2_Lattice.hpp

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
		for (int i = 0; i<dim; i++) { sizeArray[i] = size; }
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
		for (int i = 0; i<dim; i++) { sizeArray[i] = size; }
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
		for (i = 0; i<dim_; i++) size_[i] = size[i];
		halo_ = halo;


		//Calculate local size
		sizeLocal_ = new int[dim_];
		sizeLocal_[dim_ - 1] = int(ceil((parallel.grid_size()[0] - parallel.grid_rank()[0])*size_[dim_ - 1] / float(parallel.grid_size()[0])));
		sizeLocal_[dim_ - 1] -= int(ceil((parallel.grid_size()[0] - parallel.grid_rank()[0] - 1)*size_[dim_ - 1] / float(parallel.grid_size()[0])));
		sizeLocal_[dim_ - 2] = int(ceil((parallel.grid_size()[1] - parallel.grid_rank()[1])*size_[dim_ - 2] / float(parallel.grid_size()[1])));
		sizeLocal_[dim_ - 2] -= int(ceil((parallel.grid_size()[1] - parallel.grid_rank()[1] - 1)*size_[dim_ - 2] / float(parallel.grid_size()[1])));

		for (i = 0; i<dim_ - 2; i++) sizeLocal_[i] = size_[i];

		//Calculate index jumps
		jump_ = new long[dim_];
		jump_[0] = 1;
		for (i = 1; i<dim_; i++) jump_[i] = jump_[i - 1] * (sizeLocal_[i - 1] + 2 * halo_);

		//Calculate number of sites in lattice
		sitesLocal_ = 1;
		sitesLocalGross_ = 1;
		for (i = 0; i<dim_; i++)
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
		for (i = 0; i<dim_; i++)
		{
			siteFirst_ += jump_[i] * halo_;
			siteLast_ -= jump_[i] * halo_;
		}

		////calculate coordSkip



		//Get each processor to tell the others in his dim0_group its local sizeLocal_[dim-1])
		int* sizes_dim0 = new int[parallel.grid_size()[0]];
		for (i = 0; i<parallel.grid_size()[0]; i++)
		{
			if (i == parallel.grid_rank()[0]) { sizes_dim0[i] = sizeLocal_[dim_ - 1]; }
			parallel.broadcast_dim0(sizes_dim0[i], i);
		}
		//Sum up sizes for the processors of less than or equal rank
		coordSkip_[0] = 0;
		for (i = 0; i<parallel.grid_rank()[0]; i++)coordSkip_[0] += sizes_dim0[i];



		//Get each processor to tell the others in his dim1_group its local sizeLocal_[dim-2])
		int* sizes_dim1 = new int[parallel.grid_size()[1]];
		for (i = 0; i<parallel.grid_size()[1]; i++)
		{
			if (i == parallel.grid_rank()[1]) { sizes_dim1[i] = sizeLocal_[dim_ - 2]; }
			parallel.broadcast_dim1(sizes_dim1[i], i);
		}
		//Sum up sizes for the processors of less than or equal rank
		coordSkip_[1] = 0;
		for (i = 0; i<parallel.grid_rank()[1]; i++)coordSkip_[1] += sizes_dim1[i];




		////calculate sitesSkip : sitesskip used for fastread , fastload (function witch need to be coded :-) )

		sitesSkip_ = coordSkip_[0];
		for (i = 0; i<dim_ - 1; i++)sitesSkip_ *= size_[i];

		long siteSkiptemp = coordSkip_[1] * sizeLocal_[dim_ - 1];
		for (i = 0; i<dim_ - 2; i++)siteSkiptemp *= size_[i];

		sitesSkip_ += siteSkiptemp;

		//calculate sitesSkip2d : 

		int* sizes1 = new int[parallel.grid_size()[1]];
		int* sizes0 = new int[parallel.grid_size()[0]];
		long* offset1 = new long[parallel.grid_size()[1]];
		long* offset0 = new long[parallel.grid_size()[0]];

		int b = 1;
		int n;
		for (i = 0; i<dim_ - 2; i++)b *= sizeLocal_[i];

		//calulate offset in dim-2
		for (i = 0; i<parallel.grid_size()[1]; i++)sizes1[i] = sizes_dim1[i] * b;
		for (n = 0; n<parallel.grid_size()[1]; n++)offset1[n] = 0;

		for (n = 1; n<parallel.grid_size()[1]; n++)
		{
			for (i = 0; i<n; i++)offset1[n] += sizes1[i];
		}

		//calulate offset in dim-1
		for (i = 0; i<parallel.grid_size()[0]; i++)sizes0[i] = size_[dim_ - 2] * sizes_dim0[i] * b;
		for (n = 0; n<parallel.grid_size()[0]; n++)offset0[n] = 0;

		for (n = 1; n<parallel.grid_size()[0]; n++)
		{
			for (i = 0; i<n; i++)offset0[n] += sizes0[i];
		}

		sitesSkip2d_ = offset0[parallel.grid_rank()[0]] + offset1[parallel.grid_rank()[1]];

		//Set status
		status_ = status_ | initialized;

		//Free memory
		delete[] sizes_dim0;
		delete[] sizes_dim1;


	}
	#ifdef FFT3D
	void Lattice::initializeRealFFT(Lattice & lat_real, int halo)
	{

		if (lat_real.dim() != 3)
		{
			if (parallel.isRoot())
			{
				cerr << "Latfield2d::Lattice::initializeRealFFT : fft curently work only for 3d cubic lattice" << endl;
				cerr << "Latfield2d::Lattice::initializeRealFFT : coordinate lattice have not 3 dimensions" << endl;
				cerr << "Latfield2d : Abort Process Requested" << endl;

			}
			parallel.abortForce();
		}

		if (lat_real.size(0) != lat_real.size(1) | lat_real.size(2) != lat_real.size(1))
		{
			if (parallel.isRoot())
			{
				cerr << "Latfield2d::Lattice::initializeRealFFT : fft curently work only for 3d cubic lattice" << endl;
				cerr << "Latfield2d::Lattice::initializeRealFFT : coordinate lattice is not cubic" << endl;
				cerr << "Latfield2d : Abort Process Requested" << endl;

			}
			parallel.abortForce();
		}

		int lat_size[3];

		lat_size[0] = lat_real.size(0) / 2 + 1;
		lat_size[1] = lat_real.size(0);
		lat_size[2] = lat_real.size(0);

		this->initialize(3, lat_size, halo);
	}
	void Lattice::initializeComplexFFT(Lattice & lat_real, int halo)
	{

		if (lat_real.dim() != 3)
		{
			if (parallel.isRoot())
			{
				cerr << "Latfield2d::Lattice::initializeRealFFT : fft curently work only for 3d cubic lattice" << endl;
				cerr << "Latfield2d::Lattice::initializeRealFFT : coordinate lattice have not 3 dimensions" << endl;
				cerr << "Latfield2d : Abort Process Requested" << endl;

			}
			parallel.abortForce();
		}

		if (lat_real.size(0) != lat_real.size(1) | lat_real.size(2) != lat_real.size(1))
		{
			if (parallel.isRoot())
			{
				cerr << "Latfield2d::Lattice::initializeRealFFT : fft curently work only for 3d cubic lattice" << endl;
				cerr << "Latfield2d::Lattice::initializeRealFFT : coordinate lattice is not cubic" << endl;
				cerr << "Latfield2d : Abort Process Requested" << endl;

			}
			parallel.abortForce();
		}

		int lat_size[3];

		lat_size[0] = lat_real.size(0);
		lat_size[1] = lat_real.size(0);
		lat_size[2] = lat_real.size(0);

		this->initialize(3, lat_size, halo);
	}
	#endif


	void Lattice::save_arch(const string filename)
	{
		fstream file;
		int p, i;

		for (p = 0; p<parallel.size(); p++)
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
				for (i = 0; i<this->dim(); i++)
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
	const int  Lattice::size(const int i) const{ return size_[i]; }
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


	//LATfield2_Site.hpp
	//CONSTRUCTORS===================

	Site::Site() { ; }
	Site::Site(Lattice& lattice) { initialize(lattice); }
	Site::Site(Lattice& lattice, long index) { initialize(lattice, index); }

	//INITIALIZATION=================

	void Site::initialize(Lattice& lattice) { lattice_ = &lattice; index_ = 0l; }
	void Site::initialize(Lattice& lattice, long index) { lattice_ = &lattice; index_ = index; }

	//NEIGHBOURING SITE OPERATORS==

	//implement corrent periodic boundary condition; only valid for 3D 02/09/2017
	//This is not exactly true when gauge fields are considered.
	//Ideally, the (N-1) point should be identified with the 0 point; rather than the 0 point being the one-step-ahead of the (N-1) point.
	/*
	Site Site::operator+(const int direction) const
	{
		if (this->coord(direction) == lattice_->size(direction)-1) {
			Site s = *this;
			int newc[3];
			for (int i = 0; i < lattice_->dim(); ++i) {
				if (i == direction) newc[i] = 0;
				else newc[i] = this->coord(i); 
			}
			s.setCoord(newc);
			return s;
		}
		else {
			return Site(*lattice_, index_ + lattice_->jump(direction));
		}
	}
	*/
	/*
	Site Site::operator-(const int direction) const
	{
		if (this->coord(direction) == 0) {
			Site s = *this;
			int newc[3];
			for (int i = 0; i < lattice_->dim(); ++i) {
				if (i == direction) newc[i] = lattice_->size(i) - 1;
				else newc[i] = this->coord(i);
			}
			s.setCoord(newc);
			return s;
		}
		else {
			return Site(*lattice_, index_ - lattice_->jump(direction));
		}
	}
	*/

	/*New implementation of Site-moving operators*/
	Site Site::operator+(const int direction) const {
		long index = this->index_ - lattice_->jump(direction)*this->coord(direction);
		index += lattice_->jump(direction)*
		( (this->coord(direction) + 1) % lattice_->size(direction) );
		return Site(*lattice_, index);
	}

	Site Site::operator-(const int direction) const
	{
		long index = this->index_ - lattice_->jump(direction)*this->coord(direction);
		long tmp = (this->coord(direction) - 1 >= 0) ? this->coord(direction) - 1 : this->coord(direction) - 1 + lattice_->size(direction);
		index += lattice_->jump(direction)*
			(tmp % lattice_->size(direction));
		return Site(*lattice_, index);
	}

	Site Site::move(const int direction, const int steps) {
		//TODO
		long index = this->index_ - lattice_->jump(direction)*this->coord(direction);
		long tmp = (this->coord(direction) + steps >= 0) ? this->coord(direction) + steps : this->coord(direction) +steps + lattice_->size(direction);
		index += lattice_->jump(direction)*
			(tmp % lattice_->size(direction));
		this->index_ = index;
		return *this;
	}

	//LOOPING OPERATIONS====================

	void Site::first() { index_ = lattice_->siteFirst(); }

	bool Site::test() { return index_ <= lattice_->siteLast(); }

	void Site::next()
	{
		index_++;
		//If coordLocal(0) != sizeLocal(0) then next site reached
		if (coordLocal(0) != lattice_->sizeLocal(0)) { return; }
		else
		{
			index_ -= lattice_->sizeLocal(0);
			for (int i = 1; i<lattice_->dim(); i++)
			{
				index_ += lattice_->jump(i);
				//If coordLocal(i) !=sizeLocal(0) then next site reached
				if (coordLocal(i) != lattice_->sizeLocal(i)) { return; }
				index_ -= lattice_->sizeLocal(i) * lattice_->jump(i);
			}
			index_ = lattice_->siteLast() + 1;
		}
	}

	void Site::first(const int* rule) {
		int fc[3]; //coordinate of first element
		for (int i = 0; i < 3; ++i) {
			fc[rule[2 * i]] = rule[2 * i + 1] == 1 ? 0 : lattice_->size(rule[2 * i]) - 1;
		}
		this->setCoord(fc);
		return;
	}

	bool Site::test(const int* rule) {
		return true;
	}

	void Site::next(const int* rule) {
		const int nowc[3] = { this->coord(0),this->coord(1),this->coord(2) };
		int nxc[3] = { this->coord(0),this->coord(1),this->coord(2) };
		int dir[3] = { rule[0], rule[2], rule[4] };

		nxc[dir[0]] = nowc[dir[0]] + rule[1];
		for (int i = 0; i < 2; ++i) {
			if (nxc[dir[i]] >= this->lattice_->size(dir[i]) && rule[2 * i + 1] == 1) {
				nxc[dir[i]] = 0;
				nxc[dir[i + 1]] = nowc[dir[i + 1]] + rule[2 * i + 3];
			}
			if (nxc[dir[i]] < 0 && rule[2 * i + 1] == -1) {
				nxc[dir[i]] = this->lattice_->size(dir[i]) - 1;
				nxc[dir[i + 1]] = nowc[dir[i + 1]] + rule[2 * i + 3];
			}
		}
		if (nxc[2] >= this->lattice_->size(dir[2]) || nxc[2] < 0) {
			this->setIndex(this->lattice_->siteLast() +1 );
			return;
		}
		else {
			this->setCoord(nxc);
		}
		return;
	}

	//HALO OPERATIONS====================

	void Site::haloFirst() { index_ = 0; }
	bool Site::haloTest() { return index_ < lattice_->sitesLocalGross(); }

	void Site::haloNext()
	{
		index_++;

		//Can only leave boundary by reaching coord(0)=0, so otherwise done
		if (coordLocal(0) == 0)
		{
			bool is_in_halo = 0;
			for (int i = 1; i<lattice_->dim(); i++)
			{
				if (coordLocal(i)<0 || coordLocal(i) >= lattice_->sizeLocal(i)) { is_in_halo = 1; break; }
			}
			if (!is_in_halo) index_ += lattice_->sizeLocal(0);
		}
	}

	//INDEX ADVANCE======================

	void Site::indexAdvance(const long number) { index_ += number; }

	//MISCELLANEOUS======================

	long Site::index() const { return index_; }

	void Site::setIndex(const long new_index) { index_ = new_index; }

	int Site::coord(const int direction) const ////////sensible a quelle dim est scatter (seul modif a faire ici)
	{
		if (direction<lattice_->dim() - 2) { return coordLocal(direction); }
		else if (direction == lattice_->dim() - 2) { return coordLocal(direction) + lattice_->coordSkip()[1]; }
		else { return coordLocal(direction) + lattice_->coordSkip()[0]; }
	}

	int Site::coordLocal(const int direction) const
	{
		if (direction == lattice_->dim() - 1)
		{
			return index_ / lattice_->jump(direction) - lattice_->halo();
		}
		else if (direction == 0)
		{
			return index_%lattice_->jump(1) - lattice_->halo();
		}
		else
		{
			return (index_%lattice_->jump(direction + 1)) / lattice_->jump(direction) - lattice_->halo();
		}
	}

	bool Site::setCoord(int* r)
	{

		this->first();
		//Check site is local
		if (r[lattice_->dim() - 1]<this->coord(lattice_->dim() - 1) || r[lattice_->dim() - 1] >= this->coord(lattice_->dim() - 1) + lattice_->sizeLocal(lattice_->dim() - 1)
			|| r[lattice_->dim() - 2]<this->coord(lattice_->dim() - 2) || r[lattice_->dim() - 2] >= this->coord(lattice_->dim() - 2) + lattice_->sizeLocal(lattice_->dim() - 2))
		{
			return false;
			//COUT<<"LATfield::Site::setCoord(int*) - Site desired non-local!"<<endl;
		}
		else
		{

			int jump = 0;
			for (int i = 0; i<lattice_->dim(); i++)
			{
				jump += (r[i] - coord(i))*lattice_->jump(i);
			}

			this->indexAdvance(jump);
			return true;
		}
	}

	bool Site::setCoord(int x, int y = 0, int z = 0)
	{
		int* r = new int[3];
		r[0] = x;
		r[1] = y;
		r[2] = z;
		return this->setCoord(r);
		delete[] r;
	}
	Lattice& Site::lattice() { return *lattice_; }

	#ifdef FFT3D
	//class cKSite
	void cKSite::initialize(Lattice& lattice) { lattice_ = &lattice; directions_[0] = 1; directions_[1] = 2; directions_[2] = 0; }
	void cKSite::initialize(Lattice& lattice, long index) { lattice_ = &lattice; index_ = index; directions_[0] = 1; directions_[1] = 2; directions_[2] = 0; }

	cKSite cKSite::operator+(int asked_direction)
	{
		return cKSite(*lattice_, index_ + lattice_->jump(directions_[asked_direction]));
	}

	cKSite cKSite::operator-(int asked_direction)
	{
		return cKSite(*lattice_, index_ - lattice_->jump(directions_[asked_direction]));
	}
	int cKSite::coord(int asked_direction) ////////sensible a quelle dim est scatter (seul modif a faire ici)
	{
		int direction = directions_[asked_direction];

		if (direction<lattice_->dim() - 2) { return latCoordLocal(direction); }
		else if (direction == lattice_->dim() - 2) { return latCoordLocal(direction) + lattice_->coordSkip()[1]; }
		else { return latCoordLocal(direction) + lattice_->coordSkip()[0]; }
	}

	int cKSite::latCoord(int direction) ////////sensible a quelle dim est scatter (seul modif a faire ici)
	{

		if (direction<lattice_->dim() - 2) { return coordLocal(direction); }
		else if (direction == lattice_->dim() - 2) { return coordLocal(direction) + lattice_->coordSkip()[1]; }
		else { return coordLocal(direction) + lattice_->coordSkip()[0]; }
	}

	int cKSite::coordLocal(int asked_direction)
	{
		int direction = directions_[asked_direction];

		if (direction == lattice_->dim() - 1)
		{
			return index_ / lattice_->jump(direction) - lattice_->halo();
		}
		else if (direction == 0)
		{
			return index_%lattice_->jump(1) - lattice_->halo();
		}
		else
		{
			return (index_%lattice_->jump(direction + 1)) / lattice_->jump(direction) - lattice_->halo();
		}
	}
	int cKSite::latCoordLocal(int direction)
	{

		if (direction == lattice_->dim() - 1)
		{
			return index_ / lattice_->jump(direction) - lattice_->halo();
		}
		else if (direction == 0)
		{
			return index_%lattice_->jump(1) - lattice_->halo();
		}
		else
		{
			return (index_%lattice_->jump(direction + 1)) / lattice_->jump(direction) - lattice_->halo();
		}
	}

	bool cKSite::setCoord(int* r_asked)
	{
		int r[3];
		r[0] = r_asked[2];
		r[1] = r_asked[0];
		r[2] = r_asked[1];
		this->first();
		//Check site is local
		if (r[lattice_->dim() - 1]<this->latCoord(lattice_->dim() - 1) || r[lattice_->dim() - 1] >= this->latCoord(lattice_->dim() - 1) + lattice_->sizeLocal(lattice_->dim() - 1)
			|| r[lattice_->dim() - 2]<this->latCoord(lattice_->dim() - 2) || r[lattice_->dim() - 2] >= this->latCoord(lattice_->dim() - 2) + lattice_->sizeLocal(lattice_->dim() - 2))
		{
			return false;
			//COUT<<"LATfield::Site::setCoord(int*) - Site desired non-local!"<<endl;
		}
		else
		{

			int jump = 0;
			for (int i = 0; i<lattice_->dim(); i++)
			{
				jump += (r[i] - latCoord(i))*lattice_->jump(i);
			}

			this->indexAdvance(jump);
			return true;
		}
	}
	bool cKSite::setCoord(int x, int y = 0, int z = 0)
	{
		int* r = new int[3];
		r[0] = x;
		r[1] = y;
		r[2] = z;
		return this->setCoord(r);
		delete[] r;
	}

	//class rKSite
	void rKSite::initialize(Lattice& lattice) { lattice_ = &lattice; directions_[0] = 0; directions_[1] = 2; directions_[2] = 1; }
	void rKSite::initialize(Lattice& lattice, long index) { lattice_ = &lattice; index_ = index; directions_[0] = 0; directions_[1] = 2; directions_[2] = 1; }

	rKSite rKSite::operator+(int asked_direction)
	{
		return rKSite(*lattice_, index_ + lattice_->jump(directions_[asked_direction]));
	}

	rKSite rKSite::operator-(int asked_direction)
	{
		return rKSite(*lattice_, index_ - lattice_->jump(directions_[asked_direction]));
	}
	int rKSite::coord(int asked_direction) ////////sensible a quelle dim est scatter (seul modif a faire ici)
	{
		int direction = directions_[asked_direction];

		if (direction<lattice_->dim() - 2) { return latCoordLocal(direction); }
		else if (direction == lattice_->dim() - 2) { return latCoordLocal(direction) + lattice_->coordSkip()[1]; }
		else { return latCoordLocal(direction) + lattice_->coordSkip()[0]; }
	}

	int rKSite::latCoord(int direction) ////////sensible a quelle dim est scatter (seul modif a faire ici)
	{

		if (direction<lattice_->dim() - 2) { return coordLocal(direction); }
		else if (direction == lattice_->dim() - 2) { return coordLocal(direction) + lattice_->coordSkip()[1]; }
		else { return coordLocal(direction) + lattice_->coordSkip()[0]; }
	}

	int rKSite::coordLocal(int asked_direction)
	{
		int direction = directions_[asked_direction];

		if (direction == lattice_->dim() - 1)
		{
			return index_ / lattice_->jump(direction) - lattice_->halo();
		}
		else if (direction == 0)
		{
			return index_%lattice_->jump(1) - lattice_->halo();
		}
		else
		{
			return (index_%lattice_->jump(direction + 1)) / lattice_->jump(direction) - lattice_->halo();
		}
	}
	int rKSite::latCoordLocal(int direction)
	{

		if (direction == lattice_->dim() - 1)
		{
			return index_ / lattice_->jump(direction) - lattice_->halo();
		}
		else if (direction == 0)
		{
			return index_%lattice_->jump(1) - lattice_->halo();
		}
		else
		{
			return (index_%lattice_->jump(direction + 1)) / lattice_->jump(direction) - lattice_->halo();
		}
	}

	bool rKSite::setCoord(int* r_asked)
	{
		int r[3];
		r[0] = r_asked[0];
		r[1] = r_asked[2];
		r[2] = r_asked[1];
		this->first();
		//Check site is local
		if (r[lattice_->dim() - 1]<this->latCoord(lattice_->dim() - 1) || r[lattice_->dim() - 1] >= this->latCoord(lattice_->dim() - 1) + lattice_->sizeLocal(lattice_->dim() - 1)
			|| r[lattice_->dim() - 2]<this->latCoord(lattice_->dim() - 2) || r[lattice_->dim() - 2] >= this->latCoord(lattice_->dim() - 2) + lattice_->sizeLocal(lattice_->dim() - 2))
		{
			return false;
			//COUT<<"LATfield::Site::setCoord(int*) - Site desired non-local!"<<endl;
		}
		else
		{

			int jump = 0;
			for (int i = 0; i<lattice_->dim(); i++)
			{
				jump += (r[i] - latCoord(i))*lattice_->jump(i);
			}

			this->indexAdvance(jump);
			return true;
		}
	}
	bool rKSite::setCoord(int x, int y = 0, int z = 0)
	{
		int* r = new int[3];
		r[0] = x;
		r[1] = y;
		r[2] = z;
		return this->setCoord(r);
		delete[] r;
	}
	#endif

}

// LATfield2_SettingsFile.hpp

//CONSTANTS===========================
int SettingsFile::noCreate = 1;
int SettingsFile::autoCreate = 0;


//CONSTRUCTORS========================
SettingsFile::SettingsFile()
{

#ifndef SERIAL
	isRoot_ = parallel.isRoot();
#else
	isRoot_ = true;
#endif

}

SettingsFile::SettingsFile(const std::string filename, const int mode, const int argc, char** argv)
{

#ifndef SERIAL
	isRoot_ = parallel.isRoot();
#else
	isRoot_ = true;
#endif

	this->open(filename, mode, argc, argv);
}

//DESTRUCTOR==========================
SettingsFile::~SettingsFile() { this->close(); }

//OPEN================================
void SettingsFile::open(const std::string filename, const int mode, const int argc, char** argv)
{
	char c;

	filename_ = filename;
	mode_ = mode;

	if (isRoot_)
	{
		//Open file  
		file_.open(filename_.c_str(), std::fstream::in);
		if (!file_.is_open())
		{
			if ((mode_ & SettingsFile::noCreate) == 0)
			{
				std::cout << "SettingsFile: " << filename_ << " not found." << std::endl;
				std::cout << "SettingsFile: Creating..." << std::endl;
				this->create(filename_);
				std::cout << "creating ok" << std::endl;
			}
			else
			{
				std::cout << "SettingsFile: " << filename_ << " not found and auto create off." << std::endl;
				std::cout << "SettingsFile: Exiting..." << std::endl;
#ifndef SERIAL
				parallel.abortRequest();
#else
				exit(555);
#endif
			}
		}

		//Read command line into stringstream
		for (int i = 0; i<argc; i++)
		{
			for (int j = 0; argv[i][j] != '\0'; j++)
			{
				stream_ << argv[i][j];
			}
			stream_ << endl;
		}

		//Read file into stringstream
		while (!file_.eof())
		{
			c = file_.get();
			if (c == '#')
			{
				while (!file_.eof() && c != '\n') { c = file_.get(); }
				if (file_.eof()) { break; }
			}
			stream_.put(c);
		}
		file_.close();
	}


#ifndef SERIAL
	//Broadcast results to all processes 


	parallel.barrier();

	if (parallel.size()>1)
	{
		if (parallel.isRoot())
		{
			int len = stream_.str().length();
			char* streamString = new char[len + 1];
			for (int i = 0; i <= len; i++) { streamString[i] = stream_.str()[i]; }
			parallel.broadcast(len, parallel.root());
			parallel.broadcast(streamString, len + 1, parallel.root());
		}
		else
		{
			int len;
			char* streamString;
			parallel.broadcast(len, parallel.root());
			streamString = new char[len + 1];
			parallel.broadcast(streamString, len + 1, parallel.root());
			stream_ << streamString;
		}
	}

#endif
}

//FILE CLOSE============================
void SettingsFile::close()
{
	if (isRoot_) { filename_ = "."; }
}

//FILE CREATE===========================
void SettingsFile::create(const std::string filename)
{

	if (isRoot_)
	{
		filename_ = filename;
		mode_ = autoCreate;

		file_.open(filename_.c_str(), std::fstream::out);
		if (!file_.is_open())
		{
			std::cout << "SettingsFile: Cannot create: " << filename << std::endl;
			std::cout << "SettingsFile: Exiting..." << std::endl;
#ifndef SERIAL
			parallel.abortRequest();
#else
			exit(555);
#endif	
		}
		else
		{
			file_.close();
			file_.clear();
			file_.open(filename.c_str(), std::fstream::in);
		}
	}

	//parallel.barrier();
}

//PARAMETER READ===========================
template<class TemplateClass>
void SettingsFile::read(const std::string parameterName, TemplateClass &parameter)
{
	if (this->search(parameterName + '='))
	{
		stream_ >> parameter;
	}
	else
	{
		if (isRoot_)
		{
			//verifiy that the parameter name is no autocreate
			if (parameterName != "autocreate")
			{
				std::cout << "SettingsFile: " << filename_ << " has no parameter: " << parameterName << std::endl;
				std::cout << "SettingsFile: No command-line override given either" << std::endl;

				if ((mode_ & SettingsFile::noCreate) == 0)
				{
					std::cout << "SettingsFile: Adding with current value: " << parameter << std::endl;
					this->add(parameterName, parameter);
				}
				else
				{
					std::cout << "SettingsFile: Auto create off. Exiting..." << std::endl;
#ifndef SERIAL
					parallel.abortRequest();
#else
					exit(555);
#endif
				}
			}

		}
#ifndef SERIAL
		parallel.barrier();
#endif
	}
}

//PARAMETER WRITE===========================
template<class TemplateClass>
void SettingsFile::add(const std::string parameter_name, const TemplateClass &parameter)
{
	if (isRoot_)
	{
		file_.clear();
		file_.open(filename_.c_str(), std::ios::out | std::ios::app);
		file_ << parameter_name << '=' << parameter << std::endl;
		if (!file_.good())
		{
			std::cout << "SettingsFile: Could not write to file: " << filename_ << std::endl;
			std::cout << "SettingsFile: Exiting... " << std::endl;
#ifndef SERIAL
			parallel.abortRequest();
#else
			exit(555);
#endif
		}
		file_.close();
	}

}
template<class TemplateClass>
void SettingsFile::write(const std::string parameter_name, const TemplateClass &parameter)
{
	if (isRoot_)
	{
		unsigned int i = 0;
		unsigned int l = 0;
		int line_number = 0;

		char c;
		string line_temp;

		if (file_.good())file_.close();
		file_.clear();
		file_.open(filename_.c_str(), std::ios::in);


		//get number of non void line
		file_.seekg(0);
		while (!file_.eof())
		{
			getline(file_, line_temp);
			line_number++;
		}

		//creat line array and this_param array
		string lines[line_number];
		bool this_param[line_number];
		file_.seekg(0);
		file_.close();
		file_.open(filename_.c_str(), std::ios::in);

		//implement non void lines
		l = 0;
		while (!file_.eof())
		{
			getline(file_, line_temp);

			lines[l] = line_temp;
			l++;

		}
		for (l = 0; l<line_number; l++)
		{
			for (i = 0; i<parameter_name.length(); i++)
			{
				c = lines[l][i];
				if (c == parameter_name[i])this_param[l] = true;
				else
				{
					this_param[l] = false;
					i = parameter_name.length();
				}
			}
		}
		file_.close();
		file_.open(filename_.c_str(), std::ios::out | std::ios::trunc);


		bool writen;
		for (l = 0; l<line_number; l++)
		{
			if (!this_param[l])file_ << lines[l] << endl;
			else
			{
				file_ << parameter_name << '=' << parameter << endl;
				writen = true;
			}
		}

		if (!writen)file_ << parameter_name << '=' << parameter << endl;


		file_.close();
		file_.clear();

	}

#ifndef SERIAL
	parallel.barrier();
#endif

}


//SEARCH=====================================
bool SettingsFile::search(const std::string searchString)
{
	unsigned int i = 0;
	char c;
	//Set to beginning of file
	stream_.clear(); //clear any errors from having failed a previous search
	stream_.seekg(0);

	//Search
	while (stream_.good() && i<searchString.length())
	{
		c = stream_.get();
		if (c == searchString[i]) { i++; }
		else { i = 0; }
	}
	return stream_.good();
}

// LATfield2_parallel2d.hpp
Parallel2d::Parallel2d()
{


	int argc = 1;
	char** argv = new char*[argc];
	for (int i = 0; i<argc; i++) { argv[i] = new char[20]; }
#ifndef EXTERNAL_IO	
	lat_world_comm_ = MPI_COMM_WORLD;
	world_comm_ = MPI_COMM_WORLD;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(lat_world_comm_, &lat_world_rank_);
	MPI_Comm_size(lat_world_comm_, &lat_world_size_);
	MPI_Comm_rank(world_comm_, &world_rank_);
	MPI_Comm_size(world_comm_, &world_size_);
#else
	world_comm_ = MPI_COMM_WORLD;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(world_comm_, &world_rank_);
	MPI_Comm_size(world_comm_, &world_size_);

#endif

}
#ifdef EXTERNAL_IO
void Parallel2d::initialize(int proc_size0, int proc_size1, int IO_total_size, int IO_node_size)
#else
void Parallel2d::initialize(int proc_size0, int proc_size1)
#endif
{

	grid_size_[0] = proc_size0;
	grid_size_[1] = proc_size1;

	dim0_comm_ = (MPI_Comm *)malloc(grid_size_[1] * sizeof(MPI_Comm));
	dim1_comm_ = (MPI_Comm *)malloc(grid_size_[0] * sizeof(MPI_Comm));

	dim0_group_ = (MPI_Group *)malloc(grid_size_[1] * sizeof(MPI_Group));
	dim1_group_ = (MPI_Group *)malloc(grid_size_[0] * sizeof(MPI_Group));

	int rang[3], i, j, comm_rank;



#ifdef EXTERNAL_IO

	if (world_rank_ == 0)
	{
		if (proc_size0*proc_size1 + IO_total_size != world_size_)
		{
			cerr << "Latfield::Parallel2d::initialization - wrong number of process" << endl;
			cerr << "Latfield::Parallel2d::initialization - Number of total process must be equal to proc_size0*proc_size1+IO_total_size" << endl;
			cerr << "Latfield::Parallel2d::initialization - Within the call : Parallel2d(int proc_size0, int proc_size1, int IO_total_size)" << endl;
			this->abortForce();
		}



	}

	MPI_Comm_group(world_comm_, &world_group_);

	rang[0] = 0;
	rang[1] = proc_size0*proc_size1 - 1;
	rang[2] = 1;

	MPI_Group_range_incl(world_group_, 1, &rang, &lat_world_group_);
	MPI_Comm_create(world_comm_, lat_world_group_, &lat_world_comm_);

	//lat_world_comm_ = MPI_COMM_WORLD;
	//MPI_Comm_group(lat_world_comm_,&lat_world_group_);


	MPI_Group_rank(lat_world_group_, &comm_rank);
	if (comm_rank != MPI_UNDEFINED)
	{
		lat_world_rank_ = comm_rank;
		MPI_Comm_size(lat_world_comm_, &lat_world_size_);

		rang[2] = 1;
		for (j = 0; j<grid_size_[1]; j++)
		{
			rang[0] = j * grid_size_[0];
			rang[1] = grid_size_[0] - 1 + j*grid_size_[0];
			MPI_Group_range_incl(lat_world_group_, 1, &rang, &dim0_group_[j]);
			MPI_Comm_create(lat_world_comm_, dim0_group_[j], &dim0_comm_[j]);
		}


		rang[2] = grid_size_[0];
		for (i = 0; i<grid_size_[0]; i++)
		{
			rang[0] = i;
			rang[1] = i + (grid_size_[1] - 1)*grid_size_[0];
			MPI_Group_range_incl(lat_world_group_, 1, &rang, &dim1_group_[i]);
			MPI_Comm_create(lat_world_comm_, dim1_group_[i], &dim1_comm_[i]);
		}


		for (i = 0; i<grid_size_[0]; i++)
		{
			MPI_Group_rank(dim1_group_[i], &comm_rank);
			if (comm_rank != MPI_UNDEFINED)grid_rank_[1] = comm_rank;
		}

		for (j = 0; j<grid_size_[1]; j++)
		{
			MPI_Group_rank(dim0_group_[j], &comm_rank);
			if (comm_rank != MPI_UNDEFINED)grid_rank_[0] = comm_rank;
		}


		root_ = 0;
		isIO_ = false;
	}
	else
	{
		lat_world_rank_ = -1;
		root_ = 0;
		grid_rank_[1] = -1;
		grid_rank_[0] = -1;
		isIO_ = true;
	}

	if (grid_rank_[0] == grid_size_[0] - 1)last_proc_[0] = true;
	else last_proc_[0] = false;
	if (grid_rank_[1] == grid_size_[1] - 1)last_proc_[1] = true;
	else last_proc_[1] = false;




	IO_Server.initialize(proc_size0, proc_size1, IO_total_size, IO_node_size);

#else

	if (lat_world_rank_ == 0)
	{
		if (proc_size0*proc_size1 != lat_world_size_)
		{
			cerr << "Latfield::Parallel2d::initialization - wrong number of process" << endl;
			cerr << "Latfield::Parallel2d::initialization - Number of total process must be equal to proc_size0*proc_size1" << endl;
			cerr << "Latfield::Parallel2d::initialization - Within the call : Parallel2d(int proc_size0, int proc_size1)" << endl;
			this->abortForce();
		}
	}

	MPI_Comm_group(lat_world_comm_, &lat_world_group_);




	rang[2] = 1;
	for (j = 0; j<grid_size_[1]; j++)
	{
		rang[0] = j * grid_size_[0];
		rang[1] = grid_size_[0] - 1 + j*grid_size_[0];
		MPI_Group_range_incl(lat_world_group_, 1, &rang, &dim0_group_[j]);
		MPI_Comm_create(lat_world_comm_, dim0_group_[j], &dim0_comm_[j]);
	}


	rang[2] = grid_size_[0];
	for (i = 0; i<grid_size_[0]; i++)
	{
		rang[0] = i;
		rang[1] = i + (grid_size_[1] - 1)*grid_size_[0];
		MPI_Group_range_incl(lat_world_group_, 1, &rang, &dim1_group_[i]);
		MPI_Comm_create(lat_world_comm_, dim1_group_[i], &dim1_comm_[i]);
	}



	for (i = 0; i<grid_size_[0]; i++)
	{
		MPI_Group_rank(dim1_group_[i], &comm_rank);
		if (comm_rank != MPI_UNDEFINED)grid_rank_[1] = comm_rank;
	}

	for (j = 0; j<grid_size_[1]; j++)
	{
		MPI_Group_rank(dim0_group_[j], &comm_rank);
		if (comm_rank != MPI_UNDEFINED)grid_rank_[0] = comm_rank;
	}


	root_ = 0;

	if (grid_rank_[0] == grid_size_[0] - 1)last_proc_[0] = true;
	else last_proc_[0] = false;
	if (grid_rank_[1] == grid_size_[1] - 1)last_proc_[1] = true;
	else last_proc_[1] = false;



#endif

	if (root_ == lat_world_rank_)isRoot_ = true;
	else isRoot_ = false;


}

Parallel2d::~Parallel2d()
{


	/*free(dim0_comm_);
	free(dim1_comm_);
	free(dim0_group_);
	free(dim1_group_);*/
	int finalized;
	MPI_Finalized(&finalized);
	if (!finalized) { MPI_Finalize(); }
}

//ABORT AND BARRIER===============================

void Parallel2d::abortForce()
{
	MPI_Abort(world_comm_, EXIT_FAILURE);
}

void Parallel2d::abortRequest()
{
	char failure;
	if (isRoot())
	{
		failure = char(1);
		broadcast(failure, root_);
		exit(EXIT_FAILURE);
	}
	else
	{
		cout << "Parallel::abortRequest() called from non-Root process." << endl;
		cout << "Parallel::abortRequest() calling abortForce()..." << endl;
		abortForce();
	}
}

void Parallel2d::barrier()
{

	MPI_Barrier(lat_world_comm_);
}