#include"LATfield2.h"
namespace LATfield2{

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
			((this->coord(direction) + 1) % lattice_->size(direction));
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
		long tmp = (this->coord(direction) + steps >= 0) ? this->coord(direction) + steps : this->coord(direction) + steps + lattice_->size(direction);
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
			this->setIndex(this->lattice_->siteLast() + 1);
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
			return index_ % lattice_->jump(1) - lattice_->halo();
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
			return index_ % lattice_->jump(1) - lattice_->halo();
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
			return index_ % lattice_->jump(1) - lattice_->halo();
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
			return index_ % lattice_->jump(1) - lattice_->halo();
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
			return index_ % lattice_->jump(1) - lattice_->halo();
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
