#ifndef LATFIELD3_SITE_H
#define LATFIELD3_SITE_H

#include "LATfield3_Lattice.h"

namespace LATfield3{
    /*
    class Site:
    This class identifies where the current position is, w.r.t. the associated Lattice object.
    It also performs sequential looping over the local lattice, 
    or random accessing any local lattice sites.
    */
	template<int DIM>
    class Site{
    public:
        Site() = delete;
        Site(const Lattice<DIM>& lattice);
        ~Site();

		/*
		move the Site object to the first local visible site
		*/
        Site& first();

		/*
		move the Site object to the next local visible site.
		(assuming right now the object is at a local visible site)
		*/
        Site& next();

		/*
		check if the Site object is not beyond the last local 
		visible site
		*/
        bool test() const;

		/*
		move the Site object to the first local halo site.
		*/
        Site& halo_first();

		/*
		move the Site object to the next local halo site.
		(assuming right now the object is at a local halo site)
		*/
        Site& halo_next();

		/*
		check if the Site object is not beyond the last local halo site.
		*/
        bool halo_test() const;

		/*
		specific Site object move within the local_mem region,
		based on the current Site::index_.
		The object should not move out of the local_mem region.
		*/
        Site& move(const int direction, const int steps);
        Site operator+(const int direction) const;
        Site operator-(const int direction) const;
        Site& operator+(const int direction);
        Site& operator-(const int direction);
        void index_advance(const IndexType steps);
		void set_index(const IndexType new_index) {
			this->index_ = new_index;
		}
        bool set_coord(const int* r);

		/*
		get the global coordinate of the Site object.
		*/
		IndexType* coord() const;
        IndexType coord(const int direction) const;

		/*
		get the local (visible) coordinate of the Site object
		*/
        IndexType coord_local(const int direction) const;

		/*
		get the local mem coordinate of the Site object
		*/
		IndexType coord_local_mem(const int direction) const;

        template<int DIM>
        const Lattice<DIM>& get_lattice() const {return *(this->lattice_);}
        template<int DIM>
        Lattice<DIM>& get_lattice() const {return const_cast<Lattice>(*(this->lattice_));}
        const IndexType get_index() const {return this->index_;}

    private:
        const Lattice<DIM>* lattice_;
        /*
        this index_ should be the local_mem_index in the Lattice class convention.
        */
        IndexType index_;

    };

	template<int DIM>
	Site::Site(const Lattice<DIM>& lattice) 
		:
		lattice_(lattice),
		index_(0)
	{
		;
	}

	template<int DIM>
	Site::~Site() {

	}

	template<int DIM>
	Site& Site::first() {
		this->index_ = this->lattice_.get_local_visible_first();
		return *this;
	}

	template<int DIM>
	Site& Site::next() {
		this->index_ = this->lattice_.get_local_visible_next(this->index_);
		return *this;
	}

	template<int DIM>
	bool Site::test() const {
		return this->index_ < this->lattice_.get_local_visible_next_to_last();
	}
	
	template<int DIM>
	Site& Site::halo_first() {
		this->index_ = this->lattice_.get_local_halo_first();
		return *this;
	}

	template<int DIM>
	Site& Site::halo_next() {
		this->index_ = this->lattice_.get_local_halo_next(this->index_);
		return *this;
	}

	template<int DIM>
	bool Site::halo_test() const {
		return this->index_ < this->lattice_.get_local_halo_next_to_last();
	}
}

#endif