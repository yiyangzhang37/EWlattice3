#ifndef PARASITE_SITE_H
#define PARASITE_SITE_H

#include "ParaSite_Lattice.h"

namespace ParaSite{
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
		All the move operations should not have any relation with global index.
		Since the Site::index_ is not defined beyond the local_mem_region.
		Note: it is allowed to move to the halo layers.
		*/
		/*
		especially in for(x.first(); x.test(); x.next()) loops,
		the const methods should be called.
		*/
		Site move(const int direction, const int steps) const;
        Site operator+(const int direction) const;
        Site operator-(const int direction) const;

		/*
		emplace move operations.
		operator> corresponds to operator+.
		operator< corresponds to operator-.
		*/
		Site& emplace_move(const int direction, const int steps);
        Site& operator>(const int direction);
        Site& operator<(const int direction);

		void set_index(const IndexType new_index) {
			this->index_ = new_index;
		}
		/*
		set the Site object to the global coordinate.
		return true if this coordinate is in the local visible region.
		otherwise return false and index_ unchanged.
		*/
        bool set_coord(const IndexType* global_coord);

		/*
		get the global coordinate of the Site object.
		*/
		IndexType* coord(IndexType* global_coord) const;
        IndexType coord(const int direction) const;

		/*
		get the local (visible) coordinate of the Site object
		*/
		IndexType* coord_local_vis(IndexType* local_vis_coord) const;
        IndexType coord_local_vis(const int direction) const;

		/*
		get the local mem coordinate of the Site object
		*/
		IndexType* coord_local_mem(IndexType* local_mem_coord) const;
		IndexType coord_local_mem(const int direction) const;

        const Lattice<DIM>& get_lattice() const {return *(this->lattice_);}
        //Lattice<DIM>& get_lattice() const {return const_cast<Lattice<DIM>>(*(this->lattice_));}
        constexpr IndexType get_index() const {return this->index_;}

    private:
        const Lattice<DIM>* lattice_;
        /*
        this index_ should be the local_mem_index in the Lattice class convention.
        */
        IndexType index_;

    };

	template<int DIM>
	Site<DIM>::Site(const Lattice<DIM>& lattice) 
		:
		lattice_(&lattice),
		index_(0)
	{
		;
	}

	template<int DIM>
	Site<DIM>::~Site() {
		;
	}

	template<int DIM>
	Site<DIM>& Site<DIM>::emplace_move(const int direction, const int steps){
		this->index_ = this->lattice_->local_mem_move(this->index_, steps, direction);
		return *this;
	}

	template<int DIM>
	Site<DIM> Site<DIM>::move(const int direction, const int steps) const {
		Site<DIM> x = *this;
		x.set_index(this->lattice_->local_mem_move(this->index_, steps, direction));
		return x;
	}

	template<int DIM>
	Site<DIM> Site<DIM>::operator+(const int direction) const {
		return this->move(direction, 1);
	}

	template<int DIM>
	Site<DIM>& Site<DIM>::operator>(const int direction) {
		return this->emplace_move(direction, 1);
	}

	template<int DIM>
	Site<DIM> Site<DIM>::operator-(const int direction) const {
		return this->move(direction, -1);
	}

	template<int DIM>
	Site<DIM>& Site<DIM>::operator<(const int direction) {
		return this->emplace_move(direction, -1);
	}

	template<int DIM>
	bool Site<DIM>::set_coord(const IndexType* global_coord) {
		auto global_index = this->lattice_->global_coord2index(global_coord);
		auto is_local = this->lattice_->is_local(global_index);
		if (is_local) {
			IndexType lv_coord[DIM], lm_coord[DIM];
			this->lattice_->global_coord_to_local_vis_coord(global_coord, lv_coord);
			this->lattice_->local_vis_coord_to_local_mem_coord(lv_coord, lm_coord);
			this->index_ = this->lattice_->local_mem_coord2index(lm_coord);
			return true;
		}
		else {
			return false;
		}
	}

	template<int DIM>
	IndexType* Site<DIM>::coord(IndexType* global_coord) const{
		IndexType lm_coord[DIM];
		IndexType lv_coord[DIM];
		this->lattice_->local_mem_index2coord(this->index_, lm_coord);
		this->lattice_->local_mem_coord_to_local_vis_coord(lm_coord, lv_coord);
		this->lattice_->local_vis_coord_to_global_coord(lv_coord, global_coord);
		return global_coord;
	}

	template<int DIM>
	IndexType Site<DIM>::coord(const int direction) const{
		IndexType global_coord[DIM];
		coord(global_coord);
		return global_coord[direction];
	}

	template<int DIM>
	IndexType* Site<DIM>::coord_local_vis(IndexType* local_vis_coord) const{
		IndexType lm_coord[DIM];
		this->lattice_->local_mem_index2coord(this->index_, lm_coord);
		this->lattice_->local_mem_coord_to_local_vis_coord(lm_coord, local_vis_coord);
		return local_vis_coord;
	}
	
	template<int DIM>
	IndexType Site<DIM>::coord_local_vis(const int direction) const {
		IndexType lv_coord[DIM];
		coord_local_vis(lv_coord);
		return lv_coord[direction];
	}

	template<int DIM>
	IndexType* Site<DIM>::coord_local_mem(IndexType* local_mem_coord) const{
		this->lattice_->local_mem_index2coord(this->index_, local_mem_coord);
		return local_mem_coord;
	}

	template<int DIM>
	IndexType Site<DIM>::coord_local_mem(const int direction) const {
		IndexType lm_coord[DIM];
		coord_local_mem(lm_coord);
		return lm_coord[direction];
	}


	template<int DIM>
	Site<DIM>& Site<DIM>::first() {
		this->index_ = this->lattice_->get_local_visible_first();
		return *this;
	}

	template<int DIM>
	Site<DIM>& Site<DIM>::next() {
		this->index_ = this->lattice_->get_local_visible_next(this->index_);
		return *this;
	}

	template<int DIM>
	bool Site<DIM>::test() const {
		return this->index_ < this->lattice_->get_local_visible_next_to_last();
	}
	
	template<int DIM>
	Site<DIM>& Site<DIM>::halo_first() {
		this->index_ = this->lattice_->get_local_halo_first();
		return *this;
	}

	template<int DIM>
	Site<DIM>& Site<DIM>::halo_next() {
		this->index_ = this->lattice_->get_local_halo_next(this->index_);
		return *this;
	}

	template<int DIM>
	bool Site<DIM>::halo_test() const {
		return this->index_ < this->lattice_->get_local_halo_next_to_last();
	}
}

#endif