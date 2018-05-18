#ifndef PARASITE_LATTICE_H
#define PARASITE_LATTICE_H

#define _SCL_SECURE_NO_WARNINGS

#include <cassert>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <functional>
#include <vector>
#include <cstdint>
#include "./MPIWrapper2D.h"

namespace ParaSite{

    //IndexType and GridIndexType should be SIGNED integer types.
    //sizeof(GridIndexType) < sizeof(IndexType)
    typedef int32_t IndexType;
    typedef int32_t GridIndexType;
    typedef double DistanceType;

    /*
    helper function to transform Parallel2D grid_rank to grid_loc in this class.
    The grid_loc definition is independent of the Parallel2D implementation.
    If Parallel2D convention changes, we only need to change the following 2 functions.
    */
    inline GridIndexType* transform_gridrank_to_gridloc(const GridIndexType grid_rank[2], GridIndexType grid_loc[2]){
        grid_loc[0] = grid_rank[1];
        grid_loc[1] = grid_rank[0];
        return grid_loc;
    }

    inline GridIndexType* transform_gridloc_to_gridrank(const GridIndexType grid_loc[2], GridIndexType grid_rank[2]){
        grid_rank[0] = grid_loc[1];
        grid_rank[1] = grid_loc[0];
        return grid_rank;
    }

    /*
    class Lattice:
    stores the necessary information of a DIM-dimensional Cartesian lattice.
    The indexing convention is [x,y,z,...], with x changes first. (row-first convention)
    A parallalization is assumed in the last 2 dimensions.
    e.g., for 3D lattice, the y and z dimensions are scattered to a number of MPI processes.
    y will be identified with the MPI row, z will be identified with the MPI column.
    */
    template<int DIM>
    class Lattice{
    public:
        Lattice() = delete;
        Lattice(const IndexType* size, 
                const IndexType halo, 
                const GridIndexType* node_size = nullptr, 
                const GridIndexType* grid_loc = nullptr);
        Lattice(const IndexType size, 
                const IndexType halo, 
                const GridIndexType* node_size = nullptr,
                const GridIndexType* grid_loc = nullptr); // need test
		//Lattice(const IndexType* size,
		//	const IndexType halo,
		//	const Parallel2D& parallel);
        ~Lattice();
        void init(const IndexType* size, 
                const IndexType halo, 
                const GridIndexType* node_size = nullptr, 
                const GridIndexType* grid_loc = nullptr);
        
        /*
        make hash maps. to accelerate index query.
        normally, the table size should have the next-to-last site, to avoid out-of-range error.
        */
        void make_vis2mem_index_table();
        void delete_vis2mem_index_table();
       
        void make_mem2vis_index_table();
        void delete_mem2vis_index_table();

        void make_all_index_tables(){
            make_vis2mem_index_table();
            make_mem2vis_index_table();
        }
        void delete_all_index_tables(){
            delete_vis2mem_index_table();
            delete_mem2vis_index_table();
        }
        /*
        coordinate and index transformations.
        / global_index: the global index of a site, range from 0 to visible_global_sites_ - 1.
        / global_coord: the global coordinate of a site, DIM-dimensional vector.
        / local_mem_index: the index of a local (possibly non-visible) site in the full local lattice, 
          including halo sites. This gives the site location in the memory.
        / local_mem_coord: the coordinate of a local site in the full local lattice.
        / local_vis_coord: the coordinate of a local site in the visible local lattice.
        For efficiency, most of the functions will assume the input is valid.
        */
        IndexType global_coord2index(const IndexType* global_coord) const;
        IndexType* global_index2coord(const IndexType global_index, IndexType* global_coord) const;
		/*
        move the current position (global_index) with specific steps in some direction.
        periodic boundary condition is assumed.
        */
        IndexType* global_move(
            IndexType* global_coord, 
            const int steps, 
            const int direction) const; 
        IndexType global_move(
			const IndexType global_index, 
			const int steps, 
			const int direction) const;

        IndexType local_mem_coord2index(const IndexType* local_mem_coord) const;
		IndexType* local_mem_index2coord(const IndexType local_mem_index,
			IndexType* local_mem_coord) const;
        /*
        move the current local memory position (local_mem_index) with specific steps
        in some direction.
        locally periodic boundary condition is assumed.
        */
        IndexType* local_mem_move(
            IndexType* local_mem_coord,
            const int steps, 
            const int direction) const;
		IndexType local_mem_move(
            const IndexType local_mem_index,
			const int steps,
			const int direction) const;

        IndexType local_vis_coord2index(const IndexType* local_vis_coord) const;
        IndexType* local_vis_index2coord(const IndexType local_vis_index, 
            IndexType* local_vix_coord) const;

        /*
        transformation between local_mem_coord and local_vis_coord
        */
        IndexType* local_mem_coord_to_local_vis_coord(
            const IndexType* local_mem_coord,
            IndexType* local_vis_coord) const;
        IndexType* local_vis_coord_to_local_mem_coord(
            const IndexType* local_vis_coord,
            IndexType* local_mem_coord) const;

        /*
        transformation between local_vis_coord and global_coord.
        */
        IndexType* local_vis_coord_to_global_coord(
            const IndexType* local_vis_coord, 
            IndexType* global_coord) const;
        IndexType* global_coord_to_local_vis_coord(
            const IndexType* global_coord,
            IndexType* local_vis_coord) const;

        /*
        transform local_mem_coord to global_coord
        */
        IndexType* local_mem_coord_to_global_coord(
            const IndexType* local_mem_coord,
            IndexType* global_coord) const;

		/*
		transformation between three kinds of indices.
		/ when transforming from global_index to local_vis_index,
		the global_index is assumed to be in the local visible region.
		/ when transforming from local_mem_index to local_vis_index,
		the local_mem_index is assumed to be in the local visible region.
		*/
		IndexType local_vis_index_to_local_mem_index(
			const IndexType local_vis_index) const;
		IndexType local_mem_index_to_local_vis_index(
			const IndexType local_mem_index) const;

		IndexType global_index_to_local_vis_index(
			const IndexType global_index) const;
		IndexType local_vis_index_to_global_index(
			const IndexType local_vis_index) const;

        //check if the local_mem_index is a visible site.
        bool is_local_visible(const IndexType local_mem_index) const;
		bool is_local_visible(const IndexType* local_mem_coord) const;
        //check the global site belongs to which node
        GridIndexType* which_node(const IndexType global_index, GridIndexType* grid_node) const;
        //check if the global_index is a local (visible) site.
        bool is_local(const IndexType global_index) const;
        
        /*
		The following functions returns tht local_mem_index of:
		/ 1. the first local visible site
		/ 2. the next local visible site
		/ 3. the last local visible site + 1
		*/
		constexpr IndexType get_local_visible_first() const { return this->local_visible_site_first_; }
		IndexType get_local_visible_next(const IndexType local_mem_index) const;
		constexpr IndexType get_local_visible_next_to_last() const { return this->local_visible_site_next_to_last_; }
        
		/*
		The following functions returns the local_mem_index of 
		/ 1. the first local halo site. (0 if halo_==0)
		/ 2. the next local halo site
		/ 3. the last local halo site + 1. (0 if halo_==0)
		*/
		constexpr IndexType get_local_halo_first() const {return 0;}
		IndexType get_local_halo_next(const IndexType local_mem_index) const;
		constexpr IndexType get_local_halo_next_to_last() const;

        /*
        relate the halo site coordinate with the local visible coordinate being mapped (mapped_coord).
        All the coordinates are local_mem_coord.
        */
        IndexType get_mapped_coord(
           const IndexType halo_coord, const int direction) const;
        IndexType* get_mapped_coord(
            const IndexType* halo_coord, IndexType* mapped_coord) const;
        GridIndexType* get_mapped_grid_loc_from_halo_coord(
            const IndexType* halo_coord, GridIndexType grid_loc[2]) const;
        GridIndexType* get_mapped_grid_loc_from_rel_grid_loc(
            const GridIndexType rel_grid_loc[2], GridIndexType grid_loc[2]) const;
        GridIndexType* get_mapped_relative_grid_loc(
            const IndexType* halo_coord, GridIndexType rel_grid_loc[2]) const;    


        /*
        compute the lattice global distance between two global coordinates.
        assuming global periodic boundary.
        */
        DistanceType global_lat_distance(const IndexType* gc1, const IndexType* gc2) const;
        DistanceType global_lat_distance(const IndexType gidx1, const IndexType gidx2) const;
        

        //miscellaneous
        constexpr int get_dim() const {return DIM;}
        const IndexType* get_global_size() const {return this->global_size_;}
		const IndexType* get_local_size() const { return this->local_size_; }
		const IndexType* get_local_mem_size() const { return this->local_mem_size_; }
        constexpr IndexType get_halo() const {return this->halo_;}
        constexpr IndexType get_visible_global_sites() const {return this->visible_global_sites_;}
        constexpr IndexType get_visible_local_sites() const {return this->visible_local_sites_;}
        constexpr IndexType get_local_mem_sites() const {return this->local_mem_sites_;}

        constexpr GridIndexType get_total_nodes() const {return this->total_nodes_;}
        const GridIndexType* get_node_size() const {return this->node_size_;}
        const GridIndexType* get_grid_loc() const {return this->grid_loc_;}

    private:
        //global lattice size (excluding halo layers).
        IndexType global_size_[DIM];
        //local lattice size (excluding halo layers).
        IndexType local_size_[DIM];
        //local memory lattice size (including halo layers).
        IndexType local_mem_size_[DIM];
        //number of halo layers.
        IndexType halo_;

        //Total number of global lattice sites (excluding halo sites).
        IndexType visible_global_sites_;

        //Number of local sites (excluding halo sites)
        IndexType visible_local_sites_;
        //Number of gross local sites (including halo sites)
        IndexType local_mem_sites_;
        
        //global stride for global indexing
        IndexType global_stride_[DIM];
        //local stride for local indexing without halos
        IndexType local_stride_[DIM];
        //local stride for local indexing with halos
        IndexType local_mem_stride_[DIM];

        //The local_mem_index of the first visible local site
        IndexType local_visible_site_first_;
        //The local_mem_index of the next of the last visible local site
        IndexType local_visible_site_next_to_last_;

        //MPI related variables
        //The total number of MPI processes used to distribute the lattice
        GridIndexType total_nodes_ = 1;
        //The shape of the 2D MPI grid.
        //e.g., for DIM=3, 
        //grid_size_[0] partitions the y dimension; grid_size_[1] partitions the z direction.
        GridIndexType node_size_[2] = {1, 1};
        //the location of the local lattice in the global grid.
        /*
        grid_loc_[0] is the index in the grid_size_[0] direction.
        grid_loc_[1] is the index in the grid_size_[1] direction.
        Notice, if we use the Parallel2D class, then
        grid_loc_[0] = grid_rank_[1]
        grid_loc_[1] = grid_rank_[0]
        */
        GridIndexType grid_loc_[2] = {0, 0};  
        GridIndexType& loc_row_ = grid_loc_[0]; 
        GridIndexType& loc_col_ = grid_loc_[1];  

        /*
        index containers, to accelerate query.
        */ 
        bool activate_vis2mem_index_table_ = false;
        std::vector<IndexType> vis2mem_index_table_;
        
        bool activate_mem2vis_index_table_ = false;
        std::vector<IndexType> mem2vis_index_table_; /*unordered map is too slow*/
        //std::unordered_map<IndexType, IndexType> mem2vis_index_table_;

        /*TODO later: this does not add more functions; 
        just intended to improve efficiency.*/
        bool activate_vis2glb_index_table_ = false;
        std::vector<IndexType> vis2glb_index_table_;
    };

    template<int DIM>
    Lattice<DIM>::Lattice(const IndexType* size, 
                    const IndexType halo, 
                    const GridIndexType* node_size,
                    const GridIndexType* grid_loc)
    {
        this->init(size, halo, node_size, grid_loc);
    }

    template<int DIM>
    Lattice<DIM>::Lattice(const IndexType size, 
                    const IndexType halo, 
                    const GridIndexType* node_size,
                    const GridIndexType* grid_loc){
        IndexType arr_size[DIM];
        std::fill_n(arr_size, DIM, size);
        this->init(arr_size, halo, node_size, grid_loc);
    }

    template<int DIM>
    void Lattice<DIM>::init(const IndexType* size, 
                    const IndexType halo, 
                    const GridIndexType* node_size,
                    const GridIndexType* grid_loc) {
        halo_ = halo;
        std::copy_n(size, DIM, global_size_);
        if((node_size != nullptr) && (grid_loc != nullptr)){
            node_size_[0] = node_size[0];
            node_size_[1] = node_size[1];
            grid_loc_[0] = grid_loc[0];
            grid_loc_[1] = grid_loc[1];
        }
        static_assert(DIM > 1, "Dimension must be larger than 2");  
        total_nodes_ = node_size_[0] * node_size_[1];
        //first check the size of the last two dimensions are integer multiples of node sizes.
        assert(!(global_size_[DIM-1] % node_size_[1]));
        assert(!(global_size_[DIM-2] % node_size_[0]));
        //compute local_size_
        local_size_[DIM-1] = global_size_[DIM-1] / node_size_[1];
        local_size_[DIM-2] = global_size_[DIM-2] / node_size_[0];
        std::copy_n(global_size_ , DIM -2, local_size_);
        //the halo layers should be smaller than the visible size
        assert(
            std::all_of(local_size_, local_size_ + DIM, 
                    [this](IndexType x){ return this->halo_ <= x / 2;})
            );
        //compute local_mem_size_
        std::transform(local_size_, local_size_ + DIM, 
                        local_mem_size_, 
                        [this](IndexType x){ return x + 2*this->halo_;} );
        //compute visible_global_sites_
        visible_global_sites_ = std::accumulate(global_size_, global_size_ + DIM, 
                                                1, std::multiplies<IndexType>());
        visible_local_sites_ = std::accumulate(local_size_, local_size_ + DIM,
                                                1, std::multiplies<IndexType>());
        local_mem_sites_ = std::accumulate(local_mem_size_, local_mem_size_ + DIM,
                                                1, std::multiplies<IndexType>());
        //compute global_stride_
        global_stride_[0] = 1;
        std::partial_sum(global_size_, global_size_ + DIM - 1, 
                        global_stride_ + 1, std::multiplies<IndexType>());
        //compute local_stride_
        local_stride_[0] = 1;
        std::partial_sum(local_size_, local_size_ + DIM - 1, 
                        local_stride_ + 1, std::multiplies<IndexType>());
        //compute local_mem_stride_
        local_mem_stride_[0] = 1;
        std::partial_sum(local_mem_size_, local_mem_size_ + DIM - 1,
                        local_mem_stride_ + 1, std::multiplies<IndexType>());

        //compute local_visible_site_first_
        IndexType first[DIM];
        std::fill_n(first, DIM, this->halo_);
        local_visible_site_first_ = local_mem_coord2index(first);
        //compute local_visible_site_next_to_last_
        IndexType last[DIM];
        std::transform(local_mem_size_, local_mem_size_ + DIM,
                        last,
                        [this](IndexType x){ return x - 1 - this->halo_; });
        local_visible_site_next_to_last_ = local_mem_coord2index(last) + 1;
        return;
    }

    template<int DIM>
    Lattice<DIM>::~Lattice(){
        ;
    }

    template<int DIM>
    void Lattice<DIM>::make_vis2mem_index_table() {
        vis2mem_index_table_ = std::vector<IndexType>(this->visible_local_sites_ + 1);
        auto lv_idx = this->get_local_visible_first();
        IndexType count = 0;
        while(lv_idx < this->get_local_visible_next_to_last()){
            vis2mem_index_table_[count] = lv_idx;
            lv_idx = this->get_local_visible_next(lv_idx);
            count++;
        }
        *(vis2mem_index_table_.rbegin()) = this->get_local_visible_next_to_last();
        activate_vis2mem_index_table_ = true;
    }

    template<int DIM>
    void Lattice<DIM>::delete_vis2mem_index_table(){
        vis2mem_index_table_ = std::vector<IndexType>();
        activate_vis2mem_index_table_ = false;
    }

    template<int DIM>
    void Lattice<DIM>::make_mem2vis_index_table(){
        mem2vis_index_table_ = std::vector<IndexType>(this->local_mem_sites_, 0);
        auto lv_idx = this->get_local_visible_first();
        IndexType count = 0;
        while(lv_idx < this->get_local_visible_next_to_last()){
            mem2vis_index_table_[lv_idx] = count;
            lv_idx = this->get_local_visible_next(lv_idx);
            count++;
        }
        //the next-to-last element should be treated separately.
        //since get_local_visible_next(idx) is not well-defined beyond the last point. 
        mem2vis_index_table_[this->get_local_visible_next_to_last()] = count;
        activate_mem2vis_index_table_ = true;
    }

    template<int DIM>
    void Lattice<DIM>::delete_mem2vis_index_table(){
        activate_mem2vis_index_table_ = false;
        mem2vis_index_table_ = std::vector<IndexType>();
    }

    template<int DIM>
    IndexType Lattice<DIM>::global_coord2index(const IndexType* global_coord) const{
        return std::inner_product(this->global_stride_, this->global_stride_ + DIM,
                                    global_coord, 0);
    }

    template<int DIM>
    IndexType* Lattice<DIM>::global_index2coord(const IndexType global_index, 
                                            IndexType* global_coord) const{
        auto idx = global_index;
        for(int i = DIM-1; i >= 0; --i){
            global_coord[i] = idx / this->global_stride_[i];
            idx = idx % this->global_stride_[i];
        }
        return global_coord;
    }

    template<int DIM>
    IndexType* Lattice<DIM>::global_move(
        IndexType* global_coord,
		const int steps,
		const int direction) const {
        const auto& i = direction;
        global_coord[i] = 
            ( (global_coord[i] + steps) % this->global_size_[i] + this->global_size_[i] ) 
            % this->global_size_[i];
        return global_coord;
    }

	template<int DIM>
	IndexType Lattice<DIM>::global_move(
		const IndexType global_index,
		const int steps,
		const int direction) const {
        IndexType global_coord[DIM];
        global_index2coord(global_index, global_coord);
        this->global_move(global_coord, steps, direction);
        return global_coord2index(global_coord);
	}

    template<int DIM>
    IndexType Lattice<DIM>::local_mem_coord2index(const IndexType* local_mem_coord) const{
        return std::inner_product(local_mem_coord, local_mem_coord + DIM,
                                    this->local_mem_stride_, 0);
    }
    
	template<int DIM>
    IndexType* Lattice<DIM>::local_mem_index2coord(const IndexType local_mem_index,
                                    IndexType* local_mem_coord) const{
		auto idx = local_mem_index;
		for (auto i = DIM - 1; i >= 0; --i) {
			local_mem_coord[i] = idx / this->local_mem_stride_[i];
			idx = idx % this->local_mem_stride_[i];
		}
		return local_mem_coord;
    }

    template<int DIM>
    IndexType* Lattice<DIM>::local_mem_move(
        IndexType* local_mem_coord,
        const int steps,
        const int direction) const {
        const auto& i = direction;
        local_mem_coord[i] = ( (local_mem_coord[i]+steps) % this->local_mem_size_[i] + this->local_mem_size_[i] )
                    % this->local_mem_size_[i];
        return local_mem_coord;
    }

    template<int DIM>
    IndexType Lattice<DIM>::local_mem_move(const IndexType local_mem_index,
			const int steps,
			const int direction) const{
        IndexType lm_coord[DIM];
        local_mem_index2coord(local_mem_index, lm_coord);
        this->local_mem_move(lm_coord, steps, direction);
        return local_mem_coord2index(lm_coord);
    }

    template<int DIM>
    IndexType Lattice<DIM>::local_vis_coord2index(const IndexType* local_vis_coord) const {
        return std::inner_product(local_vis_coord, local_vis_coord + DIM, 
                                    this->local_stride_, 0);
    }

    template<int DIM>
    IndexType* Lattice<DIM>::local_vis_index2coord(const IndexType local_vis_index, 
                                                IndexType* local_vis_coord) const {
        auto idx = local_vis_index;
        for(auto i = DIM-1; i >= 0; --i){
            local_vis_coord[i] = idx / this->local_stride_[i];
            idx = idx % this->local_stride_[i];
        }
        return local_vis_coord;
    }

    template<int DIM>
    bool Lattice<DIM>::is_local_visible(const IndexType local_mem_index) const{
		IndexType local_mem_coord[DIM];
		return is_local_visible(
			local_mem_index2coord(
				local_mem_index,
				local_mem_coord));
    }

	template<int DIM>
	bool Lattice<DIM>::is_local_visible(const IndexType* local_mem_coord) const {
		for (auto i = 0; i < DIM; ++i) {
			if ((local_mem_coord[i] >= halo_) && (local_mem_coord[i] < local_mem_size_[i] - halo_)) continue;
			else return false;
		}
		return true;
	}

    template<int DIM>
    GridIndexType* Lattice<DIM>::which_node(const IndexType global_index, GridIndexType* grid_node) const{
        IndexType global_coord[DIM];
        global_index2coord(global_index, global_coord);
        auto& g1_coord = global_coord[DIM - 2];
        auto& g2_coord = global_coord[DIM - 1];
        grid_node[0] = g1_coord / this->local_size_[DIM - 2];
        grid_node[1] = g2_coord / this->local_size_[DIM - 1];
        return grid_node;
    }

    template<int DIM>
    bool Lattice<DIM>::is_local(const IndexType global_index) const{
        GridIndexType grid_node[2];
        this->which_node(global_index, grid_node);
        if(grid_node[0] == this->grid_loc_[0] && grid_node[1] == this->grid_loc_[1])
            return true;
        else
            return false;
    }

    template<int DIM>
    IndexType* Lattice<DIM>::local_mem_coord_to_local_vis_coord(
            const IndexType* local_mem_coord,
            IndexType* local_vis_coord) const {
        std::transform(local_mem_coord, local_mem_coord + DIM,
                        local_vis_coord,
                        [this](IndexType x){ return x - this->halo_; });
        return local_vis_coord;
    }

    template<int DIM>
    IndexType* Lattice<DIM>::local_vis_coord_to_local_mem_coord(
        const IndexType* local_vis_coord,
        IndexType* local_mem_coord) const {
        std::transform(local_vis_coord, local_vis_coord + DIM,
                        local_mem_coord,
                        [this](IndexType x){ return x + this->halo_; });
        return local_mem_coord;
    }

    template<int DIM>
    IndexType* Lattice<DIM>::local_vis_coord_to_global_coord(
            const IndexType* local_vis_coord, 
            IndexType* global_coord) const {
        global_coord[DIM - 2] = local_vis_coord[DIM - 2] + this->grid_loc_[0] * this->local_size_[DIM - 2];
        global_coord[DIM - 1] = local_vis_coord[DIM - 1] + this->grid_loc_[1] * this->local_size_[DIM - 1];
        std::copy_n(local_vis_coord, DIM - 2, global_coord);
        return global_coord;
    }
    
    template<int DIM>
    IndexType* Lattice<DIM>::global_coord_to_local_vis_coord(
        const IndexType* global_coord,
        IndexType* local_vis_coord) const{
        local_vis_coord[DIM-2] = global_coord[DIM-2] % this->local_size_[DIM-2];
        local_vis_coord[DIM-1] = global_coord[DIM-1] % this->local_size_[DIM-1];
        std::copy_n(global_coord, DIM-2, local_vis_coord);
        return local_vis_coord;
    }

    template<int DIM>
    IndexType* Lattice<DIM>::local_mem_coord_to_global_coord(
        const IndexType* local_mem_coord,
        IndexType* global_coord) const {
        IndexType lm_first_vis[DIM];
        local_mem_index2coord(this->local_visible_site_first_, lm_first_vis);
        IndexType first_vis_coord[DIM], first_vis_glb_coord[DIM];
        std::fill_n(first_vis_coord, DIM, 0);
        this->local_vis_coord_to_global_coord(first_vis_coord, first_vis_glb_coord);
        IndexType offset[DIM];
        std::transform(local_mem_coord, local_mem_coord + DIM, 
                        lm_first_vis, 
                        offset, std::minus<IndexType>());
        std::transform(first_vis_glb_coord, first_vis_glb_coord + DIM,
                        offset, 
                        global_coord, std::plus<IndexType>());
        std::transform(global_coord, global_coord + DIM, 
                        this->global_size_, 
                        global_coord, 
                        [](IndexType x, IndexType y)
                        { return (x%y + y)%y; });
        return global_coord;
    }

    template<int DIM>
	IndexType Lattice<DIM>::local_vis_index_to_local_mem_index(
		const IndexType local_vis_index) const {
        if(activate_vis2mem_index_table_ == true){
            return vis2mem_index_table_.at(local_vis_index);
        }
        else{
            IndexType lv_coord[DIM], lm_coord[DIM];
            local_vis_index2coord(local_vis_index, lv_coord);
            local_vis_coord_to_local_mem_coord(lv_coord, lm_coord);
            return local_mem_coord2index(lm_coord);
        }
	}

	template<int DIM>
	IndexType Lattice<DIM>::local_mem_index_to_local_vis_index(
        const IndexType local_mem_index) const {
        if(activate_mem2vis_index_table_){
            return mem2vis_index_table_.at(local_mem_index);
        }
        else{
            IndexType lv_coord[DIM], lm_coord[DIM];
            local_mem_index2coord(local_mem_index, lm_coord);
            local_mem_coord_to_local_vis_coord(lm_coord, lv_coord);
            return local_vis_coord2index(lv_coord);
        }
	}

	template<int DIM>
	IndexType Lattice<DIM>::global_index_to_local_vis_index(
		const IndexType global_index) const {
		IndexType global_coord[DIM], local_vis_coord[DIM];
		global_index2coord(global_index, global_coord);
		global_coord_to_local_vis_coord(global_coord, local_vis_coord);
		return local_vis_coord2index(local_vis_coord);
	}

	template<int DIM>
	IndexType Lattice<DIM>::local_vis_index_to_global_index(
		const IndexType local_vis_index) const {
		IndexType global_coord[DIM], lv_coord[DIM];
		local_vis_index2coord(local_vis_index, lv_coord);
		local_vis_coord_to_global_coord(lv_coord, global_coord);
		return global_coord2index(global_coord);
	}

    template<int DIM>
    IndexType Lattice<DIM>::get_local_visible_next(const IndexType local_mem_index) const {
		auto local_vis_index_new = local_mem_index_to_local_vis_index(local_mem_index) + 1;
		return local_vis_index_to_local_mem_index(local_vis_index_new);
    }

    template<int DIM>
    IndexType Lattice<DIM>::get_local_halo_next(const IndexType local_mem_index) const{
        IndexType local_mem_coord[DIM];
        local_mem_index2coord(local_mem_index + 1, local_mem_coord);
        if(this->is_local_visible(local_mem_coord)){
            local_mem_coord[0] = this->local_mem_size_[0] - this->halo_;
        }
		return local_mem_coord2index(local_mem_coord);
    }

    template<int DIM>
    constexpr IndexType Lattice<DIM>::get_local_halo_next_to_last() const {
        return (this->halo_ <= 0) ? 0 : this->local_mem_sites_;
    }

    template<int DIM>
    IndexType Lattice<DIM>::get_mapped_coord(
           const IndexType halo_coord, const int direction) const{
        if(halo_coord < this->halo_) 
            return this->local_mem_size_[direction] - 2*this->halo_ + halo_coord;
        if(halo_coord >= this->local_mem_size_[direction] - this->halo_)
            return 2*this->halo_ - this->local_mem_size_[direction] + halo_coord;
        return halo_coord;
    }

    template<int DIM>
    IndexType* Lattice<DIM>::get_mapped_coord(
        const IndexType* halo_coord, IndexType* mapped_coord) const{
        for(auto i = 0; i < DIM; ++i){
            mapped_coord[i] = get_mapped_coord(halo_coord[i], i);
        }
        return mapped_coord;
    }

    template<int DIM>
    GridIndexType* Lattice<DIM>::get_mapped_grid_loc_from_halo_coord(
        const IndexType* halo_coord, GridIndexType grid_loc[2]) const{
        GridIndexType rel_grid_loc[2];
        this->get_mapped_relative_grid_loc(halo_coord, rel_grid_loc);
        return get_mapped_grid_loc_from_rel_grid_loc(rel_grid_loc, grid_loc);
    }

    template<int DIM>
    GridIndexType* Lattice<DIM>::get_mapped_grid_loc_from_rel_grid_loc(
        const GridIndexType rel_grid_loc[2], GridIndexType grid_loc[2]) const{
        grid_loc[0] = ( (this->grid_loc_[0] + rel_grid_loc[0]) % this->node_size_[0]
                        + this->node_size_[0] ) % this->node_size_[0];
        grid_loc[1] = ( (this->grid_loc_[1] + rel_grid_loc[1]) % this->node_size_[1]
                        + this->node_size_[1] ) % this->node_size_[1];
        return grid_loc;
    }

    template<int DIM>
    GridIndexType* Lattice<DIM>::get_mapped_relative_grid_loc(
        const IndexType* halo_coord, GridIndexType rel_grid_loc[2]) const{
        if(halo_coord[DIM-2] < this->halo_ && this->node_size_[0] > 1){
            rel_grid_loc[0] = -1;
        } else if(halo_coord[DIM-2] >= this->local_mem_size_[DIM-2] - this->halo_ 
                    && this->node_size_[0] > 1){
            rel_grid_loc[0] = 1;
        } else{
            rel_grid_loc[0] = 0;
        }

        if(halo_coord[DIM-1] < this->halo_ && this->node_size_[1] > 1){
            rel_grid_loc[1] = -1;
        } else if(halo_coord[DIM-1] >= this->local_mem_size_[DIM-1] - this->halo_
                    && this->node_size_[1] > 1){
            rel_grid_loc[1] = 1;
        } else{
            rel_grid_loc[1] = 0;
        }
        return rel_grid_loc;
    }

    template<int DIM>
    DistanceType Lattice<DIM>::global_lat_distance(
        const IndexType* gc1, const IndexType* gc2) const {
        IndexType np_dspl[DIM]; //non-periodic displacement
        IndexType bd_dspl[DIM]; //displacement crossing border
        IndexType dspl[DIM];
        std::transform(gc1, gc1 + DIM, gc2, 
                    np_dspl,
                    [](auto x1, auto x2) {return std::abs(x2-x1); });
        std::transform(np_dspl, np_dspl + DIM, this->global_size_, 
                    bd_dspl,
                    [](auto x1, auto x2) {return std::abs(x2-x1); });
        std::transform(np_dspl, np_dspl + DIM, bd_dspl, 
                    dspl, 
                    [](auto x1, auto x2) {return x1 < x2 ? x1 : x2; });
        return sqrt( std::inner_product(dspl, dspl + DIM, dspl, 0) );
    }

    template<int DIM>
    DistanceType Lattice<DIM>::global_lat_distance(
        const IndexType gidx1, const IndexType gidx2) const {
        IndexType gc1[DIM], gc2[DIM];
        this->global_index2coord(gidx1, gc1);
        this->global_index2coord(gidx2, gc2);
        return global_lat_distance(gc1, gc2);
    }



    // Non-member functions
    template<int DIM>
    IndexType* local_vis_coord_to_global_coord(
            const IndexType* local_vis_coord, 
            IndexType* global_coord,
            const GridIndexType* grid_loc,
            const IndexType* local_size) {
        global_coord[DIM - 2] = local_vis_coord[DIM - 2] + grid_loc[0] * local_size[DIM - 2];
        global_coord[DIM - 1] = local_vis_coord[DIM - 1] + grid_loc[1] * local_size[DIM - 1];
        std::copy_n(local_vis_coord, DIM - 2, global_coord);
        return global_coord;
    }

}

#endif