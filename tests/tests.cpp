#include "tests.h"
#include <random>
#include <functional>
#include <cmath>
#include <set>
#include <array>
#include <typeinfo>
//#include <cxxabi.h>

int Test_Build_Parallel2D(const int n_rows, const int n_cols) {
	Parallel2D parallel;
	auto total_size = n_rows * n_cols;
	int status = 1;
	//check world_size
	status = compare_results(total_size, parallel.get_world_size(), "World size") && status;
	//check world rank
	if (parallel.get_world_size() == 1) {
		status = compare_results(0, parallel.get_world_rank(), "World rank (Single process)") && status;
	}
	else {
		std::cout << "Process world rank: " << parallel.get_world_rank() << std::endl;
	}
	//check root ID
	status = compare_results(0, parallel.get_root(), "Root ID") && status;
	//check is_root
	if (parallel.get_world_rank() == parallel.get_root()) {
		status = compare_results(true, parallel.is_root(), "IsRoot") && status;
	}
	else {
		status = compare_results(false , parallel.is_root(), "IsRoot") && status;
	}
	//check world group size
	int world_group_size = 0;
	MPI_Group_size(parallel.get_world_group(), &world_group_size);
	status = compare_results(parallel.get_world_size(),
		world_group_size,
		"World group size") && status;
	//check world group rank
	int world_group_rank = 0;
	MPI_Group_rank(parallel.get_world_group(), &world_group_rank);
	status = compare_results(parallel.get_world_rank(),
		world_group_rank,
		"World group rank") && status;

	return status;
}

int Test_Build_2DGrid(const int n_rows, const int n_cols) {
	Parallel2D parallel;
	parallel.InitializeGrid(n_rows, n_cols);
	auto total_size = n_rows * n_cols;
	int status = 1;
	//check grid size
	status = compare_results(n_rows, parallel.get_grid_size()[0], "Grid rows") && status;
	status = compare_results(n_cols, parallel.get_grid_size()[1], "Grid columns") && status;
	auto row_rank = parallel.get_world_rank() % n_cols;
	auto col_rank = parallel.get_world_rank() / n_cols;
	//check grid ranks
	status = compare_results(row_rank, parallel.get_grid_rank()[0], "Row rank") && status;
	status = compare_results(col_rank, parallel.get_grid_rank()[1], "Column rank") && status;
	auto rebuild_world_rank = parallel.get_grid_rank()[0] * parallel.get_stride()[1] 
							+ parallel.get_grid_rank()[1] * parallel.get_stride()[0];
	status = compare_results(rebuild_world_rank, parallel.get_world_rank(),"Rebuild world rank") && status; 
	return status;
}

int Test_Broadcast(const int n_rows, const int n_cols){
	Parallel2D parallel;
	parallel.InitializeGrid(n_rows, n_cols);
	int status = 1;

	auto bcast_number = 2 * parallel.get_world_rank();
	int bcast_array[] = {bcast_number, 3 * bcast_number, 7 * bcast_number};
	//world broadcast
	auto b_from = parallel.get_world_size() - 1;
	parallel.Broadcast(bcast_number, b_from);
	status = compare_results(2 * b_from,
							bcast_number, "World broadcast") && status;
	parallel.Broadcast(bcast_array, 3, b_from);
	int theo_vals[] = {2*b_from, 6*b_from, 14*b_from};
	status = compare_results(theo_vals, bcast_array, 3,
							"World broadcast array") && status;
	
	//row broadcast
	b_from = parallel.get_grid_size()[1] - 1;
	bcast_number = 3 * parallel.get_grid_rank()[0];
	int bcast_row_array[] = {bcast_number, 3*bcast_number, 7*bcast_number};
	parallel.Broadcast_row(bcast_number, b_from);
	status = compare_results(3 * b_from, bcast_number, "Row broadcast") && status;
	parallel.Broadcast_row(bcast_row_array, 3, b_from);
	int theo_row_vals[] = {3*b_from, 9*b_from, 21*b_from};
	status = compare_results(theo_row_vals, bcast_row_array, 3, "Row broadcast array") && status;

	//column broadcast
	b_from = parallel.get_grid_size()[0] - 1;
	bcast_number = 17 * parallel.get_grid_rank()[1];
	int bcast_col_array[] = {bcast_number, 3*bcast_number, 7*bcast_number};
	parallel.Broadcast_col(bcast_number, b_from);
	status = compare_results(17 * b_from, bcast_number, "Column broadcast") && status;
	parallel.Broadcast_col(bcast_col_array, 3, b_from);
	int theo_col_vals[] = {17*b_from, 51*b_from, 119*b_from};
	status = compare_results(theo_col_vals, bcast_col_array, 3, "Column broadcast array");
	return status;
}

int Test_Sum(const int n_rows, const int n_cols){
	Parallel2D parallel;
	parallel.InitializeGrid(n_rows, n_cols);
	int status = 1;

	//world sum
	auto number = parallel.get_world_rank() * parallel.get_world_rank();
	parallel.Sum(number);
	const int& n = parallel.get_world_size();
	status = compare_results(n*(n-1)*(2*n-1)/6, number, "World Sum") && status;
	const auto& x = parallel.get_world_rank() + 1;
	int array[] = {x, x * x * x };
	parallel.Sum(array, 2);
	int theo_vals[] = {n*(n+1)/2, n*n*(n+1)*(n+1)/4};
	status = compare_results(theo_vals, array, 2, "World sum arrays") && status;
	
	//row sum
	number = parallel.get_grid_rank()[0] * parallel.get_grid_rank()[0];
	const auto& xr = parallel.get_grid_rank()[0] + 1;
	int row_array[] = {xr, xr*xr*xr};
	parallel.Sum_row(number);
	const int& nc = parallel.get_grid_size()[1];
	status = compare_results(nc*(nc-1)*(2*nc-1)/6, number, "Row sum") && status;
	parallel.Sum_row(row_array, 2);
	int theo_row_vals[] = {nc*(nc+1)/2, nc*nc*(nc+1)*(nc+1)/4};
	status = compare_results(theo_row_vals, row_array, 2, "Row sum arrays") && status;

	//column sum
	number = parallel.get_grid_rank()[1] * parallel.get_grid_rank()[1];
	const auto& xc = parallel.get_grid_rank()[1] + 1;
	int col_array[] = {xc, xc*xc*xc };
	parallel.Sum_col(number);
	const int& nr = parallel.get_grid_size()[0];
	status = compare_results(nr*(nr-1)*(2*nr-1)/6, number, "Column sum") && status;
	parallel.Sum_col(col_array, 2);
	int theo_col_vals[] = {nr*(nr+1)/2, nr*nr*(nr+1)*(nr+1)/4};
	status = compare_results(theo_col_vals, col_array, 2, "Column sum arrays") && status;

	return status;
}

int Test_SendRecv(const int n_rows, const int n_cols){
	Parallel2D parallel;
	parallel.InitializeGrid(n_rows, n_cols);
	int status = 1;

	//world send-receive
	int send_num = (parallel.get_world_rank()+1) * 7;
	int recv_num = 0;
	int exchange = parallel.get_world_size() - 1 - parallel.get_world_rank();
	parallel.Send(send_num, exchange);
	parallel.Receive(recv_num, exchange);
	const int& n = (exchange+1)*7;
	status = compare_results( n, recv_num, "World send-recv");

	int send_array[] = {send_num*send_num, send_num*send_num*send_num};
	int recv_array[] = {0, 0};
	parallel.Send(send_array, 2, exchange);
	parallel.Receive(recv_array, 2, exchange);
	int theo_vals[] = {n*n, n*n*n};
	status = compare_results(theo_vals, recv_array, 2, "World send-recv array") && status;

	//row send-receive
	send_num = (parallel.get_grid_rank()[0] + 1) * 11;
	recv_num = 0;
	exchange = parallel.get_grid_size()[1] - 1 - parallel.get_grid_rank()[0];
	parallel.Send_row(send_num, exchange);
	parallel.Receive_row(recv_num, exchange);
	const int& nr = (exchange+1)*11;
	status = compare_results(nr, recv_num, "Row send-recv") && status;

	int send_row_array[] = {send_num*send_num, send_num*send_num*send_num};
	int recv_row_array[] = {0, 0};
	parallel.Send_row(send_row_array, 2, exchange);
	parallel.Receive_row(recv_row_array, 2, exchange);
	int theo_row_vals[] = {nr*nr, nr*nr*nr};
	status = compare_results(theo_row_vals, recv_row_array, 2, "Row send-recv array") && status;

	//column send-receive
	send_num = (parallel.get_grid_rank()[1] + 1) * 19;
	recv_num = 0;
	exchange = parallel.get_grid_size()[0] - 1 - parallel.get_grid_rank()[1];
	parallel.Send_col(send_num, exchange);
	parallel.Receive_col(recv_num, exchange);
	const int& nc = (exchange+1)*19;
	status = compare_results(nc, recv_num, "Column send-recv") && status;

	int send_col_array[] = {send_num*send_num, send_num*send_num*send_num};
	int recv_col_array[] = {0, 0};
	parallel.Send_col(send_col_array, 2, exchange);
	parallel.Receive_col(recv_col_array, 2, exchange);
	int theo_col_vals[] = {nc*nc, nc*nc*nc};
	status = compare_results(theo_col_vals, recv_col_array, 2, "Column send-recv array") && status;

	return status;
}

int Test_SendUp(const int n_rows, const int n_cols){
	Parallel2D parallel;
	parallel.InitializeGrid(n_rows, n_cols);
	int status = 1;

	//sendup-row
	const int& nr = parallel.get_grid_rank()[0];
	const int& sc = parallel.get_grid_size()[1];
	int sendup_row[] = { nr + 5, nr*nr, (nr-3)*(nr+4) };
	int recvup_row[] = {-1, -3, -6};
	int send_row[] = { nr + 5, nr*nr, (nr-3)*(nr+4) };
	int recv_row[] = {0, 0, 0};
	parallel.SendUp_row(sendup_row, recvup_row, 3);
	parallel.Send_row(send_row, 3, (nr+1) % sc );
	parallel.Receive_row(recv_row, 3, ((nr-1)%sc + sc) % sc );
	status = compare_results(recv_row, recvup_row, 3, "Row send-up array") && status;

	//sendup-column
	const int& nc = parallel.get_grid_rank()[1];
	const int& sr = parallel.get_grid_size()[0];
	int sendup_col[] = { nc + 5, nc*nc, (nc-4)*(nc+3) };
	int recvup_col[] = {-1, -3, -6};
	int send_col[] = { nc + 5, nc*nc, (nc-4)*(nc+3) };
	int recv_col[] = {0, 0, 0};
	parallel.SendUp_col(sendup_col, recvup_col, 3);
	parallel.Send_col(send_col, 3, (nc+1) % sr );
	parallel.Receive_col(recv_col, 3, ((nc-1)%sr + sr) % sr );
	status = compare_results(recv_col, recvup_col, 3, "Col send-up array") & status;

	return status;
}

int Test_SendDown(const int n_rows, const int n_cols){
	Parallel2D parallel;
	parallel.InitializeGrid(n_rows, n_cols);
	int status = 1;

	//sendup-row
	const int& nr = parallel.get_grid_rank()[0];
	const int& sc = parallel.get_grid_size()[1];
	int sendup_row[] = { nr + 5, nr*nr, (nr-3)*(nr+4) };
	int recvup_row[] = {-1, -3, -6};
	int send_row[] = { nr + 5, nr*nr, (nr-3)*(nr+4) };
	int recv_row[] = {0, 0, 0};
	parallel.SendDown_row(sendup_row, recvup_row, 3);
	parallel.Send_row(send_row, 3, ((nr-1)%sc + sc) % sc );
	parallel.Receive_row(recv_row, 3, (nr+1) % sc );
	status = compare_results(recv_row, recvup_row, 3, "Row send-down array") && status;

	//sendup-column
	const int& nc = parallel.get_grid_rank()[1];
	const int& sr = parallel.get_grid_size()[0];
	int sendup_col[] = { nc + 5, nc*nc, (nc-4)*(nc+3) };
	int recvup_col[] = {-1, -3, -6};
	int send_col[] = { nc + 5, nc*nc, (nc-4)*(nc+3) };
	int recv_col[] = {0, 0, 0};
	parallel.SendDown_col(sendup_col, recvup_col, 3);
	parallel.Send_col(send_col, 3, ((nc-1)%sr + sr) % sr );
	parallel.Receive_col(recv_col, 3, (nc+1) % sr );
	status = compare_results(recv_col, recvup_col, 3, "Col send-down array") & status;

	return status;
}

int Test_Build_Lattice3D(
				const IndexType size[3], 
                const int halo, 
                const int node_size[2], 
                const int grid_loc[2]){
	int status = 1;
	Lattice<3> lat(size, halo, node_size, grid_loc);
	//check some global properties
	status = compare_results(3, lat.get_dim(), "Lattice dimension") && status;
	status = compare_results(size, lat.get_global_size(), 3, "Global size") && status;
	IndexType global_sites = std::accumulate(size, size + 3, 1, std::multiplies<IndexType>());
	status = compare_results(global_sites, lat.get_visible_global_sites(), "Visible global sites") && status;

	//check some local memory properties
	if(node_size[0] * node_size[1] == 1){
		IndexType mem_size[3] = {size[0] + 2*halo, size[1] + 2*halo, size[2] + 2*halo};
		IndexType total_mem_sites = std::accumulate(mem_size, mem_size + 3, 
									1, std::multiplies<IndexType>());
		status = compare_results<IndexType>(total_mem_sites, lat.get_local_mem_sites(), 
											"Total local memory sites") && status;
	}

	//check properties of local visible region
	auto first_local_vis = lat.get_local_visible_first();
	auto last_local_vis = lat.get_local_visible_next_to_last() - 1;
	IndexType first_local_vis_coord[3];
	IndexType last_local_vis_coord[3];
	lat.local_mem_index2coord(first_local_vis, first_local_vis_coord);
	IndexType theo_first[] = {lat.get_halo(), lat.get_halo(), lat.get_halo()};
	status = compare_results<IndexType>(theo_first, first_local_vis_coord, 3, 
								"First local visible site (local_mem_index)") && status;
	if(node_size[0] * node_size[1] == 1){
		IndexType mem_size[3] = {size[0] + 2*halo, size[1] + 2*halo, size[2] + 2*halo};
		lat.local_mem_index2coord(last_local_vis, last_local_vis_coord);
		IndexType theo_last[] = {mem_size[0] - 1 - halo, mem_size[1] - 1 - halo, mem_size[2] - 1 - halo};
		status = compare_results<IndexType>(theo_last, last_local_vis_coord, 3, 
								"Last local visible site (local_mem_index)") && status;
	}

	IndexType local_size[3];
	for(auto i = 0; i < 3; ++i){
		local_size[i] = last_local_vis_coord[i] - first_local_vis_coord[i] + 1;
	}
	
	status = compare_results<IndexType>(lat.get_local_size(), local_size, 3,
				"Computed local size") && status;
	return status;
}

int Test_Lattice3D_Global_Coordinate_Transformation(){
	IndexType size[] = {11,16,23};
	int node_size[] = {1,1};
	int loc[] = {0,0};
	Lattice<3> lat(size, 1, node_size, loc);
	int status = 1;

	IndexType global_sites = std::accumulate(size, size + 3, 1, std::multiplies<IndexType>());
	
	//check global coordinate and index transformation
	IndexType global_coord_1[3] = {0,0,0};
	IndexType new_coord[3];
	auto global_index = lat.global_coord2index(global_coord_1);
	status = compare_results<IndexType>(0, global_index, 
							"Global coordinate to index 1") && status;
	status = compare_results<IndexType>(global_coord_1, 
				lat.global_index2coord(global_index, new_coord), 3,
				"Global index to coordinate 1") && status;

	IndexType global_coord_2[3] = {size[0] - 1, size[1] - 1, size[2] - 1};
	global_index = lat.global_coord2index(global_coord_2);
	status = compare_results<IndexType>(global_sites - 1, global_index, 
							"Global coordinate to index 2") && status;
	status = compare_results<IndexType>(global_coord_2, 
				lat.global_index2coord(global_index, new_coord), 3, 
				"Global index to coordinate 2") && status;

	IndexType global_coord_3[3] = {0, 1, 0};
	global_index = lat.global_coord2index(global_coord_3);
	status = compare_results<IndexType>(size[0], global_index, 
							"Global coordinate to index 3") && status;
	status = compare_results<IndexType>(global_coord_3, 
				lat.global_index2coord(global_index, new_coord), 3, 
				"Global index to coordinate 3") && status;

	IndexType global_coord_4[3] = {1, 1, 2};
	global_index = lat.global_coord2index(global_coord_4);
	status = compare_results<IndexType>(1 + size[0] + 2*size[0]*size[1], 
							global_index, 
							"Global coordinate to index 4") && status;
	status = compare_results<IndexType>(global_coord_4, 
				lat.global_index2coord(global_index, new_coord), 3, 
				"Global index to coordinate 4") && status;

	return status;
}

int Test_Local_Visible_Loop_Single_Process(
	const int rounds, 
	const bool use_hash_table){
	IndexType size[] = {100,100,100};
	int node_size[] = {1,1};
	int loc[] = {0,0};
	Lattice<3> lat(size, 1, node_size, loc);
	int status = 1;
	if(use_hash_table){
		lat.make_vis2mem_index_table();
		lat.make_mem2vis_index_table();
	}
	for(auto i = 0; i < rounds; ++i){
		auto first_lv = lat.get_local_visible_first();
		auto next_last_lv = lat.get_local_visible_next_to_last();
		IndexType count = 0;
		IndexType lm_coord[3];
		for(auto i = first_lv; i < next_last_lv; i = lat.get_local_visible_next(i)){
			count++;
			status = compare_results(i, 
								lat.local_mem_coord2index(lat.local_mem_index2coord(i, lm_coord)), 
								"Visible sites loop index-coordinate transformation", false) && status;
		}
		status = compare_results<IndexType>(lat.get_visible_local_sites(), count, 
					"Visible sites loop count") && status;
	}
	return status;
}

int Test_Local_Visible_Loop_Pseudo_Multi_Processes(){
	IndexType size[] = {96,128,160};
	int node_size[] = {4,8};
	int loc[] = {3,5};
	Lattice<3> lat(size, 1, node_size, loc);
	int status = 1;

	auto first_lv = lat.get_local_visible_first();
	auto next_last_lv = lat.get_local_visible_next_to_last();
	IndexType count = 0;
	IndexType lm_coord[3];
	for(auto i = first_lv; i < next_last_lv; i = lat.get_local_visible_next(i)){
		count++;
		status = compare_results(i, 
							lat.local_mem_coord2index(lat.local_mem_index2coord(i, lm_coord)), 
							"Visible sites loop index-coordinate transformation", 
							false) && status;
	}
	status = compare_results<IndexType>(lat.get_visible_local_sites(), count, 
				"Visible sites loop count") && status;
	return status;
}

int Test_Local_Visible_Loop(const int n_rows, const int n_cols){
	Parallel2D parallel;
	parallel.InitializeGrid(n_rows, n_cols);
	auto grid_rank = parallel.get_grid_rank();
	int grid_loc[2];
	int node_size[] = {n_rows, n_cols};
	int status = 1;
	IndexType lat_size[] = {64, 128, 96};
	Lattice<3> lat(lat_size, 2, node_size,
					 transform_gridrank_to_gridloc(grid_rank, grid_loc));
	Site<3> x(lat);
	int count = 0;
	std::set<IndexType> lv_idx_set;
	for(x.first(); x.test(); x.next()){
		count++;
		auto lv_idx = lat.local_mem_index_to_local_vis_index(x.get_index());
		status = lat.is_local_visible(x.get_index()) && status;
		lv_idx_set.emplace(lv_idx);
	}
	status = compare_results<IndexType>(lv_idx_set.size(), lat.get_visible_local_sites(),
				"Different Local Visible Sites Count") && status;
	status = compare_results<IndexType>(lat.get_visible_local_sites(), count, 
				"Visible Loop Count") && status;
	return status;
}

int Test_Halo_Loop(const int n_rows, const int n_cols) {
	Parallel2D parallel;
	parallel.InitializeGrid(n_rows, n_cols);
	auto grid_rank = parallel.get_grid_rank();
	int grid_loc[2];
	int node_size[] = {n_rows, n_cols};
	int status = 1;
	IndexType lat_size[] = {64, 128, 96};
	Lattice<3> lat(lat_size, 2, node_size,
					 transform_gridrank_to_gridloc(grid_rank, grid_loc));
	Site<3> x(lat);
	int count = 0;
	std::set<IndexType> halo_idx_set;
	for(x.halo_first(); x.halo_test(); x.halo_next()){
		count++;
		halo_idx_set.emplace(x.get_index());
	}
	status = compare_results<IndexType>(lat.get_local_mem_sites() - lat.get_visible_local_sites(),
					halo_idx_set.size(), "Different Local Halo Sites Count") && status;
	status = compare_results<IndexType>(lat.get_local_mem_sites() - lat.get_visible_local_sites(),
					count, "Halo Loop Count") && status;

	//Test if the union of local visible sites and halo sites are
	//the same as local memory sites.
	std::set<IndexType> lv_idx_set;
	for(x.first(); x.test(); x.next()){
		lv_idx_set.emplace(x.get_index());
	}
	for(auto iter = halo_idx_set.begin(); iter != halo_idx_set.end(); ++iter){
		lv_idx_set.emplace(*iter);
	}
	IndexType s = 0;
	for(auto iter = lv_idx_set.begin(); iter != lv_idx_set.end(); ++iter){
		status = (*iter == s) && status;
		s++;
	}
	return status;
}

int Test_Site_Move(const int n_rows, const int n_cols){
	Parallel2D parallel;
	parallel.InitializeGrid(n_rows, n_cols);
	auto grid_rank = parallel.get_grid_rank();
	int grid_loc[2];
	int node_size[] = {n_rows, n_cols};
	int status = 1;
	IndexType lat_size[] = {64, 128, 96};
	Lattice<3> lat(lat_size, 2, node_size,
					 transform_gridrank_to_gridloc(grid_rank, grid_loc));
	Site<3> x(lat);

	x.first();

	IndexType mem_coord[3], vis_coord[3], glb_coord[3];
	x.coord_local_mem(mem_coord);
	print_coordinate(mem_coord, 3, "mem", false);
	if(lat.is_local_visible(x.get_index())){
		x.coord_local_vis(vis_coord);
		x.coord(glb_coord);
		print_coordinate(vis_coord, 3, "vis", false);
		print_coordinate(glb_coord, 3, "glb", true);
	}
	else {
		std::cout << std::endl;
	}

	for(auto _ = 0; _ < 10; ++_){
		x.move(1, -2).move(2,1).move(0,3);
		x.coord_local_mem(mem_coord);
		print_coordinate(mem_coord, 3, "mem", false);
		if(lat.is_local_visible(x.get_index())){
			x.coord_local_vis(vis_coord);
			x.coord(glb_coord);
			print_coordinate(vis_coord, 3, "vis", false);
			print_coordinate(glb_coord, 3, "glb", true);
		}
		else {
			std::cout << std::endl;
		}
	}
	return 1;
}

int Test_UpdateHalo(const int n_rows, const int n_cols){
	Parallel2D parallel;
	parallel.InitializeGrid(n_rows, n_cols);
	auto grid_rank = parallel.get_grid_rank();
	int grid_loc[2];
	int node_size[] = {n_rows, n_cols};
	int status = 1;
	IndexType lat_size[] = {16, 32, 24};
	Lattice<3> lat(lat_size, 2, node_size,
					 transform_gridrank_to_gridloc(grid_rank, grid_loc));
	Site<3> x(lat);
	Field<IndexType,3> fc(lat,3);
	fc.assign_parallel_object(parallel);
	for(x.first(); x.test(); x.next()){
		fc(x, 0) = x.coord(0);
		fc(x, 1) = x.coord(1);
		fc(x, 2) = x.coord(2);
	}
	for(x.halo_first(); x.halo_test(); x.halo_next()){
		fc(x, 0) = 0;
		fc(x, 1) = 0;
		fc(x, 2) = 0;
	}
	auto tab = fc.get_update_halo_table();
	fc.update_halo();
	
	for(IndexType i = 0; i < lat.get_local_mem_sites(); ++i){
		x.set_index(i);
		IndexType lm_coord[3], glb_coord[3];
		x.coord_local_mem(lm_coord);
		lat.local_mem_coord_to_global_coord(lm_coord, glb_coord);
		if(fc(x, 0) == glb_coord[0]
			&& fc(x, 1) == glb_coord[1]
			&& fc(x, 2) == glb_coord[2]) {
			IndexType fc_coord[3] = {fc(x,0), fc(x,1), fc(x, 2)};
			print_coordinate<IndexType>(glb_coord, 3, "Global coordinate", false);
			print_coordinate<IndexType>(fc_coord, 3, "fc coord");
			continue;
		}
		else{
			IndexType fc_coord[3] = {fc(x,0), fc(x,1), fc(x, 2)};
			print_coordinate<IndexType>(glb_coord, 3, "Global coordinate", false);
			print_coordinate<IndexType>(fc_coord, 3, "fc coord");
			status = 0;
		}
	}
	return status;
}

int Test_HDF5_WriteDataset() {
	HDF5_Wrapper::HDF5Wrapper h5("test1.h5");
	int data[100];
	std::generate_n(data, 100, [n = 0]() mutable {return static_cast<int>(n++);});
	hsize_t size[2] = {6,6};
	hsize_t mem_size[2] = {10,10};
	for(auto b = 0; b < 3; ++b){
		hsize_t file_offset[2] = {static_cast<hsize_t>(2*b),0};
		hsize_t mem_offset[2] = {static_cast<hsize_t>(2*b+1),3};
		hsize_t block_size[2] = {2,6};
		h5.SaveSingleDatasetFile<int,2>(
			"test.h5", 
			"int_data", 
			size, 
			file_offset, 
			block_size,
			mem_size,
			mem_offset,
			block_size, 
			data, b == 0);
	}
	return 1;
}

int Test_SerialOperation(){
	Parallel2D parallel(MPI_COMM_WORLD);
	auto my_rank = parallel.get_world_rank();
	auto world_size = parallel.get_world_size();
	int msg_send = 1;
	int msg_recv = 0;
	if(my_rank == 0){
		std::cout << "Job finished. " << my_rank << std::endl;
		if(world_size > 1)
			parallel.Send(msg_send, 1);
	}
	else{
		parallel.Receive(msg_recv, my_rank - 1);
		std::cout << "Job finished. " << my_rank << std::endl;

		if(my_rank < world_size -1){
			parallel.Send(msg_send, my_rank + 1);
		}

	}
	return 1;
}

int Test_ParaSiteFieldWrite(const int n_rows, const int n_cols){
	Parallel2D parallel(MPI_COMM_WORLD);
	parallel.InitializeGrid(n_rows, n_cols);
	GridIndexType node_size[] = { n_rows, n_cols };
	GridIndexType grid_rank[] = { parallel.get_grid_rank()[0], parallel.get_grid_rank()[1] };

	GridIndexType grid_loc[2];
	transform_gridrank_to_gridloc(grid_rank, grid_loc);
	IndexType lat_size[3] = {12,24,36};
	Lattice<3> lat(lat_size, 2, node_size, grid_loc);
	Site<3> x(lat);
	Field<int, 3> f(lat, 3, parallel);
	for(x.first(); x.test(); x.next()){
		f(x, 0) = x.coord(0);
		f(x, 1) = x.coord(1);
		f(x, 2) = x.coord(2);
	}
	f.update_halo();
	f.write("f.h5");
	return 1;
}
