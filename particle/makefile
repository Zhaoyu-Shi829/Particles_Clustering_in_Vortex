#include ../mglet.ifc.pc
#include ../mglet.gfort.pc
#include ../mglet.ia64
# Preprocession options
CPP=$(MGLET_CPP)
CPPFLAGS=$(MGLET_CPPFLAGS)

# C compilation options
CC=$(MGLET_CC)
CFLAGS=$(MGLET_CFLAGS)

# Fortran compilation options
FC=$(MGLET_FC)
FFLAGS=$(MGLET_FFLAGS) -I..
F90FLAGS=$(MGLET_F90FLAGS)

# Linking options
LFLAGS=$(MGLET_LFLAGS)
LLIBS=$(MGLET_LLIBS)


NAME= ppack


OBJECTS= \
bounded_gauss.o \
correctvar.o \
pp_add_contactpart.o \
pp_add_resting_part.o \
pp_adm.o \
pp_basb.o \
pp_block_copy_grid.o \
pp_bound.o \
pp_bound_body.o \
pp_can_move.o \
pp_change_grid.o \
pp_change_proc.o \
pp_check_bed.o \
pp_check_collision.o \
pp_check_collision_circlepath.o \
pp_check_overlap.o \
pp_check_resting_part.o \
pp_coarse_grid.o \
pp_collisions.o \
pp_collisions_overlap.o \
pp_compute_dist.o \
pp_cutoff_filtervel.o \
pp_cylindrical_coords.o \
pp_cylindrical_coords_inv.o \
pp_cylindrical_forces.o \
pp_cylindrical_forces_inv.o \
pp_dep_rates.o \
pp_del_parts.o \
pp_diam_from_file.o \
pp_dim_arrays.o \
pp_filter_generic.o \
pp_filter_periodicx_ngrd.o \
pp_filter_periodicy_ngrd.o \
pp_filterwidth.o \
pp_find_resting_particles.o \
pp_fix_resting_part.o \
pp_forces.o \
pp_getvelgradient1d.o \
pp_getvelocities_4thorder.o \
pp_getvelocities_4thorder_onecomponent.o \
pp_getvelocities_generic.o \
pp_getvelocities_cl_cs.o \
pp_getvelocities_sl_cs.o \
pp_get_shift.o \
pp_getpartsfromcoarse.o \
pp_ibasb.o \
pp_icmf07.o \
pp_idx_to_number.o \
pp_init.o \
pp_init_single.o \
pp_inlet_randnum.o \
pp_move_part.o \
pp_mgftoc.o \
pp_newindex.o \
pp_num_vars.o \
pp_new_parts.o \
pp_ode_upart_d_rhs_t.o \
pp_ode_upart_dummy.o \
pp_ode_upart_jac.o \
pp_ode_upart_rhs.o \
pp_ode_upart_solve.o \
pp_ode_upart_store_sol.o \
pp_part_diameter.o \
pp_rand_part_pos.o \
pp_rand_part_pos_single.o \
pp_parbound.o  \
pp_read.o \
pp_read_bed.o \
pp_read_bedsum.o \
pp_read_bedsum_grid.o \
pp_read_section.o \
pp_remove_bed_part.o \
pp_remove_from_contactlist.o \
pp_resting_part_to_bedsum.o \
pp_saffman_accel.o \
pp_saffman_tabulated.o \
pp_schmeeckle_bound.o \
pp_set_shift.o \
pp_sgs.o \
pp_sgs_fluidvel.o \
pp_sommerfeld_wall.o \
pp_sort_in_grids.o \
pp_spherical_coords.o \
pp_spherical_coords_inv.o \
pp_spherical_forces.o \
pp_spherical_forces_inv.o \
pp_stokes_drag.o \
pp_terminate.o \
pp_tophatfilter.o \
pp_tophatfiltervel.o \
pp_transport.o \
pp_transport_grid.o \
pp_tstpos.o \
pp_tstvelpart.o \
pp_tst_cylindrical_pos.o \
pp_tst_spherical_pos.o \
pp_update_ds_stat.o \
pp_write.o \
pp_write_bed.o \
pp_write_bedsum.o \
pp_write_bedsum_grid.o \
pp_write_section.o \
pp_writedefs.o \
solve_quad_eq.o \
sros4.o \
pp_fields_mod.o \
pp_settings_mod.o \
pp_statistics_mod.o \
finecell_fluidvol.o



libppack.a:  $(OBJECTS) ppack.def 
	ar ruv libppack.a $(OBJECTS)
	ranlib libppack.a

bounded_gauss.o: bounded_gauss.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def bounded_gauss.F |$(CPP) $(CPPFLAGS)  > bounded_gauss.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bounded_gauss.f

correctvar.o: correctvar.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def correctvar.F |$(CPP) $(CPPFLAGS)  > correctvar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  correctvar.f

pp_add_contactpart.o: pp_add_contactpart.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_add_contactpart.F |$(CPP) $(CPPFLAGS)  > pp_add_contactpart.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_add_contactpart.f

pp_add_resting_part.o: pp_add_resting_part.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_add_resting_part.F |$(CPP) $(CPPFLAGS)  > pp_add_resting_part.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_add_resting_part.f

pp_adm.o: pp_adm.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_adm.F |$(CPP) $(CPPFLAGS)  > pp_adm.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_adm.f

pp_basb.o: pp_basb.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_basb.F |$(CPP) $(CPPFLAGS)  > pp_basb.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_basb.f

pp_block_copy_grid.o: pp_block_copy_grid.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_block_copy_grid.F |$(CPP) $(CPPFLAGS)  > pp_block_copy_grid.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_block_copy_grid.f

pp_bound.o: pp_bound.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_bound.F |$(CPP) $(CPPFLAGS)  > pp_bound.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_bound.f

pp_bound_body.o: pp_bound_body.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_bound_body.F |$(CPP) $(CPPFLAGS)  > pp_bound_body.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_bound_body.f

pp_can_move.o: pp_can_move.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_can_move.F |$(CPP) $(CPPFLAGS)  > pp_can_move.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_can_move.f

pp_change_grid.o: pp_change_grid.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_change_grid.F |$(CPP) $(CPPFLAGS)  > pp_change_grid.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_change_grid.f

pp_change_proc.o: pp_change_proc.F ppack.def ../mg41.def pp_settings_mod.o
	cat ../mg41.def ppack.def pp_change_proc.F |$(CPP) $(CPPFLAGS)  > pp_change_proc.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_change_proc.f

pp_check_bed.o: pp_check_bed.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_check_bed.F |$(CPP) $(CPPFLAGS)  > pp_check_bed.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_check_bed.f
         
pp_check_collision.o: pp_check_collision.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_check_collision.F |$(CPP) $(CPPFLAGS)  > pp_check_collision.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_check_collision.f

pp_check_collision_circlepath.o: pp_check_collision_circlepath.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_check_collision_circlepath.F |$(CPP) $(CPPFLAGS)  > pp_check_collision_circlepath.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_check_collision_circlepath.f

pp_check_overlap.o: pp_check_overlap.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_check_overlap.F |$(CPP) $(CPPFLAGS)  > pp_check_overlap.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_check_overlap.f

pp_check_resting_part.o: pp_check_resting_part.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_check_resting_part.F |$(CPP) $(CPPFLAGS)  > pp_check_resting_part.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_check_resting_part.f

pp_coarse_grid.o: pp_coarse_grid.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_coarse_grid.F |$(CPP) $(CPPFLAGS)  > pp_coarse_grid.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_coarse_grid.f

pp_collisions.o: pp_collisions.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_collisions.F |$(CPP) $(CPPFLAGS)  > pp_collisions.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_collisions.f

pp_collisions_overlap.o: pp_collisions_overlap.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_collisions_overlap.F |$(CPP) $(CPPFLAGS)  > pp_collisions_overlap.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_collisions_overlap.f

pp_compute_dist.o: pp_compute_dist.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_compute_dist.F |$(CPP) $(CPPFLAGS)  > pp_compute_dist.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_compute_dist.f

pp_cutoff_filtervel.o: pp_cutoff_filtervel.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_cutoff_filtervel.F |$(CPP) $(CPPFLAGS)  > pp_cutoff_filtervel.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_cutoff_filtervel.f
         
pp_cylindrical_coords.o: pp_cylindrical_coords.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_cylindrical_coords.F |$(CPP) $(CPPFLAGS)  > pp_cylindrical_coords.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_cylindrical_coords.f
         
pp_cylindrical_coords_inv.o: pp_cylindrical_coords_inv.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_cylindrical_coords_inv.F |$(CPP) $(CPPFLAGS)  > pp_cylindrical_coords_inv.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_cylindrical_coords_inv.f
         
pp_cylindrical_forces.o: pp_cylindrical_forces.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_cylindrical_forces.F |$(CPP) $(CPPFLAGS)  > pp_cylindrical_forces.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_cylindrical_forces.f
         
pp_cylindrical_forces_inv.o: pp_cylindrical_forces_inv.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_cylindrical_forces_inv.F |$(CPP) $(CPPFLAGS)  > pp_cylindrical_forces_inv.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_cylindrical_forces_inv.f

pp_dep_rates.o: pp_dep_rates.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_dep_rates.F |$(CPP) $(CPPFLAGS)  > pp_dep_rates.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_dep_rates.f

pp_del_parts.o: pp_del_parts.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_del_parts.F |$(CPP) $(CPPFLAGS)  > pp_del_parts.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_del_parts.f

pp_diam_from_file.o: pp_diam_from_file.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_diam_from_file.F |$(CPP) $(CPPFLAGS)  > pp_diam_from_file.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_diam_from_file.f

pp_dim_arrays.o: pp_dim_arrays.F ppack.def ../mg41.def pp_settings_mod.o ../coicells.h ibfieldscc_mod.o fields_mod.o precision_mod.o pointer_mod.o
	cat ../mg41.def ppack.def pp_dim_arrays.F |$(CPP) $(CPPFLAGS)  > pp_dim_arrays.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_dim_arrays.f

pp_filter_generic.o: pp_filter_generic.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_filter_generic.F |$(CPP) $(CPPFLAGS)  > pp_filter_generic.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_filter_generic.f

pp_filter_periodicx_ngrd.o: pp_filter_periodicx_ngrd.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_filter_periodicx_ngrd.F |$(CPP) $(CPPFLAGS)  > pp_filter_periodicx_ngrd.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_filter_periodicx_ngrd.f

pp_filter_periodicy_ngrd.o: pp_filter_periodicy_ngrd.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_filter_periodicy_ngrd.F |$(CPP) $(CPPFLAGS)  > pp_filter_periodicy_ngrd.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_filter_periodicy_ngrd.f

pp_filterwidth.o: pp_filterwidth.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_filterwidth.F |$(CPP) $(CPPFLAGS)  > pp_filterwidth.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_filterwidth.f

pp_find_resting_particles.o: pp_find_resting_particles.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_find_resting_particles.F |$(CPP) $(CPPFLAGS)  > pp_find_resting_particles.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_find_resting_particles.f

pp_fix_resting_part.o: pp_fix_resting_part.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_fix_resting_part.F |$(CPP) $(CPPFLAGS)  > pp_fix_resting_part.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_fix_resting_part.f

pp_forces.o: pp_forces.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_forces.F |$(CPP) $(CPPFLAGS)  > pp_forces.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_forces.f

pp_getvelocities_cl_cs.o: pp_getvelocities_cl_cs.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_getvelocities_cl_cs.F |$(CPP) $(CPPFLAGS)  > pp_getvelocities_cl_cs.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_getvelocities_cl_cs.f

pp_getvelgradient1d.o: pp_getvelgradient1d.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_getvelgradient1d.F |$(CPP) $(CPPFLAGS)  > pp_getvelgradient1d.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_getvelgradient1d.f

pp_getvelocities_4thorder.o: pp_getvelocities_4thorder.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_getvelocities_4thorder.F |$(CPP) $(CPPFLAGS)  > pp_getvelocities_4thorder.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_getvelocities_4thorder.f

pp_getvelocities_4thorder_onecomponent.o: pp_getvelocities_4thorder_onecomponent.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_getvelocities_4thorder_onecomponent.F |$(CPP) $(CPPFLAGS)  > pp_getvelocities_4thorder_onecomponent.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_getvelocities_4thorder_onecomponent.f

pp_getvelocities_generic.o: pp_getvelocities_generic.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_getvelocities_generic.F |$(CPP) $(CPPFLAGS)  > pp_getvelocities_generic.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_getvelocities_generic.f

pp_getvelocities_sl_cs.o: pp_getvelocities_sl_cs.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_getvelocities_sl_cs.F |$(CPP) $(CPPFLAGS)  > pp_getvelocities_sl_cs.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_getvelocities_sl_cs.f

pp_get_shift.o: pp_get_shift.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_get_shift.F |$(CPP) $(CPPFLAGS)  > pp_get_shift.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_get_shift.f

pp_getpartsfromcoarse.o: pp_getpartsfromcoarse.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_getpartsfromcoarse.F |$(CPP) $(CPPFLAGS)  > pp_getpartsfromcoarse.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_getpartsfromcoarse.f

pp_ibasb.o: pp_ibasb.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_ibasb.F |$(CPP) $(CPPFLAGS)  > pp_ibasb.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_ibasb.f

pp_icmf07.o: pp_icmf07.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_icmf07.F |$(CPP) $(CPPFLAGS)  > pp_icmf07.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_icmf07.f

pp_idx_to_number.o: pp_idx_to_number.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_idx_to_number.F |$(CPP) $(CPPFLAGS)  > pp_idx_to_number.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_idx_to_number.f

pp_init.o: pp_init.F ppack.def pp_fields_mod.o ../mg41.def pp_settings_mod.o pointer_mod.o pp_statistics_mod.o
	cat ../mg41.def ppack.def pp_init.F |$(CPP) $(CPPFLAGS)  > pp_init.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_init.f

pp_init_single.o: pp_init_single.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_init_single.F |$(CPP) $(CPPFLAGS)  > pp_init_single.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_init_single.f

pp_inlet_randnum.o: pp_inlet_randnum.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_inlet_randnum.F |$(CPP) $(CPPFLAGS)  > pp_inlet_randnum.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_inlet_randnum.f

pp_move_part.o: pp_move_part.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_move_part.F |$(CPP) $(CPPFLAGS)  > pp_move_part.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_move_part.f

pp_mgftoc.o: pp_mgftoc.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_mgftoc.F |$(CPP) $(CPPFLAGS)  > pp_mgftoc.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_mgftoc.f

pp_newindex.o: pp_newindex.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_newindex.F |$(CPP) $(CPPFLAGS)  > pp_newindex.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_newindex.f

pp_num_vars.o: pp_num_vars.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_num_vars.F |$(CPP) $(CPPFLAGS)  > pp_num_vars.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_num_vars.f

pp_new_parts.o: pp_new_parts.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_new_parts.F |$(CPP) $(CPPFLAGS)  > pp_new_parts.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_new_parts.f

pp_ode_upart_d_rhs_t.o: pp_ode_upart_d_rhs_t.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_ode_upart_d_rhs_t.F |$(CPP) $(CPPFLAGS)  > pp_ode_upart_d_rhs_t.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_ode_upart_d_rhs_t.f

pp_ode_upart_dummy.o: pp_ode_upart_dummy.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_ode_upart_dummy.F |$(CPP) $(CPPFLAGS)  > pp_ode_upart_dummy.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_ode_upart_dummy.f

pp_ode_upart_jac.o: pp_ode_upart_jac.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_ode_upart_jac.F |$(CPP) $(CPPFLAGS)  > pp_ode_upart_jac.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_ode_upart_jac.f

pp_ode_upart_rhs.o: pp_ode_upart_rhs.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_ode_upart_rhs.F |$(CPP) $(CPPFLAGS)  > pp_ode_upart_rhs.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_ode_upart_rhs.f

pp_ode_upart_solve.o: pp_ode_upart_solve.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_ode_upart_solve.F |$(CPP) $(CPPFLAGS)  > pp_ode_upart_solve.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_ode_upart_solve.f

pp_ode_upart_store_sol.o: pp_ode_upart_store_sol.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_ode_upart_store_sol.F |$(CPP) $(CPPFLAGS)  > pp_ode_upart_store_sol.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_ode_upart_store_sol.f

pp_part_diameter.o: pp_part_diameter.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_part_diameter.F |$(CPP) $(CPPFLAGS)  > pp_part_diameter.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_part_diameter.f

pp_rand_part_pos.o: pp_rand_part_pos.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_rand_part_pos.F |$(CPP) $(CPPFLAGS)  > pp_rand_part_pos.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_rand_part_pos.f

pp_rand_part_pos_single.o: pp_rand_part_pos_single.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_rand_part_pos_single.F |$(CPP) $(CPPFLAGS)  > pp_rand_part_pos_single.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_rand_part_pos_single.f

pp_parbound.o: pp_parbound.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_parbound.F |$(CPP) $(CPPFLAGS)  > pp_parbound.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_parbound.f

pp_read.o: pp_read.F ppack.def ../mg41.def pp_settings_mod.o pointer_mod.o fields_mod.o
	cat ../mg41.def ppack.def pp_read.F |$(CPP) $(CPPFLAGS)  > pp_read.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_read.f

pp_read_bed.o: pp_read_bed.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_read_bed.F |$(CPP) $(CPPFLAGS)  > pp_read_bed.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_read_bed.f

pp_read_bedsum.o: pp_read_bedsum.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_read_bedsum.F |$(CPP) $(CPPFLAGS)  > pp_read_bedsum.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_read_bedsum.f

pp_read_bedsum_grid.o: pp_read_bedsum_grid.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_read_bedsum_grid.F |$(CPP) $(CPPFLAGS)  > pp_read_bedsum_grid.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_read_bedsum_grid.f

pp_read_section.o: pp_read_section.F ppack.def ../mg41.def pp_settings_mod.o
	cat ../mg41.def ppack.def pp_read_section.F |$(CPP) $(CPPFLAGS)  > pp_read_section.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_read_section.f

pp_remove_bed_part.o: pp_remove_bed_part.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_remove_bed_part.F |$(CPP) $(CPPFLAGS)  > pp_remove_bed_part.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_remove_bed_part.f

pp_remove_from_contactlist.o: pp_remove_from_contactlist.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_remove_from_contactlist.F |$(CPP) $(CPPFLAGS)  > pp_remove_from_contactlist.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_remove_from_contactlist.f

pp_resting_part_to_bedsum.o: pp_resting_part_to_bedsum.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_resting_part_to_bedsum.F |$(CPP) $(CPPFLAGS)  > pp_resting_part_to_bedsum.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_resting_part_to_bedsum.f

pp_saffman_accel.o: pp_saffman_accel.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_saffman_accel.F |$(CPP) $(CPPFLAGS)  > pp_saffman_accel.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_saffman_accel.f

pp_saffman_tabulated.o: pp_saffman_tabulated.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_saffman_tabulated.F |$(CPP) $(CPPFLAGS)  > pp_saffman_tabulated.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_saffman_tabulated.f

pp_schmeeckle_bound.o: pp_schmeeckle_bound.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_schmeeckle_bound.F |$(CPP) $(CPPFLAGS)  > pp_schmeeckle_bound.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_schmeeckle_bound.f

pp_set_shift.o: pp_set_shift.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_set_shift.F |$(CPP) $(CPPFLAGS)  > pp_set_shift.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_set_shift.f

pp_sgs.o: pp_sgs.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_sgs.F |$(CPP) $(CPPFLAGS)  > pp_sgs.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_sgs.f

pp_sgs_fluidvel.o: pp_sgs_fluidvel.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_sgs_fluidvel.F |$(CPP) $(CPPFLAGS)  > pp_sgs_fluidvel.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_sgs_fluidvel.f

pp_sommerfeld_wall.o: pp_sommerfeld_wall.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_sommerfeld_wall.F |$(CPP) $(CPPFLAGS)  > pp_sommerfeld_wall.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_sommerfeld_wall.f

pp_sort_in_grids.o: pp_sort_in_grids.F ppack.def ../mg41.def pointer_mod.o fields_mod.o ../colevel.h
	cat ../mg41.def ppack.def pp_sort_in_grids.F |$(CPP) $(CPPFLAGS)  > pp_sort_in_grids.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_sort_in_grids.f

pp_spherical_coords.o: pp_spherical_coords.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_spherical_coords.F |$(CPP) $(CPPFLAGS)  > pp_spherical_coords.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_spherical_coords.f
         
pp_spherical_coords_inv.o: pp_spherical_coords_inv.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_spherical_coords_inv.F |$(CPP) $(CPPFLAGS)  > pp_spherical_coords_inv.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_spherical_coords_inv.f
         
pp_spherical_forces.o: pp_spherical_forces.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_spherical_forces.F |$(CPP) $(CPPFLAGS)  > pp_spherical_forces.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_spherical_forces.f
         
pp_spherical_forces_inv.o: pp_spherical_forces_inv.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_spherical_forces_inv.F |$(CPP) $(CPPFLAGS)  > pp_spherical_forces_inv.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_spherical_forces_inv.f
         
pp_stokes_drag.o: pp_stokes_drag.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_stokes_drag.F |$(CPP) $(CPPFLAGS)  > pp_stokes_drag.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_stokes_drag.f

pp_terminate.o: pp_terminate.F ppack.def ../mg41.def pp_statistics_mod.o
	cat ../mg41.def ppack.def pp_terminate.F |$(CPP) $(CPPFLAGS)  > pp_terminate.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_terminate.f

pp_transport.o: pp_transport.F pp_fields_mod.o ppack.def ../mg41.def pp_statistics_mod.o
	cat ../mg41.def ppack.def pp_transport.F |$(CPP) $(CPPFLAGS)  > pp_transport.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_transport.f

pp_transport_grid.o: pp_transport_grid.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_transport_grid.F |$(CPP) $(CPPFLAGS)  > pp_transport_grid.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_transport_grid.f

pp_tophatfilter.o: pp_tophatfilter.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_tophatfilter.F |$(CPP) $(CPPFLAGS)  > pp_tophatfilter.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_tophatfilter.f

pp_tophatfiltervel.o: pp_tophatfiltervel.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_tophatfiltervel.F |$(CPP) $(CPPFLAGS)  > pp_tophatfiltervel.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_tophatfiltervel.f

pp_tst_cylindrical_pos.o: pp_tst_cylindrical_pos.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_tst_cylindrical_pos.F |$(CPP) $(CPPFLAGS)  > pp_tst_cylindrical_pos.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_tst_cylindrical_pos.f

pp_tst_spherical_pos.o: pp_tst_spherical_pos.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_tst_spherical_pos.F |$(CPP) $(CPPFLAGS)  > pp_tst_spherical_pos.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_tst_spherical_pos.f
         
pp_tstpos.o: pp_tstpos.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_tstpos.F |$(CPP) $(CPPFLAGS)  > pp_tstpos.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_tstpos.f
         
pp_tstvelpart.o: pp_tstvelpart.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_tstvelpart.F |$(CPP) $(CPPFLAGS)  > pp_tstvelpart.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_tstvelpart.f

pp_update_ds_stat.o: pp_update_ds_stat.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_update_ds_stat.F |$(CPP) $(CPPFLAGS)  > pp_update_ds_stat.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_update_ds_stat.f

pp_write.o: pp_write.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_write.F |$(CPP) $(CPPFLAGS)  > pp_write.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_write.f

pp_write_bed.o: pp_write_bed.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_write_bed.F |$(CPP) $(CPPFLAGS)  > pp_write_bed.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_write_bed.f

pp_write_bedsum.o: pp_write_bedsum.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_write_bedsum.F |$(CPP) $(CPPFLAGS)  > pp_write_bedsum.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_write_bedsum.f

pp_write_bedsum_grid.o: pp_write_bedsum_grid.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_write_bedsum_grid.F |$(CPP) $(CPPFLAGS)  > pp_write_bedsum_grid.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_write_bedsum_grid.f

pp_write_section.o: pp_write_section.F ppack.def ../mg41.def pp_settings_mod.o
	cat ../mg41.def ppack.def pp_write_section.F |$(CPP) $(CPPFLAGS)  > pp_write_section.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_write_section.f

pp_writedefs.o: pp_writedefs.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def pp_writedefs.F |$(CPP) $(CPPFLAGS)  > pp_writedefs.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pp_writedefs.f

sros4.o: sros4.F90
	 $(F77) $(OPTIONS) $(FFLAGS)  sros4.F90

solve_quad_eq.o: solve_quad_eq.F ppack.def ../mg41.def
	cat ../mg41.def ppack.def solve_quad_eq.F |$(CPP) $(CPPFLAGS)  > solve_quad_eq.f
	 $(F77) $(OPTIONS) $(FFLAGS)  solve_quad_eq.f

pp_fields_mod.o: pp_fields_mod.F90 ppack.def ../mg41.def allocator_mod.o
	$(FC) $(FFLAGS) $(F90FLAGS) pp_fields_mod.F90

pp_settings_mod.o: pp_settings_mod.F90 ppack.def ../mg41.def precision_mod.o
	$(FC) $(FFLAGS) $(F90FLAGS) pp_settings_mod.F90

pp_statistics_mod.o: pp_statistics_mod.F90 precision_mod.o pointer_mod.o h5fieldio_mod.o
	$(FC) $(FFLAGS) $(F90FLAGS) pp_statistics_mod.F90

finecell_fluidvol.o: finecell_fluidvol.F90 precision_mod.o
	$(FC) $(FFLAGS) $(F90FLAGS) finecell_fluidvol.F90

allocator_mod.o:
	make -C .. allocator_mod.o
precision_mod.o:
	make -C .. precision_mod.o
ibfieldscc_mod.o:
	make -C .. ibfieldscc_mod.o
pointer_mod.o:
	make -C .. pointer_mod.o
h5fieldio_mod.o:
	make -C .. h5fieldio_mod.o
fields_mod.o:
	make -C .. fields_mod.o

clean :
	rm -f *.o *.f write_switches.tmp *~

clobber : clean
	rm -f libppack.a *~


