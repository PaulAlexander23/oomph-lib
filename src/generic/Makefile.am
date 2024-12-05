
# A few file definitions
#-----------------------



# Define the sources
sources =  \
	Qelement_face_coordinate_translation_schemes.cc \
	Qelements.cc \
	Qspectral_elements.cc \
	Telements.cc \
	algebraic_elements.cc \
	assembly_handler.cc \
	binary_tree.cc \
	black_box_newton_solver.cc \
	block_preconditioner.cc \
	brick_mesh.cc \
	communicator.cc \
	complex_matrices.cc \
	dg_elements.cc \
	domain.cc \
	double_multi_vector.cc \
	double_vector.cc \
	double_vector_with_halo.cc \
	eigen_solver.cc \
	elastic_problems.cc \
	element_with_external_element.cc \
	element_with_moving_nodes.cc \
	elements.cc \
	error_estimator.cc \
	explicit_timesteppers.cc \
	extruded_domain.cc \
	extruded_macro_element.cc \
	frontal_solver.cc \
	fsi.cc \
	general_purpose_block_preconditioners.cc \
	general_purpose_preconditioners.cc \
	geometric_multigrid.cc \
	geompack_scaffold_mesh.cc \
	hermite_elements.cc \
	hijacked_elements.cc \
	hp_refineable_elements.cc \
	implicit_midpoint_rule.cc \
	integral.cc \
	iterative_linear_solver.cc \
	line_mesh.cc \
	linear_algebra_distribution.cc \
	linear_solver.cc \
	macro_element.cc \
	macro_element_node_update_element.cc \
	matrices.cc \
	matrix_vector_product.cc \
	mesh.cc \
	mesh_as_geometric_object.cc \
	missing_masters.cc \
	multi_domain.cc \
	nodes.cc \
	octree.cc \
	oomph_definitions.cc \
	oomph_utilities.cc \
	orthpoly.cc \
	partitioning.cc \
	periodic_orbit_handler.cc \
	pml_meshes.cc \
	preconditioner_array.cc \
	problem.cc \
	pseudosolid_node_update_elements.cc \
	quad_mesh.cc \
	quadtree.cc \
	refineable_brick_element.cc \
	refineable_elements.cc \
	refineable_line_element.cc \
	refineable_mesh.cc \
	refineable_quad_element.cc \
	sample_point_container.cc \
	sample_point_parameters.cc \
	simple_cubic_scaffold_tet_mesh.cc \
	spines.cc \
	stored_shape_function_elements.cc \
	sum_of_matrices.cc \
	superlu.c \
	superlu_complex.c \
	tet_mesh.cc \
	tetgen_scaffold_mesh.cc \
	timesteppers.cc \
	tree.cc \
	triangle_mesh.cc \
	triangle_scaffold_mesh.cc \
	unstructured_two_d_mesh_geometry_base.cc

if OOMPH_HAS_MUMPS
sources+=mumps_fortran_solver.F mumps_solver.cc
endif

if WANT_MPI
sources += superlu_dist.c
endif


# Define the headers
headers =  \
	Qelement_face_coordinate_translation_schemes.h \
	Qelements.h \
	Qspectral_elements.h \
	SuperLU_preconditioner.h \
	Telements.h \
	Vector.h \
	algebraic_elements.h \
	assembly_handler.h \
	binary_tree.h \
	black_box_newton_solver.h \
	block_preconditioner.h \
	brick_mesh.h \
	communicator.h \
	complex_matrices.h \
	dg_elements.h \
	displacement_control_element.h \
	domain.h \
	double_multi_vector.h \
	double_vector.h \
	double_vector_with_halo.h \
	eigen_solver.h \
	elastic_problems.h \
	element_with_error.h \
	element_with_external_element.h \
	element_with_moving_nodes.h \
	elements.h \
	error_estimator.h \
	explicit_timesteppers.h \
	extruded_domain.h \
	extruded_macro_element.h \
	face_element_as_geometric_object.h \
	face_mesh_project.h \
	frontal_solver.h \
	fsi.h \
	general_purpose_block_preconditioners.h \
	general_purpose_preconditioners.h \
	generalised_newtonian_constitutive_models.h \
	generalised_timesteppers.h \
	geom_obj_with_boundary.h \
	geom_objects.h \
	geometric_multigrid.h \
	geompack_scaffold_mesh.h \
	hermite_elements.h \
	hijacked_elements.h \
	hp_refineable_elements.h \
	implicit_midpoint_rule.h \
	integral.h \
	iterative_linear_solver.h \
	line_mesh.h \
	line_visualiser.h \
	linear_algebra_distribution.h \
	linear_solver.h \
	macro_element.h \
	macro_element_node_update_element.h \
	map_matrix.h \
	matrices.h \
	matrix_vector_product.h \
	mesh.h \
	mesh_as_geometric_object.h \
	missing_masters.h \
	multi_domain.h \
	nodes.h \
	octree.h \
	oomph_definitions.h \
	oomph_utilities.h \
	orthpoly.h \
	partitioning.h \
	periodic_orbit_handler.h \
	pml_mapping_functions.h \
	pml_meshes.h \
	preconditioner.h \
	preconditioner_array.h \
	problem.h \
	projection.h \
	pseudo_buckling_ring.h \
	pseudosolid_node_update_elements.h \
	quad_mesh.h \
	quadtree.h \
	refineable_brick_element.h \
	refineable_brick_mesh.h \
	refineable_brick_spectral_element.h \
	refineable_elements.h \
	refineable_line_element.h \
	refineable_line_mesh.h \
	refineable_line_spectral_element.h \
	refineable_mesh.h \
	refineable_quad_element.h \
	refineable_quad_mesh.h \
	refineable_quad_spectral_element.h \
	sample_point_container.h \
	sample_point_parameters.h \
	shape.h \
	simple_cubic_scaffold_tet_mesh.h \
	sparse_vector.h \
	spines.h \
	stored_shape_function_elements.h \
	sum_of_matrices.h \
	tet_mesh.h \
	tetgen_scaffold_mesh.h \
	timesteppers.h \
	trapezoid_rule.h \
	tree.h \
	triangle_mesh.h \
	triangle_scaffold_mesh.h \
	unstructured_two_d_mesh_geometry_base.h \
	vector_matrix.h


if OOMPH_HAS_MUMPS
headers+=mumps_solver.h
endif

if HAVE_STACKTRACE
headers += stacktrace.h
endif

# Template-only files: Need to be included in include directory
# but should not be listed in combined header file
templates = \
tree.template.cc refineable_mesh.template.cc multi_domain.template.cc \
missing_masters.template.cc

# Do we have the HYPRE library?
if OOMPH_HAS_HYPRE

   # Add the hypre sources
   sources += hypre_solver.cc


   # Add the hypre headers
   headers += hypre_solver.h

endif

# Do we have the TRILINOS library?
if OOMPH_HAS_TRILINOS

   # Add the trilinos sources
   sources += trilinos_helpers.cc trilinos_preconditioners.cc \
              trilinos_solver.cc trilinos_eigen_solver.cc


   # Add the trilinos headers
   headers +=   trilinos_helpers.h trilinos_preconditioners.h \
                trilinos_solver.h trilinos_eigen_solver.h

endif

# Need to worry about these for dependency but they shouldn't
# be included into the include directory
extra_headers = cfortran.h frontal.h lapack_qz.h mumps.h

# Define name of library
libname = generic

# Combine headers and sources
headers_and_sources = $(headers) $(extra_headers) $(sources)


# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libgeneric.la


# Sources that the library depends on:
#-------------------------------------
libgeneric_la_SOURCES = $(headers_and_sources)


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/generic

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) $(templates)

# Required libraries
#-------------------
# Include path for library headers -- need to refer to include files
# in src as they will not have been installed yet!
AM_CPPFLAGS += -I$(top_srcdir)/external_src -I$(top_srcdir)/external_src/oomph_tetgen



# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time
$(libname).h: dummy_$(libname).h

dummy_$(libname).h:  $(headers)
	echo  $(libname) $(headers) > all_$(libname).aux
	$(AWK) -f $(top_srcdir)/bin/headers.awk < \
	       all_$(libname).aux > $(libname).h
	rm all_$(libname).aux


# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook:
	(cd $(library_includedir) && rm -f $(headers) $(templates) )
	(echo "$(headers) $(templates)" > include_files.list )
	($(top_srcdir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
	(rm -r include_files.list.aux include_files.list )
else
install-data-hook:
endif


# Cleanup: Kill automatically generated combined header file
clean-local:
	rm -f $(libname).h *~
distclean-local:
	rm -f $(libname).h *~

#EXTRA_DIST=prettyprint.h prettyprint98.h
