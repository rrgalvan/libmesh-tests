// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// <h1>Interior Penalty Discontinuous Galerkin</h1>
// \author Lorenzo Botti
// \date 2010
//
// This example is based on Adaptivity Example 3, but uses an
// Interior Penalty Discontinuous Galerkin formulation.

#include <iostream>

// LibMesh include files.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe_interface.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/discontinuity_measure.h"
#include "libmesh/string_to_enum.h"

#include "libmesh/exact_solution.h"

#include "fe_wrapper.h"
#include "dg_facecoupling.h"
// #include "dg_faceterm.h"

//#define QORDER TWENTYSIXTH

// Bring in everything from the libMesh namespace
using namespace libMesh;

Number exact_solution (const Point & p,
                       const Parameters &,
                       const std::string &,
                       const std::string &)
{
  const Real x = p(0);
  const Real y = p(1);
  // const Real z = p(2);

  // return (x+1)*(1-y*y)*(1-z*z);
  return (x+1)*(1-y*y);
}

// We now define the gradient of the exact solution, again being careful
// to obtain an angle from atan2 in the correct
// quadrant.
Gradient exact_derivative(const Point & p,
                          const Parameters &,  // es parameters, not needed
                          const std::string &,            // sys_name, not needed
                          const std::string &)            // unk_name, not needed
{
  // Gradient value to be returned.
  Gradient gradu;

  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  // const Real z = p(2);

  gradu(0) = (1-y*y);
  gradu(1) = -2*y*(x+1);
  gradu(2) = 0;

  return gradu;
}

Number exact_laplacian ( const Point & p,
			 const Parameters &,
			 const std::string &,
			 const std::string & )
{
  const Real x = p(0);
  // const Real y = p(1);
  // const Real z = p(2);

  // return -2*(x+1)*(1-z*z) - 2*(x+1)*(1-y*y);
  return -2*(x+1);
}

// We now define the matrix assembly function for the
// Laplace system.  We need to first compute element volume
// matrices, and then take into account the boundary
// conditions and the flux integrals, which will be handled
// via an interior penalty method.
void assemble_ellipticdg(EquationSystems & es,
                         const std::string & libmesh_dbg_var(system_name))
{
  libMesh::out << " assembling elliptic dg system... ";
  libMesh::out.flush();

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "EllipticDG");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();
  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & ellipticdg_system = es.get_system<LinearImplicitSystem> ("EllipticDG");
  // Get some parameters that we need during assembly
  const Real penalty = es.parameters.get<Real> ("penalty");
  std::string refinement_type = es.parameters.get<std::string> ("refinement");

  // A reference to the DofMap object for this system. The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap & dof_map = ellipticdg_system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = ellipticdg_system.variable_type(0);

  // Quadrature rules for numerical integration.
#ifdef QORDER
  QGauss qrule(dim, QORDER);
  QGauss qface(dim-1, QORDER);
#else
  QGauss qrule(dim, fe_type.default_quadrature_order());
  QGauss qface(dim-1, fe_type.default_quadrature_order());
#endif

  FE_Wrapper fe( FEBase::build(dim, fe_type), &qrule );
  FE_Wrapper fe_elem_face( FEBase::build(dim, fe_type), &qface );
  FE_Wrapper fe_neighbor_face( FEBase::build(dim, fe_type), &qface );

  // Local matrix and right-hand-side vector
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // Data structures to contain the element and neighbor boundary matrix
  // contribution. This matrices will do the coupling between the dofs of
  // the element and those of his neighbors.
  // Ken: matrix coupling elem and neighbor dofs
  DenseMatrix<Number> Kne;
  DenseMatrix<Number> Ken;
  DenseMatrix<Number> Kee;
  DenseMatrix<Number> Knn;

  // Now we will loop over all the elements in the mesh.  We will
  // compute first the element interior matrix and right-hand-side contribution
  // and then the element and neighbors boundary matrix contributions.
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // Init local dofs
      fe.init_dofs(elem, dof_map);
      fe.fe->reinit(elem);
      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.
      Ke.resize(fe.n_dofs(), fe.n_dofs());
      Fe.resize(fe.n_dofs());

      // Now we will build the element interior matrix.  This involves
      // a double loop to integrate the test functions (i) against
      // the trial functions (j).
      for (unsigned int qp=0; qp<fe.n_quad_points(); qp++)
        for (unsigned int i=0; i<fe.n_dofs(); i++)
          for (unsigned int j=0; j<fe.n_dofs(); j++)
	    {
	      Ke(i,j) += fe.JxW[qp]*(fe.dphi[i][qp]*fe.dphi[j][qp]);
	    }

      // Now we will build the RHS force
      for (unsigned int qp=0; qp<fe.n_quad_points(); qp++)
        for (unsigned int i=0; i<fe.n_dofs(); i++)
	  {
	    Fe(i) += -fe.JxW[qp]*exact_laplacian(fe.qrule_points[qp],
						 es.parameters, "null", "void")*fe.phi[i][qp];
	  }

      // Now we address boundary conditions.
      // We consider Dirichlet bc imposed via the interior penalty method <<
      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      for (auto side : elem->side_index_range())
        {
	  if (elem->neighbor_ptr(side) == libmesh_nullptr)
            {
              fe_elem_face.fe->reinit(elem, side);

              // Compute h element dimension (used for interior penalty parameter)
              std::unique_ptr<const Elem> elem_side (elem->build_side_ptr(side));
              const unsigned int elem_b_order
  		= static_cast<unsigned int> (fe_elem_face.fe->get_order());
              const double h_elem
  		= elem->volume()/elem_side->volume() * 1./pow(elem_b_order, 2.);

  	      // Integrate in element face
  	      auto const n_dofs = fe.n_dofs();
  	      auto const n_face_points = fe_elem_face.n_quad_points();
              for (unsigned int qp=0; qp<n_face_points; qp++)
                {
		  auto const& weight = fe_elem_face.JxW[qp];
  		  auto const& normal = fe_elem_face.qrule_normals[qp];
                  const Number bc_value = exact_solution(fe_elem_face.qrule_points[qp],
							 es.parameters, "null", "void");
		  for (unsigned int i=0; i<n_dofs; i++)
                    {
  		      auto const& v  = fe_elem_face.phi[i][qp];
  		      auto const& dv = fe_elem_face.dphi[i][qp];
                      // Matrix contribution
                      for (unsigned int j=0; j<n_dofs; j++)
                        {
  			  auto const& u  = fe_elem_face.phi[j][qp];
  			  auto const& du = fe_elem_face.dphi[j][qp];
  			  // stability
                          Ke(i,j) += weight * penalty/h_elem * u * v;
                          // consistency
                          Ke(i,j) -= weight * (v*(du*normal) + u*(dv*normal));
                        }
                      // RHS contributions
                      //   + stability
                      Fe(i) += weight * bc_value * penalty/h_elem * v;
                      //   + consistency
                      Fe(i) -= weight * dv * (bc_value*normal);
                    }
                }
            }

          // If the element is not on a boundary of the domain
          // we loop over his neighbors to compute the element
          // and neighbor boundary matrix contributions
          else
            {
              // Store a pointer to the neighbor we are currently
              // working on.
              const Elem * neighbor = elem->neighbor_ptr(side);

              // Get the global id of the element and the neighbor
              const unsigned int elem_id = elem->id();
              const unsigned int neighbor_id = neighbor->id();

              // If the neighbor has the same h level and is active
              // perform integration only if our global id is bigger than our neighbor id.
              // We don't want to compute twice the same contributions.
              // If the neighbor has a different h level perform integration
              // only if the neighbor is at a lower level.
              if ((neighbor->active() &&
                   (neighbor->level() == elem->level()) &&
                   (elem_id < neighbor_id)) ||
                  (neighbor->level() < elem->level()))
                {
                  // Pointer to the element side
                  std::unique_ptr<const Elem> elem_side (elem->build_side_ptr(side));

                  // h dimension to compute the interior penalty penalty parameter
                  const unsigned int elem_b_order = static_cast<unsigned int>(fe_elem_face.fe->get_order());
                  const unsigned int neighbor_b_order = static_cast<unsigned int>(fe_neighbor_face.fe->get_order());
                  const double side_order = (elem_b_order + neighbor_b_order)/2.;
                  const double h_elem = (elem->volume()/elem_side->volume()) * 1./pow(side_order,2.);

                  // The quadrature point locations on the neighbor side
                  std::vector<Point> qface_neighbor_point;

                  // The quadrature point locations on the element side
                  std::vector<Point > qface_point;

                  // Reinitialize shape functions on the element side
                  fe_elem_face.fe->reinit(elem, side);

                  // Get the physical locations of the element quadrature points
                  qface_point = fe_elem_face.fe->get_xyz();

                  // Find their locations on the neighbor (save in qface_neighbor_point)
                  unsigned int side_neighbor = neighbor->which_neighbor_am_i(elem);
                  if (refinement_type == "p")
                    fe_neighbor_face.fe->side_map (neighbor,
                                                elem_side.get(),
                                                side_neighbor,
                                                fe_elem_face.qrule->get_points(),
                                                qface_neighbor_point);
                  else
                    FEInterface::inverse_map (elem->dim(),
                                              fe.fe->get_fe_type(),
                                              neighbor,
                                              qface_point,
                                              qface_neighbor_point);

                  // Calculate the neighbor element shape functions at those locations
                  fe_neighbor_face.fe->reinit(neighbor, &qface_neighbor_point);

                  // Get the degree of freedom indices for the
                  // neighbor.  These define where in the global
                  // matrix this neighbor will contribute to.
  		  fe_neighbor_face.init_dofs(neighbor, dof_map);

		  DG_FaceCoupling face_coupling(fe_elem_face, fe_neighbor_face);
		  SIP_BilinearForm a_sip(face_coupling, penalty, h_elem);
		  FaceIntegrator<SIP_BilinearForm> a_sip_integrator(a_sip);
		  a_sip_integrator.integrate();

                  // The element and neighbor boundary matrix are now built
                  // for this side.  Add them to the global matrix
                  // The SparseMatrix::add_matrix() members do this for us.
		  a_sip_integrator.save_to_system(ellipticdg_system);
                }
            }
        }
      // The element interior matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      ellipticdg_system.matrix->add_matrix(Ke, fe.dof_indices);
      ellipticdg_system.rhs->add_vector(Fe, fe.dof_indices);
    }

  libMesh::out << "done" << std::endl;
}



int main (int argc, char** argv)
{
  LibMeshInit init(argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  //Parse the input file
  GetPot input_file("dg_laplace_err_order.in");

  //Read in parameters from the input file
  const unsigned int uniform_refinement_steps  = input_file("uniform_h_r_steps", 3);
  const Real refine_fraction                   = input_file("refine_fraction", 0.5);
  const Real coarsen_fraction                  = input_file("coarsen_fraction", 0.);
  const unsigned int max_h_level               = input_file("max_h_level", 10);
  const std::string refinement_type            = input_file("refinement_type","p");
  Order p_order                                = static_cast<Order>(input_file("p_order", 1));
  const std::string element_type               = input_file("element_type", "tensor");
  const Real penalty                           = input_file("ip_penalty", 10.);
  const bool singularity                       = input_file("singularity", true);
  const unsigned int dim                       = input_file("dimension", 3);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");


  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  if (dim == 1)
    MeshTools::Generation::build_line(mesh, 1, -1., 0.);
  else if (dim == 2)
    MeshTools::Generation::build_square(mesh, 2, 2,
					-1, 1,
					-1, 1,
					QUAD9);
  else
    MeshTools::Generation::build_cube(mesh, 2, 2, 2,
				      -1, 1,
				      -1, 1,
				      -1, 1,
				      PRISM6
				      );

  // Use triangles if the config file says so
  if (element_type == "simplex")
    MeshTools::Modification::all_tri(mesh);

  // Mesh Refinement object
  MeshRefinement mesh_refinement(mesh);
  mesh_refinement.refine_fraction() = refine_fraction;
  mesh_refinement.coarsen_fraction() = coarsen_fraction;
  mesh_refinement.max_h_level() = max_h_level;

  // // Do uniform refinement
  // for (unsigned int rstep=0; rstep<uniform_refinement_steps; rstep++)
  //   mesh_refinement.uniformly_refine(1);

  // Crate an equation system object
  EquationSystems equation_system (mesh);

  // Set parameters for the equation system and the solver
  equation_system.parameters.set<Real>("linear solver tolerance") = TOLERANCE * TOLERANCE;
  equation_system.parameters.set<unsigned int>("linear solver maximum iterations") = 1000;
  equation_system.parameters.set<Real>("penalty") = penalty;
  equation_system.parameters.set<bool>("singularity") = singularity;
  equation_system.parameters.set<std::string>("refinement") = refinement_type;

  // Create a system named ellipticdg
  LinearImplicitSystem & ellipticdg_system = equation_system.add_system<LinearImplicitSystem> ("EllipticDG");

  // Add a variable "u" to "ellipticdg" using the p_order specified in the config file
  if (on_command_line("element_type"))
    {
      std::string fe_str =
        command_line_value(std::string("element_type"),
                           std::string("MONOMIAL"));

      if (fe_str != "MONOMIAL" || fe_str != "XYZ")
        libmesh_error_msg("Error: This example must be run with MONOMIAL or XYZ element types.");

      ellipticdg_system.add_variable ("u", p_order, Utility::string_to_enum<FEFamily>(fe_str));
    }
  else
    ellipticdg_system.add_variable ("u", p_order, MONOMIAL);

  // Give the system a pointer to the matrix assembly function
  ellipticdg_system.attach_assemble_function (assemble_ellipticdg);

  // Initialize the data structures for the equation system
  equation_system.init();


  // Construct ExactSolution object and attach solution functions
  ExactSolution exact_sol(equation_system);
  exact_sol.attach_exact_value(exact_solution);
  exact_sol.attach_exact_deriv(exact_derivative);
  // exact_sol.extra_quadrature_order(1);

  // A refinement loop.
  for (unsigned int rstep=0; rstep<uniform_refinement_steps; ++rstep)
    {
      libMesh::out << "\nBeginning Solve " << rstep << std::endl;
      libMesh::out << "Number of elements: " << mesh.n_elem() << std::endl;

      // Solve the system
      ellipticdg_system.solve();

      libMesh::out << "System has: "
                   << equation_system.n_active_dofs()
                   << " degrees of freedom."
                   << std::endl;

      libMesh::out << "Linear solver converged at step: "
                   << ellipticdg_system.n_linear_iterations()
                   << ", final residual: "
                   << ellipticdg_system.final_linear_residual()
                   << std::endl;

      // Compute the error
      exact_sol.compute_error("EllipticDG", "u");

      // Print out the error values
      Real l2_error = exact_sol.l2_error("EllipticDG", "u");
      Real h1_error = exact_sol.h1_error("EllipticDG", "u");
      libMesh::out << "L2-Error is: " << l2_error << std::endl;
      libMesh::out << "H1-Error is: " << h1_error << std::endl;

      // Compute error order
      Real old_l2_error, old_h1_error, l2_error_order, h1_error_order;
      if(rstep>0) {
	l2_error_order = log(old_l2_error/l2_error)/log(2);
	h1_error_order = log(old_h1_error/h1_error)/log(2);
	libMesh::out << "  L2-Error-Order is: " << l2_error_order << std::endl;
	libMesh::out << "  H1-Error-Order is: " << h1_error_order << std::endl;
      }
      old_l2_error = l2_error;
      old_h1_error = h1_error;

      // // Take the error in error and decide which elements will be coarsened or refined
      // mesh_refinement.flag_elements_by_error_fraction(error);
      // if (refinement_type == "p")
      //   mesh_refinement.switch_h_to_p_refinement();
      // if (refinement_type == "hp")
      //   mesh_refinement.add_p_to_h_refinement();

      // Refine and coarsen the flagged elements
      mesh_refinement.uniformly_refine(1);
      equation_system.reinit();
    }

  // Write out the solution
  // After solving the system write the solution
  // to a ExodusII-formatted plot file.
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO (mesh).write_discontinuous_exodusII("output.e", equation_system);
#endif

#endif // #ifndef LIBMESH_ENABLE_AMR

  // All done.
  return 0;
}
