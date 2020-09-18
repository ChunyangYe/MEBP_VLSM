#include <iostream>
#include <string>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/preprocessor/seq/for_each.hpp>


#if defined(SOLVER_BACKEND_VEXCL)
#  include <amgcl/value_type/static_matrix.hpp>
#  include <amgcl/adapter/block_matrix.hpp>
#  include <amgcl/backend/vexcl.hpp>
#  include <amgcl/backend/vexcl_static_matrix.hpp>
typedef amgcl::backend::vexcl<double> Backend;
#elif defined(SOLVER_BACKEND_VIENNACL)
#  include <amgcl/backend/viennacl.hpp>
typedef amgcl::backend::viennacl< viennacl::compressed_matrix<double> > Backend;
#elif defined(SOLVER_BACKEND_CUDA)
#  include <amgcl/backend/cuda.hpp>
#  include <amgcl/relaxation/cusparse_ilu0.hpp>
typedef amgcl::backend::cuda<double> Backend;
#elif defined(SOLVER_BACKEND_EIGEN)
#  include <amgcl/backend/eigen.hpp>
typedef amgcl::backend::eigen<double> Backend;
#elif defined(SOLVER_BACKEND_BLAZE)
#  include <amgcl/backend/blaze.hpp>
typedef amgcl::backend::blaze<double> Backend;
#else
#  ifndef SOLVER_BACKEND_BUILTIN
#    define SOLVER_BACKEND_BUILTIN
#  endif
#  include <amgcl/backend/builtin.hpp>
#  include <amgcl/value_type/static_matrix.hpp>
#  include <amgcl/adapter/block_matrix.hpp>
typedef amgcl::backend::builtin<double> Backend;
#endif

#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/coarsening/rigid_body_modes.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/preconditioner/runtime.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/reorder.hpp>
#include <amgcl/io/mm.hpp>
#include <amgcl/io/binary.hpp>

#include <amgcl/profiler.hpp>


#ifndef AMGCL_BLOCK_SIZES
#  define AMGCL_BLOCK_SIZES (3)(4)
#endif

namespace amgcl { profiler<> prof; }
using amgcl::prof;
using amgcl::precondition;

#ifdef SOLVER_BACKEND_BUILTIN
//---------------------------------------------------------------------------
template <int B>
std::tuple<size_t, double> block_solve(
	const boost::property_tree::ptree &prm,
	size_t rows,
	std::vector<ptrdiff_t> const &ptr,
	std::vector<ptrdiff_t> const &col,
	std::vector<double>    const &val,
	std::vector<double>    const &rhs,
	std::vector<double>          &x,
	bool reorder
)
{
	typedef amgcl::static_matrix<double, B, B> value_type;
	typedef amgcl::static_matrix<double, B, 1> rhs_type;
	typedef amgcl::backend::builtin<value_type> BBackend;

	typedef amgcl::make_solver<
		amgcl::runtime::preconditioner<BBackend>,
		amgcl::runtime::solver::wrapper<BBackend>
	> Solver;

	auto A = amgcl::adapter::block_matrix<value_type>(std::tie(rows, ptr, col, val));

	std::tuple<size_t, double> info;

	if (reorder) {
		prof.tic("reorder");
		amgcl::adapter::reorder<> perm(A);
		prof.toc("reorder");

		prof.tic("setup");
		Solver solve(perm(A), prm);
		prof.toc("setup");

		std::cout << solve << std::endl;

		rhs_type const * fptr = reinterpret_cast<rhs_type const *>(&rhs[0]);
		rhs_type       * xptr = reinterpret_cast<rhs_type       *>(&x[0]);

		amgcl::backend::numa_vector<rhs_type> F(perm(amgcl::make_iterator_range(fptr, fptr + rows / B)));
		amgcl::backend::numa_vector<rhs_type> X(perm(amgcl::make_iterator_range(xptr, xptr + rows / B)));

		prof.tic("solve");
		info = solve(F, X);
		prof.toc("solve");

		perm.inverse(X, xptr);
	}
	else {
		prof.tic("setup");
		Solver solve(A, prm);
		prof.toc("setup");

		std::cout << solve << std::endl;

		rhs_type const * fptr = reinterpret_cast<rhs_type const *>(&rhs[0]);
		rhs_type       * xptr = reinterpret_cast<rhs_type       *>(&x[0]);

		amgcl::backend::numa_vector<rhs_type> F(fptr, fptr + rows / B);
		amgcl::backend::numa_vector<rhs_type> X(xptr, xptr + rows / B);

		prof.tic("solve");
		info = solve(F, X);
		prof.toc("solve");

		std::copy(X.data(), X.data() + X.size(), xptr);
	}

	return info;
}
#endif


//---------------------------------------------------------------------------
std::tuple<size_t, double> scalar_solve(
	const boost::property_tree::ptree &prm,
	size_t rows,
	std::vector<ptrdiff_t> const &ptr,
	std::vector<ptrdiff_t> const &col,
	std::vector<double>    const &val,
	std::vector<double>    const &rhs,
	std::vector<double>          &x,
	bool reorder
)
{
	Backend::params bprm;

#if defined(SOLVER_BACKEND_VEXCL)
	vex::Context ctx(vex::Filter::Env);
	std::cout << ctx << std::endl;
	bprm.q = ctx;
#elif defined(SOLVER_BACKEND_VIENNACL)
	std::cout
		<< viennacl::ocl::current_device().name()
		<< " (" << viennacl::ocl::current_device().vendor() << ")\n\n";
#elif defined(SOLVER_BACKEND_CUDA)
	cusparseCreate(&bprm.cusparse_handle);
	{
		int dev;
		cudaGetDevice(&dev);

		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, dev);
		std::cout << prop.name << std::endl << std::endl;
	}
#endif

	typedef amgcl::make_solver<
		amgcl::runtime::preconditioner<Backend>,
		amgcl::runtime::solver::wrapper<Backend>
	> Solver;

	std::tuple<size_t, double> info;

	if (reorder) {
		prof.tic("reorder");
		amgcl::adapter::reorder<> perm(std::tie(rows, ptr, col, val));
		prof.toc("reorder");

		prof.tic("setup");
		Solver solve(perm(std::tie(rows, ptr, col, val)), prm, bprm);
		prof.toc("setup");

		std::cout << solve << std::endl;

		std::vector<double> tmp(rows);

		perm.forward(rhs, tmp);
		auto f_b = Backend::copy_vector(tmp, bprm);

		perm.forward(x, tmp);
		auto x_b = Backend::copy_vector(tmp, bprm);

		prof.tic("solve");
		info = solve(*f_b, *x_b);
		prof.toc("solve");

#if defined(SOLVER_BACKEND_VEXCL)
		vex::copy(*x_b, tmp);
#elif defined(SOLVER_BACKEND_VIENNACL)
		viennacl::fast_copy(*x_b, tmp);
#elif defined(SOLVER_BACKEND_CUDA)
		thrust::copy(x_b->begin(), x_b->end(), tmp.begin());
#else
		std::copy(&(*x_b)[0], &(*x_b)[0] + rows, &tmp[0]);
#endif

		perm.inverse(tmp, x);
	}
	else {
		prof.tic("setup");
		Solver solve(std::tie(rows, ptr, col, val), prm, bprm);
		prof.toc("setup");

		std::cout << solve << std::endl;

		auto f_b = Backend::copy_vector(rhs, bprm);
		auto x_b = Backend::copy_vector(x, bprm);

		prof.tic("solve");
		info = solve(*f_b, *x_b);
		prof.toc("solve");

#if defined(SOLVER_BACKEND_VEXCL)
		vex::copy(*x_b, x);
#elif defined(SOLVER_BACKEND_VIENNACL)
		viennacl::fast_copy(*x_b, x);
#elif defined(SOLVER_BACKEND_CUDA)
		thrust::copy(x_b->begin(), x_b->end(), x.begin());
#else
		std::copy(&(*x_b)[0], &(*x_b)[0] + rows, &x[0]);
#endif
	}

	return info;
}

#define AMGCL_CALL_BLOCK_SOLVER(z, data, B)                                    \
  case B:                                                                      \
    return block_solve<B>(prm, rows, ptr, col, val, rhs, x, reorder);

//---------------------------------------------------------------------------
std::tuple<size_t, double> solve(
	const boost::property_tree::ptree &prm,
	size_t rows,
	std::vector<ptrdiff_t> const &ptr,
	std::vector<ptrdiff_t> const &col,
	std::vector<double>    const &val,
	std::vector<double>    const &rhs,
	std::vector<double>          &x,
	int block_size,
	bool reorder
)
{
	switch (block_size) {
	case 1:
		return scalar_solve(prm, rows, ptr, col, val, rhs, x, reorder);
#if defined(SOLVER_BACKEND_BUILTIN) || defined(SOLVER_BACKEND_VEXCL)
		BOOST_PP_SEQ_FOR_EACH(AMGCL_CALL_BLOCK_SOLVER, ~, AMGCL_BLOCK_SIZES)
#endif
	default:
		precondition(false, "Unsupported block size");
		return std::make_tuple(0, 0.0);
	}
}

