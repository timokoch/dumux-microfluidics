// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Linear
 * \brief Preconditioned iterative solver for the incompressible Stokes problem
 */
#ifndef DUMUX_INCOMPRESSIBLE_STOKES_SOLVER_HH
#define DUMUX_INCOMPRESSIBLE_STOKES_SOLVER_HH

#include <type_traits>
#include <memory>
#include <tuple>

#include <dune/common/parametertree.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>

#include <dumux/common/math.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/preconditioners.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/linearsolvertraits.hh>

#include <dumux/assembly/coloring.hh>
#include <dumux/parallel/parallel_for.hh>
#include <dumux/discretization/extrusion.hh>

#include "symmetrizedirichlet.hh"

namespace Dumux {

namespace Detail {

// compute coloring for parallel sweep
template<class M>
void computeColorsForMatrixSweep_(const M& matrix, std::deque<std::vector<std::size_t>>& coloredIndices)
{
    // allocate some temporary memory
    std::vector<int> colors(matrix.N(), -1);
    std::vector<int> neighborColors; neighborColors.reserve(30);
    std::vector<bool> colorUsed; colorUsed.reserve(30);

    // a row that has my row index in their column set cannot have the same color
    // TODO this assumes a symmetric matrix pattern
    for (std::size_t i = 0; i < colors.size(); ++i)
    {
        neighborColors.clear();
        auto& row = matrix[i];
        const auto endColIt = row.end();
        for (auto colIt = row.begin(); colIt != endColIt; ++colIt)
            neighborColors.push_back(colors[colIt.index()]);

        const auto c = Detail::smallestAvailableColor(neighborColors, colorUsed);
        colors[i] = c;

        // add element to the set of elements with the same color
        if (c < coloredIndices.size())
            coloredIndices[c].push_back(i);
        else
            coloredIndices.emplace_back(1, i);
    }
}

// parallel SOR kernel (relaxed Gauss-Seidel)
template<bool forward, int l, class M, class X, class Y, class K>
void parallelBlockSOR_(const M& A, X& update, const Y& b, const K& relaxationFactor,
                       const std::deque<std::vector<std::size_t>>& coloredIndices)
{
    for (int color = 0; color < coloredIndices.size(); ++color)
    {
        const auto c = forward ? color : coloredIndices.size()-1-color;
        const auto& rowIndices = coloredIndices[c];
        Dumux::parallelFor(rowIndices.size(), [&](const std::size_t k)
        {
            const auto i = rowIndices[k];
            auto& row = A[i];
            auto v = update[i];
            auto rhs = b[i];
            const auto endColIt = row.end();
            auto colIt = row.begin();

            for (; colIt.index()<i; ++colIt)
                colIt->mmv(update[colIt.index()], rhs);
            const auto diag = colIt;
            for (; colIt != endColIt; ++colIt)
                colIt->mmv(update[colIt.index()], rhs);

            if constexpr (forward)
                Dune::algmeta_itsteps<l-1,typename M::block_type>::bsorf(
                    *diag, v, rhs, relaxationFactor
                );
            else
                Dune::algmeta_itsteps<l-1,typename M::block_type>::bsorb(
                    *diag, v, rhs, relaxationFactor
                );

            update[i].axpy(relaxationFactor, v);
        });
    }
}

} // end namespace Detail

/*! \brief Multihreaded SOR preconditioner
 *
 *  \tparam M The matrix type to operate on
 *  \tparam X Type of the update
 *  \tparam Y Type of the defect
 *  \tparam l The block level to invert. Default is 1
 */
template<class M, class X, class Y, int l=1>
class ParMTSOR : public Dune::Preconditioner<X,Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    typedef M matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;
    //! \brief scalar type underlying the field_type
    typedef Dune::Simd::Scalar<field_type> scalar_field_type;
    //! \brief real scalar type underlying the field_type
    typedef typename Dune::FieldTraits<scalar_field_type>::real_type real_field_type;

    //! \brief Constructor.
    ParMTSOR (const M& A, int n, real_field_type w)
    : _A_(A), numSteps_(n), relaxationFactor_(w)
    {
        Dune::CheckIfDiagonalPresent<M,l>::check(_A_);
        Detail::computeColorsForMatrixSweep_(_A_, colors_);
    }

    ParMTSOR (const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& A, const Dune::ParameterTree& configuration)
    : ParMTSOR(A->getmat(), configuration)
    {}

    ParMTSOR (const M& A, const Dune::ParameterTree& configuration)
    : ParMTSOR(A, configuration.get<int>("iterations",1), configuration.get<real_field_type>("relaxation",1.0))
    {}

    void pre (X&, Y&) override {}

    void apply (X& v, const Y& d) override
    {
        this->template apply<true>(v,d);
    }

    template<bool forward>
    void apply(X& v, const Y& d)
    {
        for (int i=0; i<numSteps_; i++)
            Detail::parallelBlockSOR_<forward, l>(_A_, v, d, relaxationFactor_, colors_);
    }

    void post (X&) override {}

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

private:
    const M& _A_;
    int numSteps_;
    real_field_type relaxationFactor_;

    //! for each color a vector of row indices that can be dealt with in parallel
    std::deque<std::vector<std::size_t>> colors_;
};

/*! \brief Multithreaded SSOR preconditioner
 *
 *  \tparam M The matrix type to operate on
 *  \tparam X Type of the update
 *  \tparam Y Type of the defect
 *  \tparam l The block level to invert. Default is 1
 */
template<class M, class X, class Y, int l=1>
class ParMTSSOR : public Dune::Preconditioner<X,Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    typedef M matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;
    //! \brief scalar type underlying the field_type
    typedef Dune::Simd::Scalar<field_type> scalar_field_type;
    //! \brief real scalar type underlying the field_type
    typedef typename Dune::FieldTraits<scalar_field_type>::real_type real_field_type;

    //! \brief Constructor.
    ParMTSSOR (const M& A, int n, real_field_type w)
    : _A_(A), numSteps_(n), relaxationFactor_(w)
    {
        Dune::CheckIfDiagonalPresent<M,l>::check(_A_);
        Detail::computeColorsForMatrixSweep_(_A_, colors_);
    }

    ParMTSSOR (const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& A, const Dune::ParameterTree& configuration)
    : ParMTSSOR(A->getmat(), configuration)
    {}

    ParMTSSOR (const M& A, const Dune::ParameterTree& configuration)
    : ParMTSSOR(A, configuration.get<int>("iterations",1), configuration.get<real_field_type>("relaxation",1.0))
    {}

    void pre (X&, Y&) override {}

    void apply (X& v, const Y& d) override
    {
        for (int i=0; i<numSteps_; i++)
        {
            Detail::parallelBlockSOR_<true, l>(_A_, v, d, relaxationFactor_, colors_);
            Detail::parallelBlockSOR_<false, l>(_A_, v, d, relaxationFactor_, colors_);
        }
    }

    void post (X&) override {}

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

private:
    const M& _A_;
    int numSteps_;
    real_field_type relaxationFactor_;

    //! for each color a vector of row indices that can be dealt with in parallel
    std::deque<std::vector<std::size_t>> colors_;
};

} // end namespace Dumux

namespace Dune::Amg {

// make it possible to use Dumux::ParMTSOR as AMG smoother
template<class M, class X, class Y, int l>
struct ConstructionTraits<Dumux::ParMTSOR<M,X,Y,l> >
{
    using Arguments = DefaultConstructionArgs<SeqSOR<M,X,Y,l> >;
    static inline std::shared_ptr<Dumux::ParMTSOR<M,X,Y,l>> construct(Arguments& args)
    {
        return std::make_shared<Dumux::ParMTSOR<M,X,Y,l>>(
            args.getMatrix(), args.getArgs().iterations, args.getArgs().relaxationFactor
        );
    }
};

template<class M, class X, class Y, int l>
struct SmootherApplier<Dumux::ParMTSOR<M,X,Y,l> >
{
    typedef Dumux::ParMTSOR<M,X,Y,l> Smoother;
    typedef typename Smoother::range_type Range;
    typedef typename Smoother::domain_type Domain;

    static void preSmooth(Smoother& smoother, Domain& v, Range& d)
    { smoother.template apply<true>(v,d); }

    static void postSmooth(Smoother& smoother, Domain& v, Range& d)
    { smoother.template apply<false>(v,d); }
};

// make it possible to use Dumux::ParMTSSOR as AMG smoother
template<class M, class X, class Y, int l>
struct ConstructionTraits<Dumux::ParMTSSOR<M,X,Y,l> >
{
    using Arguments = DefaultConstructionArgs<SeqSSOR<M,X,Y,l> >;
    static inline std::shared_ptr<Dumux::ParMTSSOR<M,X,Y,l>> construct(Arguments& args)
    {
        return std::make_shared<Dumux::ParMTSSOR<M,X,Y,l>>(
            args.getMatrix(), args.getArgs().iterations, args.getArgs().relaxationFactor
        );
    }
};

} // end namespace Dune::Amg

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A preconditioner based on the preconditioned Uzawa algorithm for saddle-point problems of the form
 * \f$
 \begin{pmatrix} A & B \\ C & 0 \end{pmatrix}
 \begin{pmatrix} u \\ p \end{pmatrix} =
 \begin{pmatrix} f \\ g \end{pmatrix}
 * \f$
 *
 * This preconditioner is especially suited for solving the incompressible Stokes equations.
 *
 * \tparam M Type of the matrix.
 * \tparam X Type of the update.
 * \tparam Y Type of the defect.
 * \tparam l Preconditioner block level (for compatibility reasons, unused).
 */
template<class M, class X, class Y, int l = 1>
class IncompressibleStokesPreconditioner : public Dune::Preconditioner<X,Y>
{
    static_assert(Dumux::isMultiTypeBlockMatrix<M>::value && M::M() == 2 && M::N() == 2, "Expects a 2x2 MultiTypeBlockMatrix.");
    static_assert(l== 1, "IncompressibleStokesPreconditioner expects a block level of 1.");

    using A = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_0][Dune::Indices::_0])>;
    using U = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_0])>;

    using P = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_1][Dune::Indices::_1])>;
    using V = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_1])>;

public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = M;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;
    //! \brief Scalar type underlying the field_type.
    using scalar_field_type = Dune::Simd::Scalar<field_type>;
    //! \brief the type of the pressure operator
    using PressureLinearOperator = Dune::MatrixAdapter<P,V,V>;

    IncompressibleStokesPreconditioner(
        const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& op,
        const std::shared_ptr<const Dune::AssembledLinearOperator<P,V,V>>& pop,
        const Dune::ParameterTree& params
    )
    : matrix_(op->getmat())
    , pmatrix_(pop->getmat())
    , verbosity_(params.get<int>("verbosity"))
    , paramGroup_(params.get<std::string>("ParameterGroup"))
    {
        initPreconditioner_(params);
    }

    /*!
     * \brief Prepare the preconditioner.
     */
    virtual void pre(X& x, Y& b) {}

    /*!
     * \brief Apply the preconditioner
     *
     * \param update The update to be computed.
     * \param currentDefect The current defect.
     *
     * The currentDefect has be be in a consistent representation,
     * Definition 2.3 Blatt and Bastian (2009) https://doi.org/10.1504/IJCSE.2008.021112
     * The update is initially zero. At exit the update has to be
     * in a consistent representation. This usually requires communication.
     */
    virtual void apply(X& update, const Y& currentDefect)
    {
        using namespace Dune::Indices;

        [[maybe_unused]] const auto& B = matrix_[_0][_1];
        [[maybe_unused]] const auto& C = matrix_[_1][_0];
        const auto& f = currentDefect[_0];
        const auto& g = currentDefect[_1];
        auto& u = update[_0];
        auto& p = update[_1];

        // The preconditioned Uzawa iteration
        //
        // in comparison to the classical SeqUzawa, we precondition the
        // pressure update. Also, here we assume D = 0 and that the initial
        // update is zero.
        // Depending on the problem, the pressure preconditioner might be an
        // inverse weighted mass matrix or a (compatible) Laplacian, or
        // a weighted combination, see https://doi.org/10.1002/fld.1650080802

        // u_k+1 = u_k + Q_A^−1*f,
        auto uRhs = f;
        applyPreconditionerForA_(u, uRhs);

        // p_k+1 = p_k + Q_P^-1*(g - C*u_k+1)
        auto pRhs = g;
        C.mmv(u, pRhs);
        applyPreconditionerForP_(p, pRhs);

        // full Schur complement style update (symmetric)
        // requires another application of Q_A^−1
        // u_k+2 = u_k+1 - Q_A^−1*(Bp_k+1)
        // uRhs = 0;
        // B.mv(p, uRhs);
        // auto uIncrement = uRhs; uIncrement = 0;
        // applyPreconditionerForA_(uIncrement, uRhs);
        // u -= uIncrement; // uIncrement = -(u_k+2 - u_k+1)

        // instead we only do a projection
        // project velocity into divergence-free space
        // u_k+2 = u_k+1 + Bp_k+1
        // auto uHalf = u;
        // when using this in parallel:
        // uHalf is already consistent
        // we have to make u consistent after applying B
        // i.e. sum up over the border entities
        // u = 0.0;
        // B.umv(p, u);
        // comms_[0]->addOwnerCopyToOwnerCopy(u, u);
        // u += uHalf;

        // // update pressure
        // // p_k+2 = p_k+1 + Cu_k+1 - g
        // C.umv(uHalf, p);
        // p -= g;
    }

    /*!
     * \brief Clean up.
     */
    virtual void post(X& x) {}

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::sequential;
        //return Dune::SolverCategory::overlapping;
    }

private:
    void initPreconditioner_(const Dune::ParameterTree& params)
    {
        using namespace Dune::Indices;

        if (getParamFromGroup<bool>(paramGroup_, "LinearSolver.DirectSolverForVelocity", false))
        {
            directSolver_ = std::make_shared<Dune::UMFPack<A>>(matrix_[_0][_0], verbosity_);
            using Wrap = Dune::InverseOperator2Preconditioner<Dune::InverseOperator<U, U>>;
            preconditionerForA_ = std::make_shared<Wrap>(*directSolver_);
        }
        else
        {
            using VelLinearOperator = Dune::MatrixAdapter<A, U, U>;
            auto lopV = std::make_shared<VelLinearOperator>(matrix_[_0][_0]);
            // using AMGForA = Dune::Amg::AMG<VelLinearOperator, U, Dumux::ParMTSOR<A, U, U>, Comm>;
            // using VelLinearOperator = Dune::NonoverlappingSchwarzOperator<A, U, U, Comm>;
            // auto lopV = std::make_shared<VelLinearOperator>(matrix_[_0][_0], *comms_[0]);
            // auto creator = Dune::AMGCreator();
            // preconditionerForA_ = creator.makeAMG(lopV, "ssor", params);
            preconditionerForA_ = std::make_shared<Dune::Amg::AMG<VelLinearOperator, U, Dumux::ParMTSOR<A,U,U>>>(lopV, params);
        }

        using PressLinearOperator = Dune::MatrixAdapter<P, V, V>;
        auto lopP = std::make_shared<PressLinearOperator>(pmatrix_);
        using PressJacobi = Dune::SeqJac<P, V, V>;
        // using ParPressJacobi = Dune::BlockPreconditioner<V, V, Comm, PressJacobi>;
        // auto seqPre = std::make_shared<PressJacobi>(lopP, params);
        // preconditionerForP_ = std::make_shared<ParPressJacobi>(seqPre, *comms_[1]);
        preconditionerForP_ = std::make_shared<PressJacobi>(lopP, params);
    }

    template<class Sol, class Rhs>
    void applyPreconditionerForA_(Sol& sol, Rhs& rhs) const
    {
        preconditionerForA_->pre(sol, rhs);
        preconditionerForA_->apply(sol, rhs);
        preconditionerForA_->post(sol);
    }

    template<class Sol, class Rhs>
    void applyPreconditionerForP_(Sol& sol, Rhs& rhs) const
    {
        preconditionerForP_->pre(sol, rhs);
        preconditionerForP_->apply(sol, rhs);
        preconditionerForP_->post(sol);
    }

    //! \brief The matrix we operate on.
    const M& matrix_;
    //! \brief The matrix we operate on.
    const P& pmatrix_;
    //! \brief The verbosity level
    const int verbosity_;

    //std::array<std::shared_ptr<const Comm>, 2> comms_;

    std::shared_ptr<Dune::Preconditioner<U, U>> preconditionerForA_;
    std::shared_ptr<Dune::Preconditioner<V, V>> preconditionerForP_;
    std::shared_ptr<Dune::InverseOperator<U, U>> directSolver_;
    const std::string paramGroup_;
};

/*!
 * \ingroup Linear
 * \brief Preconditioned iterative solver for the incompressible Stokes problem
 * \note Uses IncompressibleStokesPreconditioner as preconditioner
 */
template<class Matrix, class Vector, class VelocityGG, class PressureGG>
class IncompressibleStokesSolver
: public LinearSolver
{
    // "[Scalar and matrix-vector products] are easy to compute in a non-overlapping distributed setting,
    // if matrix and residuals are kept in _additive [unique]_ representation,
    // and iterates and directions are kept in _consistent_ representation." (Sander 2020 Dune)
    // "The consistent representation of a vector is obtained from its additive representation by taking the latter,
    // and for each vertex[/face] in the border partition, add the corresponding entries from the other subdomains" (Sander 2020 Dune)
    using Communication = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
    using Preconditioner = IncompressibleStokesPreconditioner<Matrix, Vector, Vector>;
public:
    /*!
     * \brief Constructor
     * \param vGridGeometry grid geometry of the velocity discretization
     * \param pGridGeometry grid geometry of the pressure discretization
     * \param params a parameter tree passed along to the iterative solvers
     *
     * The parameter tree requires the (constant) viscosity in the key "Component.LiquidDynamicViscosity"
     * and the (constant) density in the key "Component.LiquidDensity"
     */
    IncompressibleStokesSolver(std::shared_ptr<const VelocityGG> vGridGeometry,
                               std::shared_ptr<const PressureGG> pGridGeometry,
                               const Vector& dirichletDofs)
    : LinearSolver()
    , vGridGeometry_(std::move(vGridGeometry))
    , pGridGeometry_(std::move(pGridGeometry))
    , dirichletDofs_(dirichletDofs)
    {
        params_ = LinearSolverParameters<LinearSolverTraits<VelocityGG>>::createParameterTree(this->paramGroup());
        density_ = getParamFromGroup<double>(this->paramGroup(), "Component.LiquidDensity");
        viscosity_ = getParamFromGroup<double>(this->paramGroup(), "Component.LiquidDynamicViscosity");
        weight_ = getParamFromGroup<double>(this->paramGroup(), "LinearSolver.Preconditioner.MassMatrixWeight", 1.0);
        solverType_ = getParamFromGroup<std::string>(this->paramGroup(), "LinearSolver.Type", "gmres");

        // pHelperVelocity_ = std::make_unique<ParallelISTLHelper<LinearSolverTraits<VelocityGG>>>(
        //     vGridGeometry_->gridView(), vGridGeometry_->dofMapper()
        // );
        // pHelperPressure_ = std::make_unique<ParallelISTLHelper<LinearSolverTraits<PressureGG>>>(
        //     pGridGeometry_->gridView(), pGridGeometry_->dofMapper()
        // );

        // commVelocity_ = std::make_shared<Communication>(pHelperVelocity_->gridView().comm(), Dune::SolverCategory::nonoverlapping);
        // pHelperVelocity_->createParallelIndexSet(*commVelocity_);
        // commPressure_ = std::make_shared<Communication>(pHelperPressure_->gridView().comm(), Dune::SolverCategory::overlapping);
        // pHelperPressure_->createParallelIndexSet(*commPressure_);

        // const auto comms = std::array<std::shared_ptr<const Communication>, 2>({ commVelocity_, commPressure_ });
        // scalarProduct_ = std::make_shared<Dumux::ParallelScalarProduct<Vector, Communication, Vector::size()>>(comms);

        scalarProduct_ = std::make_shared<Dune::ScalarProduct<Vector>>();
    }

    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        auto bTmp = b;
        auto ATmp = A;

        return applyIterativeSolver_(ATmp, x, bTmp);
    }

    Scalar norm(const Vector& b) const
    {
        return scalarProduct_->norm(b);
    }

    std::string name() const
    {
        return "Block-preconditioned incompressible Stokes solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    bool applyIterativeSolver_(Matrix& A, Vector& x, Vector& b)
    {
        // make Dirichlet boundary conditions symmetric
        // do this first because its easier if the operator is still purely local
        // otherwise we might sum up wrong contributions at border face dofs adjacent to the domain boundary
        if (getParamFromGroup<bool>(this->paramGroup(), "LinearSolver.SymmetrizeDirichlet", false))
            symmetrizeDirichlet(A, b, dirichletDofs_);

        // // make matrix and right-hand side consistent
        // {
        //     using namespace Dune::Indices;
        //     using M = std::decay_t<decltype(std::declval<Matrix>()[_0][_0])>;
        //     using U = std::decay_t<decltype(std::declval<Vector>()[_0])>;
        //     using PTraits = typename LinearSolverTraits<VelocityGG>::template ParallelNonoverlapping<M, U>;
        //     prepareLinearAlgebraParallel<LinearSolverTraits<VelocityGG>, PTraits>(A[_0][_0], b[_0], *pHelperVelocity_);
        // }

        // make Matrix symmetric on the block-scale
        using namespace Dune::Indices;
        A[_1] *= -1.0/density_;
        b[_1] *= -1.0/density_;

        // const auto comms = std::array<std::shared_ptr<const Communication>, 2>({ commVelocity_, commPressure_ });
        // auto op = std::make_shared<Dumux::ParallelStokesOperator<Matrix, Vector, Vector, Communication>>(A, comms);
        auto op = std::make_shared<Dune::MatrixAdapter<Matrix, Vector, Vector>>(A);

        // auto pop = makeTpfaLaplaceOperator<typename Preconditioner::PressureLinearOperator>(pGridGeometry_->gridView(), *commPressure_);
        // auto pop2 = makePressureLinearOperator_<typename Preconditioner::PressureLinearOperator>(*commPressure_);
        auto pop = makePressureLinearOperator_<typename Preconditioner::PressureLinearOperator>();
        //auto preconditioner = std::make_shared<Preconditioner>(op, pop, pop2, comms, params_.sub("preconditioner"));
        auto preconditioner = std::make_shared<Preconditioner>(op, pop, params_.sub("preconditioner"));

        params_["verbose"] = pGridGeometry_->gridView().comm().rank() == 0 ? params_["verbose"] : "0";
        std::unique_ptr<Dune::InverseOperator<Vector, Vector>> solver;
        if (solverType_ == "minres")
            solver = std::make_unique<Dune::MINRESSolver<Vector>>(op, scalarProduct_, preconditioner, params_);
        else if (solverType_ == "bicgstab")
            solver = std::make_unique<Dune::BiCGSTABSolver<Vector>>(op, scalarProduct_, preconditioner, params_);
        else if (solverType_ == "cg")
            solver = std::make_unique<Dune::CGSolver<Vector>>(op, scalarProduct_, preconditioner, params_);
        else
            solver = std::make_unique<Dune::RestartedGMResSolver<Vector>>(op, scalarProduct_, preconditioner, params_);

        solver->apply(x, b, result_);

        return result_.converged;
    }

    template<class LinearOperator, class C>
    std::shared_ptr<LinearOperator> makePressureLinearOperator_(const C& comm)
    {
        using M = typename LinearOperator::matrix_type;
        auto massMatrix = createMassMatrix_<M>();
        return std::make_shared<LinearOperator>(massMatrix, comm);
    }

    template<class LinearOperator>
    std::shared_ptr<LinearOperator> makePressureLinearOperator_()
    {
        using M = typename LinearOperator::matrix_type;
        auto massMatrix = createMassMatrix_<M>();
        return std::make_shared<LinearOperator>(massMatrix);
    }

    template<class M>
    std::shared_ptr<M> createMassMatrix_()
    {
        auto massMatrix = std::make_shared<M>();
        massMatrix->setBuildMode(M::random);
        const auto numDofs = pGridGeometry_->numDofs();

        Dune::MatrixIndexSet pattern;
        pattern.resize(numDofs, numDofs);
        for (unsigned int globalI = 0; globalI < numDofs; ++globalI)
            pattern.add(globalI, globalI);
        pattern.exportIdx(*massMatrix);

        const auto& gv = pGridGeometry_->gridView();
        auto fvGeometry = localView(*pGridGeometry_);
        for (const auto& element : elements(gv))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                using Extrusion = Extrusion_t<PressureGG>;
                const auto dofIndex = scv.dofIndex();
                if (element.partitionType() == Dune::GhostEntity) // do not modify ghosts
                    (*massMatrix)[dofIndex][dofIndex] = 1.0;
                else
                    (*massMatrix)[dofIndex][dofIndex] += weight_*Extrusion::volume(fvGeometry, scv)/(2.0*viscosity_);
            }
        }

        return massMatrix;
    }

    void printOutMatrix_(const Matrix& A) const
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Matrix::N()>(), [&](const auto i) {
            forEach(std::make_index_sequence<Matrix::M()>(), [&](const auto j) {
                if (A[i][j].nonzeroes() > 0)
                    Dune::printmatrix(std::cout, A[i][j], "A" + std::to_string(i()) + std::to_string(j()), "");
                else
                    std::cout << "A" << i() << j() << ": 0" << std::endl;
            });
        });
    }

    void printOutResidual_(const Vector& b) const
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Vector::size()>(), [&](const auto i) {
            Dune::printvector(std::cout, b[i], "b" + std::to_string(i()), "");
        });
    }

    double density_, viscosity_, weight_;
    Dune::InverseOperatorResult result_;
    Dune::ParameterTree params_;
    std::shared_ptr<const VelocityGG> vGridGeometry_;
    std::shared_ptr<const PressureGG> pGridGeometry_;
    const Vector& dirichletDofs_;
    std::string solverType_;

// #if HAVE_MPI
//     std::unique_ptr<ParallelISTLHelper<LinearSolverTraits<VelocityGG>>> pHelperVelocity_;
//     std::unique_ptr<ParallelISTLHelper<LinearSolverTraits<PressureGG>>> pHelperPressure_;

//     std::shared_ptr<Communication> commVelocity_;
//     std::shared_ptr<Communication> commPressure_;
// #endif
    std::shared_ptr<Dune::ScalarProduct<Vector>> scalarProduct_;
};

} // end namespace Dumux

#endif
